#include "bound_pressure_backend.h"

#include <cstddef>
#include <cstdint>

#include <cuda_runtime.h>

#include "cutools.h"
#include "errr.h"
#include "f_arr_view.h"

namespace mglet::backend
{

namespace
{

// ---------------------------------------------------------------------
// Address helpers. All offsets (ip3, ipx, ipy, ipz, ipbb) are assumed
// already 0-based (== Fortran get_ip3/get_ip1x/.../get_ipbb value - 1).
// (k,j,i) arguments here are Fortran 1-based, matching the original code.
// ---------------------------------------------------------------------

__device__ __forceinline__ std::int64_t idx3(std::int64_t base, int k, int j, int i, int kk, int jj)
{
    return base + (k - 1) + (std::int64_t)(j - 1) * kk + (std::int64_t)(i - 1) * kk * jj;
}

__device__ __forceinline__ std::int64_t idx1(std::int64_t base, int i)
{
    return base + (i - 1);
}

// pbuffer(kk, jj) layout -- used by bfront/bfront_bp
__device__ __forceinline__ std::int64_t idxbuf_kj(std::int64_t base, int k, int j, int kk)
{
    return base + (k - 1) + (std::int64_t)(j - 1) * kk;
}

// pbuffer(kk, ii) layout -- used by bright/bright_bp
__device__ __forceinline__ std::int64_t idxbuf_ki(std::int64_t base, int k, int i, int kk)
{
    return base + (k - 1) + (std::int64_t)(i - 1) * kk;
}

// pbuffer(jj, ii) layout -- used by bbottom/bbottom_bp
__device__ __forceinline__ std::int64_t idxbuf_ji(std::int64_t base, int j, int i, int jj)
{
    return base + (j - 1) + (std::int64_t)(i - 1) * jj;
}

// Trip count for a Fortran "DO x = start, stop, 2" loop.
__device__ __forceinline__ int step2_count(int start, int stop)
{
    if (stop < start)
        return 0;
    return (stop - start) / 2 + 1;
}

// ---------------------------------------------------------------------
// pressureftocone_A / _B
// ---------------------------------------------------------------------

__device__ __forceinline__ mgletreal pressureftocone_A(
    const mgletreal* __restrict__ p,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    int kk,
    int jj,
    int k,
    int j,
    int i)
{
    mgletreal sump = 0.0, sumvol = 0.0;
    for (int l = 0; l <= 1; ++l)
    {
        for (int m = 0; m <= 1; ++m)
        {
            for (int n = 0; n <= 1; ++n)
            {
                const mgletreal vol =
                    ddz[idx1(ipz, k + n)] * ddy[idx1(ipy, j + m)] * ddx[idx1(ipx, i + l)];
                sump += p[idx3(ip3, k + n, j + m, i + l, kk, jj)] * vol;
                sumvol += vol;
            }
        }
    }
    return sump / sumvol;
}

// Returns pc; writes bpc via pointer (matches Fortran's dual-output signature).
__device__ __forceinline__ mgletreal pressureftocone_B(
    const mgletreal* __restrict__ p,
    const mgletreal* __restrict__ bp,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    int kk,
    int jj,
    int k,
    int j,
    int i,
    mgletreal* bpc_out)
{
    mgletreal sumbp = 0.0;
    for (int l = 0; l <= 1; ++l)
        for (int m = 0; m <= 1; ++m)
            for (int n = 0; n <= 1; ++n)
                sumbp += bp[idx3(ip3, k + n, j + m, i + l, kk, jj)];

    const mgletreal bpc = sumbp < mgletreal(1.0) ? sumbp : mgletreal(1.0);
    *bpc_out = bpc;

    if (bpc < mgletreal(0.5))
    {
        return mgletreal(0.0);
    }

    mgletreal sump = 0.0, sumvol = 0.0;
    for (int l = 0; l <= 1; ++l)
    {
        for (int m = 0; m <= 1; ++m)
        {
            for (int n = 0; n <= 1; ++n)
            {
                const mgletreal vol = bp[idx3(ip3, k + n, j + m, i + l, kk, jj)] *
                                      ddz[idx1(ipz, k + n)] * ddy[idx1(ipy, j + m)] *
                                      ddx[idx1(ipx, i + l)];
                sump += p[idx3(ip3, k + n, j + m, i + l, kk, jj)] * vol;
                sumvol += vol;
            }
        }
    }
    return sump / sumvol;
}

// ---------------------------------------------------------------------
// Non-BP face kernels (bfront / bright / bbottom)
// ---------------------------------------------------------------------

__device__ void bfront_task(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletreal* __restrict__ dx,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    std::int64_t ipbb,
    int kk,
    int jj,
    int i2,
    int i3,
    int i4,
    int istag2)
{
    const int i = i3 < i4 ? i3 : i4;
    const int nj = step2_count(3, jj - 2);
    const int nk = step2_count(3, kk - 2);
    if (nj <= 0 || nk <= 0)
        return;

    const std::int64_t n = (std::int64_t)nj * nk;
    for (std::int64_t lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const int hj = (int)(lin / nk);
        const int hk = (int)(lin % nk);
        const int j = 3 + 2 * hj;
        const int k = 3 + 2 * hk;

        const mgletreal pcnew =
            pressureftocone_A(p, ddx, ddy, ddz, ip3, ipx, ipy, ipz, kk, jj, k, j, i);
        const mgletreal delta = dx[idx1(ipx, istag2)] / (ddx[idx1(ipx, i3)] + ddx[idx1(ipx, i2)]) *
                                (pbuffer[idxbuf_kj(ipbb, k, j, kk)] - pcnew);

        for (int m = 0; m <= 1; ++m)
        {
            for (int nn = 0; nn <= 1; ++nn)
            {
                p[idx3(ip3, k + nn, j + m, i2, kk, jj)] =
                    p[idx3(ip3, k + nn, j + m, i3, kk, jj)] + delta;
            }
        }
    }
}

__device__ void bright_task(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletreal* __restrict__ dy,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    std::int64_t ipbb,
    int kk,
    int jj,
    int ii,
    int j2,
    int j3,
    int j4,
    int jstag2)
{
    const int j = j3 < j4 ? j3 : j4;
    const int ni = step2_count(3, ii - 2);
    const int nk = step2_count(3, kk - 2);
    if (ni <= 0 || nk <= 0)
        return;

    const std::int64_t n = (std::int64_t)ni * nk;
    for (std::int64_t lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const int hi = (int)(lin / nk);
        const int hk = (int)(lin % nk);
        const int i = 3 + 2 * hi;
        const int k = 3 + 2 * hk;

        const mgletreal pcnew =
            pressureftocone_A(p, ddx, ddy, ddz, ip3, ipx, ipy, ipz, kk, jj, k, j, i);
        const mgletreal delta = dy[idx1(ipy, jstag2)] / (ddy[idx1(ipy, j3)] + ddy[idx1(ipy, j2)]) *
                                (pbuffer[idxbuf_ki(ipbb, k, i, kk)] - pcnew);

        for (int l = 0; l <= 1; ++l)
        {
            for (int nn = 0; nn <= 1; ++nn)
            {
                p[idx3(ip3, k + nn, j2, i + l, kk, jj)] =
                    p[idx3(ip3, k + nn, j3, i + l, kk, jj)] + delta;
            }
        }
    }
}

__device__ void bbottom_task(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletreal* __restrict__ dz,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    std::int64_t ipbb,
    int kk,
    int jj,
    int ii,
    int k2,
    int k3,
    int k4,
    int kstag2)
{
    const int k = k3 < k4 ? k3 : k4;
    const int ni = step2_count(3, ii - 2);
    const int nj = step2_count(3, jj - 2);
    if (ni <= 0 || nj <= 0)
        return;

    const std::int64_t n = (std::int64_t)ni * nj;
    for (std::int64_t lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const int hi = (int)(lin / nj);
        const int hj = (int)(lin % nj);
        const int i = 3 + 2 * hi;
        const int j = 3 + 2 * hj;

        const mgletreal pcnew =
            pressureftocone_A(p, ddx, ddy, ddz, ip3, ipx, ipy, ipz, kk, jj, k, j, i);
        const mgletreal delta = dz[idx1(ipz, kstag2)] / (ddz[idx1(ipz, k3)] + ddz[idx1(ipz, k2)]) *
                                (pbuffer[idxbuf_ji(ipbb, j, i, jj)] - pcnew);

        for (int l = 0; l <= 1; ++l)
        {
            for (int m = 0; m <= 1; ++m)
            {
                p[idx3(ip3, k2, j + m, i + l, kk, jj)] =
                    p[idx3(ip3, k3, j + m, i + l, kk, jj)] + delta;
            }
        }
    }
}

// ---------------------------------------------------------------------
// BP face kernels (bfront_bp / bright_bp / bbottom_bp)
// ---------------------------------------------------------------------

__device__ void bfront_bp_task(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ bp,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletreal* __restrict__ dx,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    std::int64_t ipbb,
    int kk,
    int jj,
    int i2,
    int i3,
    int i4,
    int istag2)
{
    const int i = i3 < i4 ? i3 : i4;
    const int nj = step2_count(3, jj - 2);
    const int nk = step2_count(3, kk - 2);
    if (nj <= 0 || nk <= 0)
        return;

    const std::int64_t n = (std::int64_t)nj * nk;
    for (std::int64_t lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const int hj = (int)(lin / nk);
        const int hk = (int)(lin % nk);
        const int j = 3 + 2 * hj;
        const int k = 3 + 2 * hk;

        mgletreal bpc;
        mgletreal pcnew =
            pressureftocone_B(p, bp, ddx, ddy, ddz, ip3, ipx, ipy, ipz, kk, jj, k, j, i, &bpc);
        if (bpc < mgletreal(0.5))
        {
            pcnew = pbuffer[idxbuf_kj(ipbb, k, j, kk)];
        }

        const mgletreal sb11 = bp[idx3(ip3, k, j, i2, kk, jj)] * bp[idx3(ip3, k, j, i3, kk, jj)];
        const mgletreal sb12 =
            bp[idx3(ip3, k, j + 1, i2, kk, jj)] * bp[idx3(ip3, k, j + 1, i3, kk, jj)];
        const mgletreal sb13 =
            bp[idx3(ip3, k + 1, j, i2, kk, jj)] * bp[idx3(ip3, k + 1, j, i3, kk, jj)];
        const mgletreal sb14 =
            bp[idx3(ip3, k + 1, j + 1, i2, kk, jj)] * bp[idx3(ip3, k + 1, j + 1, i3, kk, jj)];

        mgletreal fak = (sb11 * ddy[idx1(ipy, j)] * ddz[idx1(ipz, k)] +
                         sb12 * ddy[idx1(ipy, j + 1)] * ddz[idx1(ipz, k)] +
                         sb13 * ddy[idx1(ipy, j)] * ddz[idx1(ipz, k + 1)] +
                         sb14 * ddy[idx1(ipy, j + 1)] * ddz[idx1(ipz, k + 1)]) /
                        ((ddy[idx1(ipy, j)] + ddy[idx1(ipy, j + 1)]) *
                         (ddz[idx1(ipz, k)] + ddz[idx1(ipz, k + 1)]));
        if (fak < mgletreal(0.1))
            fak = mgletreal(1.0);
        fak = mgletreal(1.0) / fak;

        const mgletreal delta = fak * dx[idx1(ipx, istag2)] /
                                (ddx[idx1(ipx, i3)] + ddx[idx1(ipx, i2)]) *
                                (pbuffer[idxbuf_kj(ipbb, k, j, kk)] - pcnew);

        for (int m = 0; m <= 1; ++m)
        {
            for (int nn = 0; nn <= 1; ++nn)
            {
                p[idx3(ip3, k + nn, j + m, i2, kk, jj)] =
                    p[idx3(ip3, k + nn, j + m, i3, kk, jj)] + delta;
            }
        }
    }
}

__device__ void bright_bp_task(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ bp,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletreal* __restrict__ dy,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    std::int64_t ipbb,
    int kk,
    int jj,
    int ii,
    int j2,
    int j3,
    int j4,
    int jstag2)
{
    const int j = j3 < j4 ? j3 : j4;
    const int ni = step2_count(3, ii - 2);
    const int nk = step2_count(3, kk - 2);
    if (ni <= 0 || nk <= 0)
        return;

    const std::int64_t n = (std::int64_t)ni * nk;
    for (std::int64_t lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const int hi = (int)(lin / nk);
        const int hk = (int)(lin % nk);
        const int i = 3 + 2 * hi;
        const int k = 3 + 2 * hk;

        mgletreal bpc;
        mgletreal pcnew =
            pressureftocone_B(p, bp, ddx, ddy, ddz, ip3, ipx, ipy, ipz, kk, jj, k, j, i, &bpc);
        if (bpc < mgletreal(0.5))
        {
            pcnew = pbuffer[idxbuf_ki(ipbb, k, i, kk)];
        }

        const mgletreal sb11 = bp[idx3(ip3, k, j2, i, kk, jj)] * bp[idx3(ip3, k, j3, i, kk, jj)];
        const mgletreal sb12 =
            bp[idx3(ip3, k, j2, i + 1, kk, jj)] * bp[idx3(ip3, k, j3, i + 1, kk, jj)];
        const mgletreal sb13 =
            bp[idx3(ip3, k + 1, j2, i, kk, jj)] * bp[idx3(ip3, k + 1, j3, i, kk, jj)];
        const mgletreal sb14 =
            bp[idx3(ip3, k + 1, j2, i + 1, kk, jj)] * bp[idx3(ip3, k + 1, j3, i + 1, kk, jj)];

        mgletreal fak = (sb11 * ddx[idx1(ipx, i)] * ddz[idx1(ipz, k)] +
                         sb12 * ddx[idx1(ipx, i + 1)] * ddz[idx1(ipz, k)] +
                         sb13 * ddx[idx1(ipx, i)] * ddz[idx1(ipz, k + 1)] +
                         sb14 * ddx[idx1(ipx, i + 1)] * ddz[idx1(ipz, k + 1)]) /
                        ((ddx[idx1(ipx, i)] + ddx[idx1(ipx, i + 1)]) *
                         (ddz[idx1(ipz, k)] + ddz[idx1(ipz, k + 1)]));
        if (fak < mgletreal(0.1))
            fak = mgletreal(1.0);
        fak = mgletreal(1.0) / fak;

        const mgletreal delta = fak * dy[idx1(ipy, jstag2)] /
                                (ddy[idx1(ipy, j3)] + ddy[idx1(ipy, j2)]) *
                                (pbuffer[idxbuf_ki(ipbb, k, i, kk)] - pcnew);

        for (int l = 0; l <= 1; ++l)
        {
            for (int nn = 0; nn <= 1; ++nn)
            {
                p[idx3(ip3, k + nn, j2, i + l, kk, jj)] =
                    p[idx3(ip3, k + nn, j3, i + l, kk, jj)] + delta;
            }
        }
    }
}

__device__ void bbottom_bp_task(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ bp,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletreal* __restrict__ dz,
    std::int64_t ip3,
    std::int64_t ipx,
    std::int64_t ipy,
    std::int64_t ipz,
    std::int64_t ipbb,
    int kk,
    int jj,
    int ii,
    int k2,
    int k3,
    int k4,
    int kstag2)
{
    const int k = k3 < k4 ? k3 : k4;
    const int ni = step2_count(3, ii - 2);
    const int nj = step2_count(3, jj - 2);
    if (ni <= 0 || nj <= 0)
        return;

    const std::int64_t n = (std::int64_t)ni * nj;
    for (std::int64_t lin = threadIdx.x; lin < n; lin += blockDim.x)
    {
        const int hi = (int)(lin / nj);
        const int hj = (int)(lin % nj);
        const int i = 3 + 2 * hi;
        const int j = 3 + 2 * hj;

        mgletreal bpc;
        mgletreal pcnew =
            pressureftocone_B(p, bp, ddx, ddy, ddz, ip3, ipx, ipy, ipz, kk, jj, k, j, i, &bpc);
        if (bpc < mgletreal(0.5))
        {
            pcnew = pbuffer[idxbuf_ji(ipbb, j, i, jj)];
        }

        const mgletreal sb11 = bp[idx3(ip3, k2, j, i, kk, jj)] * bp[idx3(ip3, k3, j, i, kk, jj)];
        const mgletreal sb12 =
            bp[idx3(ip3, k2, j, i + 1, kk, jj)] * bp[idx3(ip3, k3, j, i + 1, kk, jj)];
        const mgletreal sb13 =
            bp[idx3(ip3, k2, j + 1, i, kk, jj)] * bp[idx3(ip3, k3, j + 1, i, kk, jj)];
        const mgletreal sb14 =
            bp[idx3(ip3, k2, j + 1, i + 1, kk, jj)] * bp[idx3(ip3, k3, j + 1, i + 1, kk, jj)];

        mgletreal fak = (sb11 * ddx[idx1(ipx, i)] * ddy[idx1(ipy, j)] +
                         sb12 * ddx[idx1(ipx, i + 1)] * ddy[idx1(ipy, j)] +
                         sb13 * ddx[idx1(ipx, i)] * ddy[idx1(ipy, j + 1)] +
                         sb14 * ddx[idx1(ipx, i + 1)] * ddy[idx1(ipy, j + 1)]) /
                        ((ddx[idx1(ipx, i)] + ddx[idx1(ipx, i + 1)]) *
                         (ddy[idx1(ipy, j)] + ddy[idx1(ipy, j + 1)]));
        if (fak < mgletreal(0.1))
            fak = mgletreal(1.0);
        fak = mgletreal(1.0) / fak;

        const mgletreal delta = fak * dz[idx1(ipz, kstag2)] /
                                (ddz[idx1(ipz, k3)] + ddz[idx1(ipz, k2)]) *
                                (pbuffer[idxbuf_ji(ipbb, j, i, jj)] - pcnew);

        for (int l = 0; l <= 1; ++l)
        {
            for (int m = 0; m <= 1; ++m)
            {
                p[idx3(ip3, k2, j + m, i + l, kk, jj)] =
                    p[idx3(ip3, k3, j + m, i + l, kk, jj)] + delta;
            }
        }
    }
}

__global__ void bound_pressure_bp_kernel(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ bp,
    const mgletreal* __restrict__ dx,
    const mgletreal* __restrict__ dy,
    const mgletreal* __restrict__ dz,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletint* __restrict__ boundtasks_lvl,
    mgletint nboundtasks,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d,
    const mgletint* __restrict__ ip1dx,
    const mgletint* __restrict__ ip1dy,
    const mgletint* __restrict__ ip1dz,
    const mgletint* __restrict__ ipbb_arr)
{
    const auto itask = blockIdx.x;
    if (itask >= nboundtasks)
    {
        return;
    }

    // Column major nd Fortran array access...
    const auto igrid = boundtasks_lvl[2 * itask + 0] - 1;
    const auto iface = boundtasks_lvl[2 * itask + 1];

    if (igrid < 0)
    {
        // dummy/unused slow
        return;
    }

    const auto kk = kkk[igrid];
    const auto jj = jjj[igrid];
    const auto ii = iii[igrid];
    const auto ip3 = ip3d[igrid] - 1;
    const auto ipx = ip1dx[igrid] - 1;
    const auto ipy = ip1dy[igrid] - 1;
    const auto ipz = ip1dz[igrid] - 1;

    // Column major nd Fortran array access...
    const auto ipbb = ipbb_arr[(iface - 1) + igrid * 6] - 1;

    switch (iface)
    {
        case 1:
            bfront_bp_task(
                p, pbuffer, bp, ddx, ddy, ddz, dx, ip3, ipx, ipy, ipz, ipbb, kk, jj, 2, 3, 4, 2);
            break;
        case 2:
            bfront_bp_task(
                p,
                pbuffer,
                bp,
                ddx,
                ddy,
                ddz,
                dx,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii - 1,
                ii - 2,
                ii - 3,
                ii - 2);
            break;
        case 3:
            bright_bp_task(
                p,
                pbuffer,
                bp,
                ddx,
                ddy,
                ddz,
                dy,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii,
                2,
                3,
                4,
                2);
            break;
        case 4:
            bright_bp_task(
                p,
                pbuffer,
                bp,
                ddx,
                ddy,
                ddz,
                dy,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii,
                jj - 1,
                jj - 2,
                jj - 3,
                jj - 2);
            break;
        case 5:
            bbottom_bp_task(
                p,
                pbuffer,
                bp,
                ddx,
                ddy,
                ddz,
                dz,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii,
                2,
                3,
                4,
                2);
            break;
        case 6:
            bbottom_bp_task(
                p,
                pbuffer,
                bp,
                ddx,
                ddy,
                ddz,
                dz,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii,
                kk - 1,
                kk - 2,
                kk - 3,
                kk - 2);
            break;
        default: return;
    }
}

__global__ void bound_pressure_nobp_kernel(
    mgletreal* __restrict__ p,
    const mgletreal* __restrict__ pbuffer,
    const mgletreal* __restrict__ dx,
    const mgletreal* __restrict__ dy,
    const mgletreal* __restrict__ dz,
    const mgletreal* __restrict__ ddx,
    const mgletreal* __restrict__ ddy,
    const mgletreal* __restrict__ ddz,
    const mgletint* __restrict__ boundtasks_lvl,
    int nboundtasks,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d,
    const mgletint* __restrict__ ip1dx,
    const mgletint* __restrict__ ip1dy,
    const mgletint* __restrict__ ip1dz,
    const mgletint* __restrict__ ipbb_arr)
{
    const int itask = blockIdx.x;
    if (itask >= nboundtasks)
    {
        return;
    }

    // Column major nd Fortran array access...
    const auto igrid = boundtasks_lvl[2 * itask + 0] - 1;
    const auto iface = boundtasks_lvl[2 * itask + 1];

    if (igrid < 0)
    {
        return;
    }

    const auto kk = kkk[igrid];
    const auto jj = jjj[igrid];
    const auto ii = iii[igrid];
    const auto ip3 = ip3d[igrid] - 1;
    const auto ipx = ip1dx[igrid] - 1;
    const auto ipy = ip1dy[igrid] - 1;
    const auto ipz = ip1dz[igrid] - 1;
    const auto ipbb = ipbb_arr[(iface - 1) + igrid * 6] - 1;

    switch (iface)
    {
        case 1:
            bfront_task(
                p, pbuffer, ddx, ddy, ddz, dx, ip3, ipx, ipy, ipz, ipbb, kk, jj, 2, 3, 4, 2);
            break;
        case 2:
            bfront_task(
                p,
                pbuffer,
                ddx,
                ddy,
                ddz,
                dx,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii - 1,
                ii - 2,
                ii - 3,
                ii - 2);
            break;
        case 3:
            bright_task(
                p, pbuffer, ddx, ddy, ddz, dy, ip3, ipx, ipy, ipz, ipbb, kk, jj, ii, 2, 3, 4, 2);
            break;
        case 4:
            bright_task(
                p,
                pbuffer,
                ddx,
                ddy,
                ddz,
                dy,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii,
                jj - 1,
                jj - 2,
                jj - 3,
                jj - 2);
            break;
        case 5:
            bbottom_task(
                p, pbuffer, ddx, ddy, ddz, dz, ip3, ipx, ipy, ipz, ipbb, kk, jj, ii, 2, 3, 4, 2);
            break;
        case 6:
            bbottom_task(
                p,
                pbuffer,
                ddx,
                ddy,
                ddz,
                dz,
                ip3,
                ipx,
                ipy,
                ipz,
                ipbb,
                kk,
                jj,
                ii,
                kk - 1,
                kk - 2,
                kk - 3,
                kk - 2);
            break;
        default: return;
    }
}

} // namespace

void bound_pressure_bp_backend(
    FArrView<mgletreal> p,
    FArrView<const mgletreal> pbuffer,
    FArrView<const mgletreal> bp,
    FArrView<const mgletreal> dx,
    FArrView<const mgletreal> dy,
    FArrView<const mgletreal> dz,
    FArrView<const mgletreal> ddx,
    FArrView<const mgletreal> ddy,
    FArrView<const mgletreal> ddz,
    mgletint nboundtasks,
    FArrView<const mgletint> boundtasks_lvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d,
    FArrView<const mgletint> ip1dx,
    FArrView<const mgletint> ip1dy,
    FArrView<const mgletint> ip1dz,
    FArrView<const mgletint> ipbb)
{
    if (nboundtasks == 0)
    {
        return;
    }

    const auto threads = ::dim3{256};
    const auto blocks = ::dim3{static_cast<unsigned>(nboundtasks)};

    bound_pressure_bp_kernel<<<blocks, threads>>>(
        p.device_ptr(),
        pbuffer.device_ptr(),
        bp.device_ptr(),
        dx.device_ptr(),
        dy.device_ptr(),
        dz.device_ptr(),
        ddx.device_ptr(),
        ddy.device_ptr(),
        ddz.device_ptr(),
        boundtasks_lvl.device_ptr(),
        nboundtasks,
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr(),
        ip1dx.device_ptr(),
        ip1dy.device_ptr(),
        ip1dz.device_ptr(),
        ipbb.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void bound_pressure_nobp_backend(
    FArrView<mgletreal> p,
    FArrView<const mgletreal> pbuffer,
    FArrView<const mgletreal> dx,
    FArrView<const mgletreal> dy,
    FArrView<const mgletreal> dz,
    FArrView<const mgletreal> ddx,
    FArrView<const mgletreal> ddy,
    FArrView<const mgletreal> ddz,
    mgletint nboundtasks,
    FArrView<const mgletint> boundtasks_lvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d,
    FArrView<const mgletint> ip1dx,
    FArrView<const mgletint> ip1dy,
    FArrView<const mgletint> ip1dz,
    FArrView<const mgletint> ipbb)
{
    if (nboundtasks == 0)
    {
        return;
    }

    const auto threads = ::dim3{256};
    const auto blocks = ::dim3{static_cast<unsigned>(nboundtasks)};

    bound_pressure_nobp_kernel<<<blocks, threads>>>(
        p.device_ptr(),
        pbuffer.device_ptr(),
        dx.device_ptr(),
        dy.device_ptr(),
        dz.device_ptr(),
        ddx.device_ptr(),
        ddy.device_ptr(),
        ddz.device_ptr(),
        boundtasks_lvl.device_ptr(),
        nboundtasks,
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr(),
        ip1dx.device_ptr(),
        ip1dy.device_ptr(),
        ip1dz.device_ptr(),
        ipbb.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

} // namespace mglet::backend
