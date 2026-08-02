#include <cstddef>

#include <cuda_runtime.h>

#include "cutools.h"
#include "errr.h"
#include "f_arr_view.h"
#include <iostream>
namespace mglet::backend
{

namespace
{

__global__ void maxabscal_kernel(
    mgletreal* __restrict__ maxabsgrid,
    const mgletreal* __restrict__ phi,
    const mgletint* __restrict__ mygrids,
    int nmygrids,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygrids)
    {
        return;
    }

    const auto igrid = mygrids[block_idx] - 1;

    const auto kk = kkk[igrid];
    const auto jj = jjj[igrid];
    const auto ii = iii[igrid];
    const auto ip3 = ip3d[igrid] - 1;

    auto local_max = mgletreal{0.0};

    for (int i = 2; i <= ii - 3; ++i)
    {
        for (int j = 2 + threadIdx.y; j <= jj - 3; j += blockDim.y)
        {
            for (int k = 2 + threadIdx.x; k <= kk - 3; k += blockDim.x)
            {
                const auto idx = ip3 + k + (std::int64_t)j * kk + (std::int64_t)i * kk * jj;
                const auto v = fabs(phi[idx]);
                local_max = v > local_max ? v : local_max;
            }
        }
    }

    extern __shared__ mgletreal sdata[];
    const auto tid = threadIdx.y * blockDim.x + threadIdx.x;
    sdata[tid] = local_max;
    __syncthreads();

    for (int s = (blockDim.x * blockDim.y) / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] = sdata[tid] > sdata[tid + s] ? sdata[tid] : sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0)
    {
        maxabsgrid[block_idx] = sdata[0];
    }
}

__global__ void rescal_kernel(
    mgletreal* __restrict__ rhs,
    const mgletreal* __restrict__ res,
    mgletint nmygrids,
    const mgletint* __restrict__ mygrids,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygrids)
    {
        return;
    }

    const auto igrid = mygrids[block_idx] - 1; // C is 0-based

    const auto kk = kkk[igrid];
    const auto jj = jjj[igrid];
    const auto ii = iii[igrid];
    const auto ip3 = ip3d[igrid] - 1; // C is 0-based

    for (mgletint i = 2 + blockIdx.y; i <= ii - 3; i += gridDim.y)
    {
        for (mgletint j = 2 + threadIdx.y; j <= jj - 3; j += blockDim.y)
        {
            for (mgletint k = 2 + threadIdx.x; k <= kk - 3; k += blockDim.x)
            {
                const auto idx = ip3 + k + j * kk + i * kk * jj;
                rhs[idx] = rhs[idx] + res[idx];
            }
        }
    }
}

__global__ void sipiter1_hyperplane_level_kernel(
    mgletreal* __restrict__ res,
    const mgletreal* __restrict__ rhs,
    const mgletreal* __restrict__ lw,
    const mgletreal* __restrict__ ls,
    const mgletreal* __restrict__ lb,
    const mgletreal* __restrict__ lpr,
    const mgletifk* __restrict__ mip,
    const mgletifk* __restrict__ idxsip,
    mgletint nmygridsonlvl,
    const mgletint* __restrict__ mygridsonlvl,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygridsonlvl)
    {
        return;
    }

    const auto igrid = mygridsonlvl[block_idx] - 1;

    const auto kk = kkk[igrid];
    const auto jj = jjj[igrid];
    const auto ii = iii[igrid];
    const auto ip3 = ip3d[igrid] - 1;

    const auto n3dmin = 3 + 3 + 3;
    const auto n3dmax = (ii - 2) + (jj - 2) + (kk - 2);

    __shared__ mgletifk s_lm, s_lp;

    for (mgletint m = n3dmin; m <= n3dmax; ++m)
    {

        if (threadIdx.x == 0)
        {
            const mgletint lm = static_cast<mgletint>(mip[ip3 + m - 1]);
            const mgletint lm1 = static_cast<mgletint>(mip[ip3 + m]);
            s_lm = lm;
            s_lp = lm1 - lm;
        }
        __syncthreads();

        const mgletint lm = s_lm;
        const mgletint lp = s_lp;

        for (mgletint ipp = 1 + threadIdx.x; ipp <= lp; ipp += blockDim.x)
        {
            const mgletint iacc = lm + ipp;

            const mgletint sip_idx = ip3 + iacc - 1;

            const mgletint local_idx = static_cast<mgletint>(idxsip[sip_idx]);
            const mgletint idx = ip3 + local_idx - 1;
            const mgletint idx_km = idx - 1;
            const mgletint idx_jm = idx - kk;
            const mgletint idx_im = idx - kk * jj;

            mgletreal val = (rhs[idx] + res[idx]) * lpr[sip_idx];
            val -= lb[sip_idx] * res[idx_km];
            val -= ls[sip_idx] * res[idx_jm];
            val -= lw[sip_idx] * res[idx_im];

            res[idx] = val;
        }

        __syncthreads();
    }
}

__global__ void sipiter2_hyperplane_level_kernel(
    mgletreal* __restrict__ dp,
    mgletreal* __restrict__ res,
    const mgletreal* __restrict__ ue,
    const mgletreal* __restrict__ un,
    const mgletreal* __restrict__ ut,
    const mgletifk* __restrict__ mip,
    const mgletifk* __restrict__ idxsip,
    const mgletint* __restrict__ mygridsonlvl,
    mgletint nmygridsonlvl,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygridsonlvl) return;

    const auto igrid = mygridsonlvl[block_idx] - 1;

    const auto kk = kkk[igrid];
    const auto jj = jjj[igrid];
    const auto ii = iii[igrid];
    const auto ip3 = ip3d[igrid] - 1;

    const mgletint n3dmin = 3 + 3 + 3;
    const mgletint n3dmax = (ii - 2) + (jj - 2) + (kk - 2);

    __shared__ mgletifk s_lm, s_lp;

    for (mgletint m = n3dmax; m >= n3dmin; --m) {

        if (threadIdx.x == 0) {
            const mgletint lm  = static_cast<mgletint>(mip[ip3 + m - 1]);
            const mgletint lm1 = static_cast<mgletint>(mip[ip3 + m]);
            s_lm = lm;
            s_lp = lm1 - lm;
        }
        __syncthreads();

        const mgletint lm = s_lm;
        const mgletint lp = s_lp;

        for (mgletint ipp = 1 + threadIdx.x; ipp <= lp; ipp += blockDim.x) {
            const mgletint iacc = lm + ipp;
            const mgletint sip_idx = ip3 + iacc - 1;

            const mgletint local_idx = static_cast<mgletint>(idxsip[sip_idx]);
            const mgletint idx    = ip3 + local_idx - 1;
            const mgletint idx_kp = idx + 1;
            const mgletint idx_jp = idx + kk;
            const mgletint idx_ip = idx + kk * jj;

            mgletreal val = res[idx];
            val -= ut[sip_idx] * res[idx_kp];
            val -= un[sip_idx] * res[idx_jp];
            val -= ue[sip_idx] * res[idx_ip];

            res[idx] = val;
        }

        __syncthreads();
    }

    const auto ni = ii - 4;
    const auto nj = jj - 4;
    const auto nk = kk - 4;
    if (ni > 0 && nj > 0 && nk > 0) {
        const auto n = ni * nj * nk;
        for (mgletint lin = threadIdx.x; lin < n; lin += blockDim.x) {
            const int k = 2 + lin % nk;
            const int j = 2 + (lin / nk) % nj;
            const int i = 2 + lin / ((std::int64_t)nk * nj);

            const auto idx = ip3 + k + j * kk + i * kk * jj;
            dp[idx] = dp[idx] + res[idx];
        }
    }
}

} // namespace

void maxabscal_backend(
    FArrView<mgletreal> maxabsgrid,
    FArrView<const mgletreal> phi,
    FArrView<const mgletint> mygrids,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d)
{
    const auto nmygrids = mygrids.flat_size();

    if (nmygrids == 0)
    {
        return;
    }

    const auto threads = ::dim3{32, 8};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygrids)};
    const std::size_t shmem_bytes = threads.x * threads.y * sizeof(mgletreal);

    maxabscal_kernel<<<blocks, threads, shmem_bytes>>>(
        maxabsgrid.device_ptr(),
        phi.device_ptr(),
        mygrids.device_ptr(),
        nmygrids,
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void rescal_backend(
    FArrView<mgletreal> rhs,
    FArrView<const mgletreal> res,
    FArrView<const mgletint> mygrids,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d)
{
    const auto nmygrids = mygrids.flat_size();

    const auto n_rhs = rhs.flat_size();
    const auto n_res = res.flat_size();

    if (n_rhs != n_res)
    {
        MGLET_ERRR();
    }

    if (n_rhs == 0)
    {
        return;
    }

    const auto threads = ::dim3{32, 8};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygrids)};

    rescal_kernel<<<blocks, threads>>>(
        rhs.device_ptr(),
        res.device_ptr(),
        nmygrids,
        mygrids.device_ptr(),
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void sipiter1_hyperplane_level_backend(
    FArrView<mgletreal>(res),
    FArrView<const mgletreal> rhs,
    FArrView<const mgletreal> siplw,
    FArrView<const mgletreal> sipls,
    FArrView<const mgletreal> siplb,
    FArrView<const mgletreal> siplpr,
    FArrView<const mgletifk> miphp,
    FArrView<const mgletifk> idxhp,
    FArrView<const mgletint> mygridsonlvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d)
{
    const auto nmygridsonlvl = mygridsonlvl.flat_size();

    const auto threads = ::dim3{256};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygridsonlvl)};

    sipiter1_hyperplane_level_kernel<<<blocks, threads>>>(
        res.device_ptr(),
        rhs.device_ptr(),
        siplw.device_ptr(),
        sipls.device_ptr(),
        siplb.device_ptr(),
        siplpr.device_ptr(),
        miphp.device_ptr(),
        idxhp.device_ptr(),
        nmygridsonlvl,
        mygridsonlvl.device_ptr(),
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void sipiter2_hyperplane_level_backend(
    FArrView<mgletreal> dp,
    FArrView<mgletreal> res,
    FArrView<const mgletreal> sipue,
    FArrView<const mgletreal> sipun,
    FArrView<const mgletreal> siput,
    FArrView<const mgletifk> miphp,
    FArrView<const mgletifk> idxhp,
    FArrView<const mgletint> mygridsonlvl,
    FArrView<const mgletint> kkk,
    FArrView<const mgletint> jjj,
    FArrView<const mgletint> iii,
    FArrView<const mgletint> ip3d)
{
    const auto nmygridsonlvl = mygridsonlvl.flat_size();

    const auto threads = ::dim3{256};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygridsonlvl)};

    sipiter2_hyperplane_level_kernel<<<blocks, threads>>>(
        dp.device_ptr(),
        res.device_ptr(),
        sipue.device_ptr(),
        sipun.device_ptr(),
        siput.device_ptr(),
        miphp.device_ptr(),
        idxhp.device_ptr(),
        mygridsonlvl.device_ptr(),
        nmygridsonlvl,
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}


} // namespace mglet::backend
