#include "fieldhelper_backend.h"

#include <cstddef>

#include <cuda_runtime.h>

#include "cutools.h"
#include "errr.h"
#include "f_arr_view.h"

namespace mglet::backend
{

namespace
{

__global__ void laplacephi_kernel(
    mgletreal* __restrict__ res,
    const mgletreal* __restrict__ phi,
    const mgletreal* __restrict__ aw,
    const mgletreal* __restrict__ ae,
    const mgletreal* __restrict__ as_,
    const mgletreal* __restrict__ an,
    const mgletreal* __restrict__ ab,
    const mgletreal* __restrict__ at,
    const mgletreal* __restrict__ ap,
    const mgletreal* __restrict__ bp,
    const mgletint* __restrict__ mygrids,
    mgletint nmygrids,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d,
    const mgletint* __restrict__ ip1dx,
    const mgletint* __restrict__ ip1dy,
    const mgletint* __restrict__ ip1dz)
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
    const auto ip3 = ip3d[igrid] - 1;
    const auto ipx = ip1dx[igrid] - 1;
    const auto ipy = ip1dy[igrid] - 1;
    const auto ipz = ip1dz[igrid] - 1;

    for (int i = 2 + blockIdx.y; i <= ii - 3; i += gridDim.y)
    {
        for (int j = 2 + threadIdx.y; j <= jj - 3; j += blockDim.y)
        {
            for (int k = 2 + threadIdx.x; k <= kk - 3; k += blockDim.x)
            {
                const auto c = ip3 + k + j * kk + i * kk * jj;
                const auto cw = c - kk * jj;
                const auto ce = c + kk * jj;
                const auto cs = c - kk;
                const auto cn = c + kk;
                const auto cb = c - 1;
                const auto ct = c + 1;

                const auto bp_c = bp[c];

                res[c] = -aw[ipx + i] * phi[cw] * bp[cw] * bp_c - ae[ipx + i] * phi[ce] * bp_c * bp[ce] -
                         as_[ipy + j] * phi[cs] * bp[cs] * bp_c - an[ipy + j] * phi[cn] * bp_c * bp[cn] -
                         ab[ipz + k] * phi[cb] * bp[cb] * bp_c - at[ipz + k] * phi[ct] * bp_c * bp[ct] - ap[c] * phi[c];
            }
        }
    }
}

} // namespace

void laplacephi(
    FArrView<mgletreal> res,
    FArrView<mgletreal> phi,
    FArrView<mgletreal> gsaw,
    FArrView<mgletreal> gsae,
    FArrView<mgletreal> gsas,
    FArrView<mgletreal> gsan,
    FArrView<mgletreal> gsab,
    FArrView<mgletreal> gsat,
    FArrView<mgletreal> gsap,
    FArrView<mgletreal> bp,
    FArrView<mgletint> mygrids,
    FArrView<mgletint> kkk,
    FArrView<mgletint> jjj,
    FArrView<mgletint> iii,
    FArrView<mgletint> ip3d,
    FArrView<mgletint> ip1dx,
    FArrView<mgletint> ip1dy,
    FArrView<mgletint> ip1dz)
{
    const auto nmygrids = mygrids.flat_size();

    const auto threads = ::dim3{32, 8};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygrids)};

    laplacephi_kernel<<<blocks, threads>>>(
        res.device_ptr(),
        phi.device_ptr(),
        gsaw.device_ptr(),
        gsae.device_ptr(),
        gsas.device_ptr(),
        gsan.device_ptr(),
        gsab.device_ptr(),
        gsat.device_ptr(),
        gsap.device_ptr(),
        bp.device_ptr(),
        mygrids.device_ptr(),
        nmygrids,
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr(),
        ip1dx.device_ptr(),
        ip1dy.device_ptr(),
        ip1dz.device_ptr());

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());
}

} // namespace mglet::backend
