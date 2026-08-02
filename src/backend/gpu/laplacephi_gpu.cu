#include <cstdint>

#include "errr.h"
#include "gpu_check.h"
#include "gpu_include.h"
#include "mapped_arr_view.h"

namespace mglet::gpu
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

__global__ void laplacephi_level_kernel(
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
    const mgletint* __restrict__ mygridsonlvl,
    mgletint nmygridsonlvl,
    const mgletint* __restrict__ kkk,
    const mgletint* __restrict__ jjj,
    const mgletint* __restrict__ iii,
    const mgletint* __restrict__ ip3d,
    const mgletint* __restrict__ ip1dx,
    const mgletint* __restrict__ ip1dy,
    const mgletint* __restrict__ ip1dz)
{
    const auto block_idx = blockIdx.x;
    if (block_idx >= nmygridsonlvl)
    {
        return;
    }

    // mygridsonlvl is zero-padded past this level's true grid count, see
    // mygridslvl in grids_mod.F90. Skip the padding.
    const auto grid_id = mygridsonlvl[block_idx];
    if (grid_id == 0)
    {
        return;
    }

    const auto igrid = grid_id - 1; // C is 0-based

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
    MappedArrView<mgletreal> res,
    MappedArrView<const mgletreal> phi,
    MappedArrView<const mgletreal> gsaw,
    MappedArrView<const mgletreal> gsae,
    MappedArrView<const mgletreal> gsas,
    MappedArrView<const mgletreal> gsan,
    MappedArrView<const mgletreal> gsab,
    MappedArrView<const mgletreal> gsat,
    MappedArrView<const mgletreal> gsap,
    MappedArrView<const mgletreal> bp,
    MappedArrView<const mgletint> mygrids,
    MappedArrView<const mgletint> kkk,
    MappedArrView<const mgletint> jjj,
    MappedArrView<const mgletint> iii,
    MappedArrView<const mgletint> ip3d,
    MappedArrView<const mgletint> ip1dx,
    MappedArrView<const mgletint> ip1dy,
    MappedArrView<const mgletint> ip1dz)
{
    const auto nmygrids = mygrids.flat_size();

    if (nmygrids == 0)
    {
        return;
    }

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

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void laplacephi_level(
    MappedArrView<mgletreal> res,
    MappedArrView<const mgletreal> phi,
    MappedArrView<const mgletreal> gsaw,
    MappedArrView<const mgletreal> gsae,
    MappedArrView<const mgletreal> gsas,
    MappedArrView<const mgletreal> gsan,
    MappedArrView<const mgletreal> gsab,
    MappedArrView<const mgletreal> gsat,
    MappedArrView<const mgletreal> gsap,
    MappedArrView<const mgletreal> bp,
    MappedArrView<const mgletint> mygridsonlvl,
    MappedArrView<const mgletint> kkk,
    MappedArrView<const mgletint> jjj,
    MappedArrView<const mgletint> iii,
    MappedArrView<const mgletint> ip3d,
    MappedArrView<const mgletint> ip1dx,
    MappedArrView<const mgletint> ip1dy,
    MappedArrView<const mgletint> ip1dz)
{
    const auto nmygridsonlvl = mygridsonlvl.flat_size();

    if (nmygridsonlvl == 0)
    {
        return;
    }

    const auto threads = ::dim3{32, 8};
    const auto blocks = ::dim3{static_cast<unsigned>(nmygridsonlvl)};

    laplacephi_level_kernel<<<blocks, threads>>>(
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
        mygridsonlvl.device_ptr(),
        nmygridsonlvl,
        kkk.device_ptr(),
        jjj.device_ptr(),
        iii.device_ptr(),
        ip3d.device_ptr(),
        ip1dx.device_ptr(),
        ip1dy.device_ptr(),
        ip1dz.device_ptr());

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

} // namespace mglet::gpu

#ifdef _MGLET_USE_BACKEND_

extern "C"
{

void laplacephi_c(
    CFI_cdesc_t* res,
    CFI_cdesc_t* phi,
    CFI_cdesc_t* gsaw,
    CFI_cdesc_t* gsae,
    CFI_cdesc_t* gsas,
    CFI_cdesc_t* gsan,
    CFI_cdesc_t* gsab,
    CFI_cdesc_t* gsat,
    CFI_cdesc_t* gsap,
    CFI_cdesc_t* bp,
    CFI_cdesc_t* mygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d,
    CFI_cdesc_t* ip1dx,
    CFI_cdesc_t* ip1dy,
    CFI_cdesc_t* ip1dz)
{
    mglet::gpu::laplacephi(
        mglet::gpu::MappedArrView<mgletreal>(res),
        mglet::gpu::MappedArrView<const mgletreal>(phi),
        mglet::gpu::MappedArrView<const mgletreal>(gsaw),
        mglet::gpu::MappedArrView<const mgletreal>(gsae),
        mglet::gpu::MappedArrView<const mgletreal>(gsas),
        mglet::gpu::MappedArrView<const mgletreal>(gsan),
        mglet::gpu::MappedArrView<const mgletreal>(gsab),
        mglet::gpu::MappedArrView<const mgletreal>(gsat),
        mglet::gpu::MappedArrView<const mgletreal>(gsap),
        mglet::gpu::MappedArrView<const mgletreal>(bp),
        mglet::gpu::MappedArrView<const mgletint>(mygrids),
        mglet::gpu::MappedArrView<const mgletint>(kkk),
        mglet::gpu::MappedArrView<const mgletint>(jjj),
        mglet::gpu::MappedArrView<const mgletint>(iii),
        mglet::gpu::MappedArrView<const mgletint>(ip3d),
        mglet::gpu::MappedArrView<const mgletint>(ip1dx),
        mglet::gpu::MappedArrView<const mgletint>(ip1dy),
        mglet::gpu::MappedArrView<const mgletint>(ip1dz));
}


void laplacephi_level_c(
    CFI_cdesc_t* res,
    CFI_cdesc_t* phi,
    CFI_cdesc_t* gsaw,
    CFI_cdesc_t* gsae,
    CFI_cdesc_t* gsas,
    CFI_cdesc_t* gsan,
    CFI_cdesc_t* gsab,
    CFI_cdesc_t* gsat,
    CFI_cdesc_t* gsap,
    CFI_cdesc_t* bp,
    mgletint nmygridsonlvl,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d,
    CFI_cdesc_t* ip1dx,
    CFI_cdesc_t* ip1dy,
    CFI_cdesc_t* ip1dz)
{
    mglet::gpu::laplacephi_level(
        mglet::gpu::MappedArrView<mgletreal>(res),
        mglet::gpu::MappedArrView<const mgletreal>(phi),
        mglet::gpu::MappedArrView<const mgletreal>(gsaw),
        mglet::gpu::MappedArrView<const mgletreal>(gsae),
        mglet::gpu::MappedArrView<const mgletreal>(gsas),
        mglet::gpu::MappedArrView<const mgletreal>(gsan),
        mglet::gpu::MappedArrView<const mgletreal>(gsab),
        mglet::gpu::MappedArrView<const mgletreal>(gsat),
        mglet::gpu::MappedArrView<const mgletreal>(gsap),
        mglet::gpu::MappedArrView<const mgletreal>(bp),
        mglet::gpu::MappedArrView<const mgletint>(mygridsonlvl),
        mglet::gpu::MappedArrView<const mgletint>(kkk),
        mglet::gpu::MappedArrView<const mgletint>(jjj),
        mglet::gpu::MappedArrView<const mgletint>(iii),
        mglet::gpu::MappedArrView<const mgletint>(ip3d),
        mglet::gpu::MappedArrView<const mgletint>(ip1dx),
        mglet::gpu::MappedArrView<const mgletint>(ip1dy),
        mglet::gpu::MappedArrView<const mgletint>(ip1dz));
}

}

#endif // _MGLET_USE_BACKEND_
