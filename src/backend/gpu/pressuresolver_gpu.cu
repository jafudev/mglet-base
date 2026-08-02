#include <cstdint>

#if defined(_MGLET_CUDA_)
#include <cuda_runtime.h>
#elif defined(_MGLET_HIP_)
#include <hip/hip_runtime.h>
#endif

#include "errr.h"
#include "mapped_arr_view.h"
#include "gpu_check.h"
#include "arr_tools_interface.h"

namespace mglet::gpu
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
    MappedArrView<mgletreal> maxabsgrid,
    MappedArrView<const mgletreal> phi,
    MappedArrView<const mgletint> mygrids,
    MappedArrView<const mgletint> kkk,
    MappedArrView<const mgletint> jjj,
    MappedArrView<const mgletint> iii,
    MappedArrView<const mgletint> ip3d)
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

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void rescal_backend(
    MappedArrView<mgletreal> rhs,
    MappedArrView<const mgletreal> res,
    MappedArrView<const mgletint> mygrids,
    MappedArrView<const mgletint> kkk,
    MappedArrView<const mgletint> jjj,
    MappedArrView<const mgletint> iii,
    MappedArrView<const mgletint> ip3d)
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

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void sipiter1_hyperplane_level_backend(
    MappedArrView<mgletreal>(res),
    MappedArrView<const mgletreal> rhs,
    MappedArrView<const mgletreal> siplw,
    MappedArrView<const mgletreal> sipls,
    MappedArrView<const mgletreal> siplb,
    MappedArrView<const mgletreal> siplpr,
    MappedArrView<const mgletifk> miphp,
    MappedArrView<const mgletifk> idxhp,
    MappedArrView<const mgletint> mygridsonlvl,
    MappedArrView<const mgletint> kkk,
    MappedArrView<const mgletint> jjj,
    MappedArrView<const mgletint> iii,
    MappedArrView<const mgletint> ip3d)
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

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void sipiter2_hyperplane_level_backend(
    MappedArrView<mgletreal> dp,
    MappedArrView<mgletreal> res,
    MappedArrView<const mgletreal> sipue,
    MappedArrView<const mgletreal> sipun,
    MappedArrView<const mgletreal> siput,
    MappedArrView<const mgletifk> miphp,
    MappedArrView<const mgletifk> idxhp,
    MappedArrView<const mgletint> mygridsonlvl,
    MappedArrView<const mgletint> kkk,
    MappedArrView<const mgletint> jjj,
    MappedArrView<const mgletint> iii,
    MappedArrView<const mgletint> ip3d)
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

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}


} // namespace mglet::gpu

#ifdef _MGLET_USE_BACKEND_

extern "C"
{

void maxabscal_c(
    CFI_cdesc_t* maxabsgrid,
    CFI_cdesc_t* phi,
    CFI_cdesc_t* mygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* idim3d)
{
    mglet::gpu::maxabscal_backend(
        mglet::gpu::MappedArrView<mgletreal>(maxabsgrid),
        mglet::gpu::MappedArrView<const mgletreal>(phi),
        mglet::gpu::MappedArrView<const mgletint>(mygrids),
        mglet::gpu::MappedArrView<const mgletint>(kkk),
        mglet::gpu::MappedArrView<const mgletint>(jjj),
        mglet::gpu::MappedArrView<const mgletint>(iii),
        mglet::gpu::MappedArrView<const mgletint>(idim3d));
}

void accumulate_pcorr_c(CFI_cdesc_t* dp, CFI_cdesc_t* hilf)
{
    mglet::gpu::add_farr_realk(
        mglet::gpu::MappedArrView<mgletreal>(dp), mglet::gpu::MappedArrView<const mgletreal>(hilf));
}

void rescal_c(
    CFI_cdesc_t* rhs,
    CFI_cdesc_t* res,
    CFI_cdesc_t* mygrids,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* idim3d)
{
    mglet::gpu::rescal_backend(
        mglet::gpu::MappedArrView<mgletreal>(rhs),
        mglet::gpu::MappedArrView<const mgletreal>(res),
        mglet::gpu::MappedArrView<const mgletint>(mygrids),
        mglet::gpu::MappedArrView<const mgletint>(kkk),
        mglet::gpu::MappedArrView<const mgletint>(jjj),
        mglet::gpu::MappedArrView<const mgletint>(iii),
        mglet::gpu::MappedArrView<const mgletint>(idim3d));
}

void sipiter1_hyperplane_level_c(
    CFI_cdesc_t* res,
    CFI_cdesc_t* rhs,
    CFI_cdesc_t* siplw,
    CFI_cdesc_t* sipls,
    CFI_cdesc_t* siplb,
    CFI_cdesc_t* siplpr,
    CFI_cdesc_t* miphp,
    CFI_cdesc_t* idxhp,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::gpu::sipiter1_hyperplane_level_backend(
        mglet::gpu::MappedArrView<mgletreal>(res),
        mglet::gpu::MappedArrView<const mgletreal>(rhs),
        mglet::gpu::MappedArrView<const mgletreal>(siplw),
        mglet::gpu::MappedArrView<const mgletreal>(sipls),
        mglet::gpu::MappedArrView<const mgletreal>(siplb),
        mglet::gpu::MappedArrView<const mgletreal>(siplpr),
        mglet::gpu::MappedArrView<const mgletifk>(miphp),
        mglet::gpu::MappedArrView<const mgletifk>(idxhp),
        mglet::gpu::MappedArrView<const mgletint>(mygridsonlvl),
        mglet::gpu::MappedArrView<const mgletint>(kkk),
        mglet::gpu::MappedArrView<const mgletint>(jjj),
        mglet::gpu::MappedArrView<const mgletint>(iii),
        mglet::gpu::MappedArrView<const mgletint>(ip3d));
}

void sipiter2_hyperplane_level_c(
    CFI_cdesc_t* dp,
    CFI_cdesc_t* res,
    CFI_cdesc_t* sipue,
    CFI_cdesc_t* sipun,
    CFI_cdesc_t* siput,
    CFI_cdesc_t* miphp,
    CFI_cdesc_t* idxhp,
    CFI_cdesc_t* mygridsonlvl,
    CFI_cdesc_t* kkk,
    CFI_cdesc_t* jjj,
    CFI_cdesc_t* iii,
    CFI_cdesc_t* ip3d)
{
    mglet::gpu::sipiter2_hyperplane_level_backend(
        mglet::gpu::MappedArrView<mgletreal>(dp),
        mglet::gpu::MappedArrView<mgletreal>(res),
        mglet::gpu::MappedArrView<const mgletreal>(sipue),
        mglet::gpu::MappedArrView<const mgletreal>(sipun),
        mglet::gpu::MappedArrView<const mgletreal>(siput),
        mglet::gpu::MappedArrView<const mgletifk>(miphp),
        mglet::gpu::MappedArrView<const mgletifk>(idxhp),
        mglet::gpu::MappedArrView<const mgletint>(mygridsonlvl),
        mglet::gpu::MappedArrView<const mgletint>(kkk),
        mglet::gpu::MappedArrView<const mgletint>(jjj),
        mglet::gpu::MappedArrView<const mgletint>(iii),
        mglet::gpu::MappedArrView<const mgletint>(ip3d));
}
}

#endif // _MGLET_USE_BACKEND_
