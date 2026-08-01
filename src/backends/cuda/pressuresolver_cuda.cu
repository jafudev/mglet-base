#include "fieldhelper_backend.h"

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

__global__ void accumulate_pcorr_kernel(mgletreal* __restrict__ dp_arr, mgletreal* __restrict__ hilf_arr, std::size_t n)
{
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
    {
        return;
    }

    dp_arr[i] = dp_arr[i] + hilf_arr[i];
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

} // namespace

void maxabscal_backend(
    FArrView<mgletreal> maxabsgrid,
    const FArrView<mgletreal> phi,
    const FArrView<mgletint> mygrids,
    mgletint nmygrids,
    const FArrView<mgletint> kkk,
    const FArrView<mgletint> jjj,
    const FArrView<mgletint> iii,
    const FArrView<mgletint> ip3d)
{
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

void accumulate_pcorr_backend(FArrView<mgletreal> dp_view, const FArrView<mgletreal> hilf_view)
{
    const auto n_dp = dp_view.flat_size();
    const auto n_hilf = hilf_view.flat_size();

    if (n_dp != n_hilf)
    {
        MGLET_ERRR();
    }

    if (n_dp == 0)
    {
        return;
    }

    const unsigned threads = 256;
    const unsigned blocks = (n_dp + threads - 1) / threads;

    accumulate_pcorr_kernel<<<blocks, threads>>>(dp_view.device_ptr(), hilf_view.device_ptr(), n_dp);

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void rescal_backend(
    FArrView<mgletreal> rhs_view,
    const FArrView<mgletreal> res_view,
    mgletint nmygrids,
    const FArrView<mgletint> mygrids_view,
    const FArrView<mgletint> kkk_view,
    const FArrView<mgletint> jjj_view,
    const FArrView<mgletint> iii_view,
    const FArrView<mgletint> ip3d_view)
{
    const auto n_rhs = rhs_view.flat_size();
    const auto n_res = res_view.flat_size();

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
        rhs_view.device_ptr(),
        res_view.device_ptr(),
        nmygrids,
        mygrids_view.device_ptr(),
        kkk_view.device_ptr(),
        jjj_view.device_ptr(),
        iii_view.device_ptr(),
        ip3d_view.device_ptr());

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

} // namespace mglet::backend
