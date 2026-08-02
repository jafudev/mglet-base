#include <cstdint>

#if defined(_MGLET_CUDA_)
#include <cuda_runtime.h>
#elif defined(_MGLET_HIP_)
#include <hip/hip_runtime.h>
#endif

#include "errr.h"
#include "gpu_check.h"
#include "gpu_tools_interface.h"
#include "mapped_arr_view.h"

namespace mglet::gpu
{

namespace
{

__global__ void set_farr_realk_kernel(mgletreal* __restrict__ arr, mgletint n, mgletreal val)
{
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n)
    {
        return;
    }

    arr[i] = val;
}

__global__ void set_farr_ifk_kernel(mgletifk* __restrict__ arr, mgletint n, mgletifk val)
{
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n)
    {
        return;
    }

    arr[i] = val;
}

__global__ void add_farr_realk_kernel(
    mgletreal* __restrict__ lhs,
    const mgletreal* __restrict__ rhs,
    mgletint n)
{
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n)
    {
        return;
    }

    lhs[i] = lhs[i] + rhs[i];
}

} // namespace

void set_farr_realk(MappedArrView<mgletreal> farr, mgletreal val)
{
    const auto n = farr.flat_size();

    if (n == 0)
    {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n + block_size - 1) / block_size;

    set_farr_realk_kernel<<<grid_size, block_size>>>(farr.device_ptr(), n, val);

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void set_farr_ifk(MappedArrView<mgletifk> farr, mgletifk val)
{
    const auto n = farr.flat_size();

    if (n == 0)
    {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n + block_size - 1) / block_size;

    set_farr_ifk_kernel<<<grid_size, block_size>>>(farr.device_ptr(), n, val);

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

void add_farr_realk(MappedArrView<mgletreal> lhs,  MappedArrView<const mgletreal> rhs)
{
    const auto n_lhs = lhs.flat_size();
    const auto n_rhs = rhs.flat_size();

    if (n_lhs != n_rhs)
    {
        MGLET_ERRR();
    }

    if (n_lhs == 0)
    {
        return;
    }

    const unsigned threads = 256;
    const unsigned blocks = (n_lhs + threads - 1) / threads;

    add_farr_realk_kernel<<<blocks, threads>>>(lhs.device_ptr(), rhs.device_ptr(), n_lhs);

    GPU_CHECK(gpuGetLastError());
    GPU_CHECK(gpuDeviceSynchronize());
}

} // namespace mglet::gpu

#ifdef _MGLET_USE_BACKEND_

extern "C" void set_field_arr_realk_c(CFI_cdesc_t* farr, mgletreal val)
{
    mglet::gpu::set_farr_realk(mglet::gpu::MappedArrView<mgletreal>(farr), val);
}

extern "C" void set_field_arr_ifk_c(CFI_cdesc_t* farr, mgletifk val)
{
    mglet::gpu::set_farr_ifk(mglet::gpu::MappedArrView<mgletifk>(farr), val);
}

extern "C" void accumulate_pcorr_c(CFI_cdesc_t* dp, CFI_cdesc_t* hilf)
{
    mglet::gpu::add_farr_realk(
        mglet::gpu::MappedArrView<mgletreal>(dp), mglet::gpu::MappedArrView<const mgletreal>(hilf));
}

#endif // _MGLET_USE_BACKEND_
