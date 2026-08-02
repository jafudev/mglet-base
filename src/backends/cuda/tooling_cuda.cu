#include <cstddef>

#include <cuda_runtime.h>

#include "cutools.h"
#include "f_arr_view.h"

namespace mglet::backend
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

void set_farr_realk(FArrView<mgletreal> farr, mgletreal val)
{
    const auto n = farr.flat_size();

    if (n == 0)
    {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n + block_size - 1) / block_size;

    set_farr_realk_kernel<<<grid_size, block_size>>>(farr.device_ptr(), n, val);

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());
}

void set_farr_ifk(FArrView<mgletifk> farr, mgletifk val)
{
    const auto n = farr.flat_size();

    if (n == 0)
    {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n + block_size - 1) / block_size;

    set_farr_ifk_kernel<<<grid_size, block_size>>>(farr.device_ptr(), n, val);

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaGetLastError());
}

void add_farr_realk(FArrView<mgletreal> lhs,  FArrView<const mgletreal> rhs)
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

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

} // namespace mglet::backend
