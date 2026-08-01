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

__global__ void accumulate_pcorr_kernel(mgletreal* __restrict__ dp_arr, mgletreal* __restrict__ hilf_arr, std::size_t n)
{
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
    {
        return;
    }

    dp_arr[i] = dp_arr[i] + hilf_arr[i];
}

} // namespace

void accumulate_pcorr_backend(FArrView<mgletreal> dp_arr, const FArrView<mgletreal> hilf_arr)
{
    const auto n_dp = dp_arr.flat_size();
    const auto n_hilf = hilf_arr.flat_size();

    if (n_dp != n_hilf)
    {
        MGLET_ERRR();
    }

    if (n_dp == 0)
    {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n_dp + block_size - 1) / block_size;

    accumulate_pcorr_kernel<<<grid_size, block_size>>>(dp_arr.device_ptr(), hilf_arr.device_ptr(), n_dp);

    CUDA_CHECK(cudaGetLastError());
}

} // namespace mglet::backend
