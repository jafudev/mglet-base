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

__global__ void set_field_arr_realk_kernel(mgletreal* __restrict__ arr, std::size_t n, mgletreal val)
{
    const auto i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n)
    {
        return;
    }

    arr[i] = val;
}

} // namespace

void set_field_arr_realk_backend(FArrView<mgletreal> arr, mgletreal val)
{
    const auto n = arr.flat_size();

    if (n == 0)
    {
        return;
    }

    const unsigned block_size = 256;
    const unsigned grid_size = (n + block_size - 1) / block_size;

    set_field_arr_realk_kernel<<<grid_size, block_size>>>(arr.device_ptr(), n, val);

    CUDA_CHECK(cudaGetLastError());
}

} // namespace mglet::backend
