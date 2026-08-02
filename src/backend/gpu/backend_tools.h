#pragma once

#include <iostream>

#if defined(_MGLET_CUDA_)
#include <cuda_runtime.h>
using gpuError_t = cudaError_t;
#define gpuSuccess cudaSuccess
#define gpuGetErrorString cudaGetErrorString
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuGetLastError cudaGetLastError
#elif defined(_MGLET_HIP_)
#include <hip/hip_runtime.h>
using gpuError_t = hipError_t;
#define gpuSuccess hipSuccess
#define gpuGetErrorString hipGetErrorString
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuGetLastError hipGetLastError
#endif

#include "errr.h"

#define GPU_CHECK(func)                                                                                                \
    do                                                                                                                 \
    {                                                                                                                  \
        gpuError_t rt = (func);                                                                                        \
        if (rt != gpuSuccess)                                                                                          \
        {                                                                                                              \
            std::cout << "GPU API call failure: " << rt << " at " << #func << std::endl;                               \
            std::cout << gpuGetErrorString(rt) << std::endl;                                                           \
            MGLET_ERRR();                                                                                              \
        }                                                                                                              \
    } while (0)
