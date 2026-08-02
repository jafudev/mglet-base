#pragma once

#include <iostream>

#include "errr.h"
#include "gpu_include.h"

#if defined(_MGLET_CUDA_)
using gpuError_t = cudaError_t;
#define gpuSuccess cudaSuccess
#define gpuGetErrorString cudaGetErrorString
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuGetLastError cudaGetLastError
#elif defined(_MGLET_HIP_)
using gpuError_t = hipError_t;
#define gpuSuccess hipSuccess
#define gpuGetErrorString hipGetErrorString
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuGetLastError hipGetLastError
#endif

#ifdef _MGLET_USE_BACKEND_
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
#endif
