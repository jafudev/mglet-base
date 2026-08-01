#pragma once
#include <iostream>

#include <cuda_runtime.h>

#include "errr.h"

// https://github.com/NVIDIA/CUDALibrarySample
#ifndef CUDA_CHECK
#define CUDA_CHECK(func)                                                       \
do {                                                                           \
    cudaError_t rt = (func);                                                   \
    if (rt != cudaSuccess) {                                                   \
        std::cout << "CUDA API call failure: " << rt << " at " << #func        \
            << std::endl;                                                      \
        MGLET_ERRR();                                                          \
    }                                                                          \
} while (0)
#endif
