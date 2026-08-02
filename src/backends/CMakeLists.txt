add_library(backends
    bound_pressure.cpp
    ctof2.cpp
    conn2.cpp
    fieldhelper.cpp
    laplacephi.cpp
    parent2.cpp
    pressuresolver.cpp
)

if(MGLET_OFFLOAD_BACKEND STREQUAL "cuda")
    target_sources(backends
        PRIVATE cuda/bound_pressure_cuda.cu
        PRIVATE cuda/ctof2_cuda.cu
        PRIVATE cuda/conn2_cuda.cu
        PRIVATE cuda/tooling_cuda.cu
        PRIVATE cuda/laplacephi_cuda.cu
        PRIVATE cuda/parent2_cuda.cu
        PRIVATE cuda/pressuresolver_cuda.cu
    )
    target_compile_options(backends
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>
    )
elseif(MGLET_OFFLOAD_BACKEND STREQUAL "hip")
    target_sources(backends PRIVATE hip/fieldhelper_hip.hip)
endif()

target_include_directories(backends
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(backends
    PRIVATE corecore
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_Fortran>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},roctx>:rocprofiler-sdk-roctx::rocprofiler-sdk-roctx>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},nvtx>:CUDA::nvtx3>
    PRIVATE $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},cuda>:CUDA::cudart>
    PRIVATE $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},hip>:hip::device>
)

set_property(TARGET backends PROPERTY CXX_STANDARD 17)
