add_library(backends
    fieldhelper.cpp
    parent2.cpp
    pressuresolver.cpp
)

if(MGLET_OFFLOAD_BACKEND STREQUAL "cuda")
    target_sources(backends 
        PRIVATE cuda/fieldhelper_cuda.cu
        PRIVATE cuda/parent2_cuda.cu
        PRIVATE cuda/pressuresolver_cuda.cu
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
