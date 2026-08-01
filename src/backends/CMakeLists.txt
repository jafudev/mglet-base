add_library(backends
    fieldhelper_backend.cpp
)

target_include_directories(backends
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(backends
    PRIVATE corecore
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_Fortran>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},roctx>:rocprofiler-sdk-roctx::rocprofiler-sdk-roctx>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},nvtx>:CUDA::nvtx3>
    PRIVATE $<$<STREQUAL:${MGLET_GPU_BACKEND},CUDA>:CUDA::cudart>
    PRIVATE $<$<STREQUAL:${MGLET_GPU_BACKEND},HIP>:hip::device>
)

set_property(TARGET backends PROPERTY CXX_STANDARD 17)
