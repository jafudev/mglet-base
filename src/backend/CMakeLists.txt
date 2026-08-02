add_library(backend STATIC
    bound_pressure.cpp
    ctof2.cpp
    conn2.cpp
    fieldhelper.cpp
    laplacephi.cpp
    parent2.cpp
    pressuresolver.cpp
)

set(GPU_SOURCES
    gpu/bound_pressure_cuda.cu
    gpu/ctof2_cuda.cu
    gpu/conn2_cuda.cu
    gpu/tooling_cuda.cu
    gpu/laplacephi_cuda.cu
    gpu/parent2_cuda.cu
    gpu/pressuresolver_cuda.cu
)

if(MGLET_OFFLOAD_BACKEND STREQUAL "cuda")
    set_source_files_properties(${GPU_SOURCES} PROPERTIES LANGUAGE CUDA)
    target_sources(backend
        PRIVATE ${GPU_SOURCES}
    )
    target_compile_options(backend
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>
    )
elseif(MGLET_OFFLOAD_BACKEND STREQUAL "hip")
    set_source_files_properties(${GPU_SOURCES} PROPERTIES LANGUAGE HIP)
    target_sources(backend
        PRIVATE ${GPU_SOURCES}
    )
endif()

target_include_directories(backend
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(backend
    PRIVATE corecore
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_Fortran>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},roctx>:rocprofiler-sdk-roctx::rocprofiler-sdk-roctx>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},nvtx>:CUDA::nvtx3>
    PRIVATE $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},cuda>:CUDA::cudart>
    PRIVATE $<$<STREQUAL:${MGLET_OFFLOAD_BACKEND},hip>:hip::device>
)

set_property(TARGET backend PROPERTY CXX_STANDARD 17)
