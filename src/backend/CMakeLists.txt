add_library(backend STATIC)

set(GPU_SOURCES
    gpu/bound_pressure_gpu.cu
    gpu/ctof2_gpu.cu
    gpu/conn2_gpu.cu
    gpu/arr_tools.cu
    gpu/laplacephi_gpu.cu
    gpu/parent2_gpu.cu
    gpu/pressuresolver_gpu.cu
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
else()
    # No GPU backend selected: gpu/*.cu above is never compiled, so give the
    # STATIC target the one source CMake needs to determine a link language.
    target_sources(backend
        PRIVATE empty.cpp
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
