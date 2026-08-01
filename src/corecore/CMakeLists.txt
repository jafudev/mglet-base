add_library(corecore STATIC
    comms_mod.F90
    err_mod.F90
    precision_mod.F90
)

target_include_directories(corecore
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${MPI_Fortran_INCLUDE_PATH}
    PRIVATE ${HDF5_Fortran_INCLUDE_DIRS}
)

target_link_libraries(corecore
    PRIVATE ${MPI_Fortran_LIBRARIES}
    PRIVATE ${HDF5_Fortran_LIBRARIES}
    PRIVATE ZLIB::ZLIB
    PRIVATE $<$<BOOL:${MGLET_OFFLOAD}>:OpenMP::OpenMP_Fortran>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},roctx>:rocprofiler-sdk-roctx::rocprofiler-sdk-roctx>
    PRIVATE $<$<STREQUAL:${MGLET_PROFILE_ANNOTATIONS},nvtx>:CUDA::nvtx3>
)

set_property(TARGET corecore PROPERTY CXX_STANDARD 17)

set_source_files_properties(err_mod.F90 PROPERTIES COMPILE_FLAGS
    $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fall-intrinsics>
)
