find_package(MPI)
add_library(mpipp SHARED "mpi++.cpp" "mpi_request.cpp" "mpi_datatype.cpp" "mpi_p2p.cpp" "mpi_collective.cpp")
target_include_directories(mpipp PUBLIC ${MPI_CXX_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR} )
target_compile_options(mpipp PUBLIC "${MPI_CXX_COMPILE_OPTIONS}")
target_link_libraries(mpipp PUBLIC ${MPI_CXX_LIBRARIES})
separate_arguments(MPI_CXX_LINK_OPTIONS NATIVE_COMMAND ${MPI_CXX_LINK_FLAGS}) # a workaround since findMPI doesn't provide this yet
target_link_options(mpipp PUBLIC ${MPI_CXX_LINK_OPTIONS})

if (${ENABLE_TESTS})
  test_gen(mpipp mpi++ mpipp) # DONE
endif()
