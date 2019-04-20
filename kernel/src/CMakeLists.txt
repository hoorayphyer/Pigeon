find_package(MPI)
add_library(mpipp SHARED "parallel/mpi++.cpp" "parallel/mpi_request.cpp" "parallel/mpi_datatype.cpp" "parallel/mpi_p2p.cpp" "parallel/mpi_collective.cpp")
target_include_directories(mpipp PRIVATE ${MPI_CXX_INCLUDE_DIRS} )
target_include_directories(mpipp PUBLIC ${KERNEL_DIR}/include )
target_compile_options(mpipp PUBLIC "${MPI_CXX_COMPILE_OPTIONS}")
target_link_libraries(mpipp PUBLIC ${MPI_CXX_LIBRARIES})
target_link_libraries(mpipp PRIVATE picgen) # TODO
separate_arguments(MPI_CXX_LINK_OPTIONS NATIVE_COMMAND ${MPI_CXX_LINK_FLAGS}) # a workaround since findMPI doesn't provide this yet
target_link_options(mpipp PUBLIC ${MPI_CXX_LINK_OPTIONS})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/io )
find_package(Silo REQUIRED)
add_library(silo_lib SHARED "io/silo++.cpp" "io/silo_optlist.cpp" "io/silo_operations.cpp")
target_include_directories(silo_lib PRIVATE ${Silo_INCLUDE_DIR} ${KERNEL_DIR}/include )
target_link_libraries( silo_lib PRIVATE ${Silo_LIBRARIES} mpipp)
