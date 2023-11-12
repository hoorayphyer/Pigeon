add_library(io SHARED "io.cpp")
target_include_directories(io PRIVATE ${KERNEL_DIR} )
target_compile_definitions(io PRIVATE APPARENT_DATA_DIR="${CMAKE_SOURCE_DIR}/Data")
target_link_libraries( io PRIVATE silopp dye filesys field_sync mpi )

if (${ENABLE_TESTS})
  test_gen(io io mpi filesys silopp) # DONE
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries( test_io_io PRIVATE stdc++fs)
  endif()
endif()
