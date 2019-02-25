include_directories(${CMAKE_SOURCE_DIR}/.subtrees/Catch2/single_include)

add_library(catch_main STATIC "catch_main.cpp")

set(CMAKE_CXX_FLAGS "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")

macro(test_gen category name)
  set(test_target "test_${category}_${name}")

  add_executable( ${test_target} "${category}/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/include")
  target_include_directories(${test_target} PRIVATE "../src")
  target_link_libraries(${test_target} PRIVATE catch_main)
  set_target_properties( ${test_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../test_bin"
    )

  add_test( NAME TEST_${category}_${name} COMMAND ${test_target} )

  unset(test_target)
endmacro()

function(link_mpi target)
  target_include_directories(${target} PUBLIC ${MPI_INCLUDE_PATH})
  target_link_libraries(${target} PRIVATE ${MPI_LIBRARIES})
  set_target_properties(${target} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  set_target_properties(${target} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endfunction()
add_library(catch_main_mpi STATIC "catch_main_mpi.cpp")
link_mpi(catch_main_mpi)

macro(mpi_test_gen category name)
  set(test_target "test_${category}_${name}")

  add_executable( ${test_target} "${category}/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/include")
  target_include_directories(${test_target} PRIVATE "../src")
  target_link_libraries(${test_target} PRIVATE catch_main_mpi)
  set_target_properties( ${test_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../test_bin"
    )

  add_test( NAME TEST_${category}_${name} COMMAND ${test_target} )

  unset(test_target)
endmacro()


test_gen(apt vec)
test_gen(apt virtual_vec)
test_gen(apt numeric)
test_gen(apt handle)

test_gen(kernel coordinate) # TODO
# test_gen(kernel grid_shape_interplay)

mpi_test_gen(parallel mpi++)
# link_mpi(test_parallel_mpi++)
mpi_test_gen(parallel locale)
# link_mpi(test_parallel_locale)

test_gen(particle state)
# test_gen(particle particle)
# test_gen(particle virtual_particle)
# test_gen(particle c_particle)
test_gen(particle array)
