include_directories(${CMAKE_SOURCE_DIR}/.subtrees/Catch2/single_include)

add_library(catch_main STATIC "catch_main.cpp")

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

test_gen(apt vec)
test_gen(apt virtual_vec)
test_gen(apt numeric)

# test_gen(particle state)
# test_gen(particle particle)
# test_gen(particle virtual_particle)
# test_gen(particle c_particle)
# test_gen(particle array)
