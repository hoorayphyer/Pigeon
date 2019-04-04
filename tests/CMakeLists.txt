set(CMAKE_CXX_FLAGS "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC") # -Wmissing-braces
add_library(catch_main STATIC "catch_main.cpp")
target_include_directories(catch_main PUBLIC ${CMAKE_SOURCE_DIR}/.subtrees/Catch2/single_include)

add_library(catch_main_mpipp STATIC "catch_main_mpi++.cpp")
target_include_directories(catch_main_mpipp PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_include_directories(catch_main_mpipp PUBLIC ${CMAKE_SOURCE_DIR}/.subtrees/Catch2/single_include)
target_link_libraries(catch_main_mpipp PRIVATE mpipp_lib)

function(test_gen category name)
  if ( ${ARGC} GREATER 2 )
    set( link_main ${ARGV2} ) # ARGV0, ARGV1, ARGV2..
  else()
    set( link_main catch_main )
  endif()

  set(test_target "test_${category}_${name}")

  add_executable( ${test_target} "${category}/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/include")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/src")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/tests")

  target_link_libraries(${test_target} PRIVATE ${link_main})
  set_target_properties( ${test_target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../test_bin"
    )

  add_test( NAME TEST_${category}_${name} COMMAND ${test_target} )

  unset(test_target)
endfunction()

# TODO all uses of std::abs in the code, see if anywhere std::labs should be used
# TODO all mpi test modules should have larger tests

# test_gen(apt array) # DONE
# test_gen(apt pair) # DONE
# test_gen(apt "foreach") # DONE
# test_gen(apt block) # DONE
# test_gen(apt vec) # DONE
# test_gen(apt virtual_vec) # TODO
# test_gen(apt numeric) # DONE
# test_gen(apt handle) # DONE


# test_gen(kernel coordinate) # TODO slight mismatched in phi
# test_gen(kernel grid) # TODO

# test_gen(parallel mpi++ catch_main_mpipp) # DONE

# test_gen(field mesh) # DONE
# test_gen(field mesh_shape_interplay) # TODO
# test_gen(field communication catch_main_mpipp) # DONE
# test_gen(field current_deposition catch_main_mpipp) # TODO

# test_gen(io silo++ catch_main_mpipp)
# target_link_libraries( test_io_silo++ PRIVATE silo_lib )

# test_gen(particle state) # DONE
# test_gen(particle particle) # TODO some tests on this are in particle_array
# test_gen(particle virtual_particle) # TODO
# test_gen(particle c_particle) # TODO
# test_gen(particle array) # DONE

# test_gen(particle pusher) # TODO
# test_gen(particle migration catch_main_mpipp) # TODO
# test_gen(particle pair_producer)# TODO

# test_gen(aperture ensemble catch_main_mpipp)
# test_gen(aperture dynamic_balance catch_main_mpipp)
