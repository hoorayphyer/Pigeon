set(CMAKE_CXX_FLAGS "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC") # -Wmissing-braces
add_library(catch_main SHARED "catch_main.cpp")
target_include_directories(catch_main PUBLIC ${CMAKE_SOURCE_DIR}/.subtrees/Catch2/single_include)

add_library(catch_main_mpipp SHARED "catch_main_mpi++.cpp")
target_include_directories(catch_main_mpipp PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_include_directories(catch_main_mpipp PUBLIC ${CMAKE_SOURCE_DIR}/.subtrees/Catch2/single_include)
target_link_libraries(catch_main_mpipp PRIVATE mpipp)

function(test_gen category name)
  set(test_target "test_${category}_${name}")

  add_executable( ${test_target} "${category}/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/include")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/src")
  target_include_directories(${test_target} PRIVATE "${CMAKE_SOURCE_DIR}/tests")

  set(link_main catch_main)
  foreach( lib IN LISTS ARGN)
    target_link_libraries(${test_target} PRIVATE ${lib})
    if( ${lib} STREQUAL "mpipp" )
      set(link_main catch_main_mpipp)
    endif()
  endforeach()
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
# TODO apt particle has newer modules

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

# test_gen(parallel mpi++ mpipp) # DONE

# test_gen(field field) # DONE
# test_gen(field mesh) # DONE
# test_gen(field mesh_shape_interplay) # TODO
# test_gen(field communication mpipp) # DONE
# test_gen(field current_deposition mpipp) # DONE # TODOL test reduce

# test_gen(io silo++ mpipp silo_lib) # TODO

# test_gen(particle state) # DONE
# test_gen(particle particle) # TODO some tests on this are in particle_array
# test_gen(particle virtual_particle) # TODO
# test_gen(particle c_particle) # TODO
# test_gen(particle array) # DONE

# test_gen(particle pusher) # TODO
# test_gen(particle pair_producer)# TODO
test_gen(particle migration mpipp dye) # TODO

test_gen(dye ensemble catch_main_mpipp dye) # TODO
# test_gen(aperture dynamic_balance mpipp)
