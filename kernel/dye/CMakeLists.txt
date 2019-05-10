add_library(dye SHARED "ensemble.cpp" "dynamic_balance.cpp")
target_include_directories(dye PRIVATE ${KERNEL_DIR})
target_link_libraries(dye PRIVATE particle_array mpipp picgen)

test_gen(dye ensemble mpipp dye) # DONE
test_gen(dye dynamic_balance mpipp dye) # TODO
target_include_directories( test_dye_dynamic_balance PRIVATE ${CMAKE_SOURCE_DIR})
