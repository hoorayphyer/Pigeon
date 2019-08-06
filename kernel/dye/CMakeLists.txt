add_library(dye SHARED "ensemble.cpp" "dynamic_balance.cpp")
target_include_directories(dye PRIVATE ${KERNEL_DIR})
target_link_libraries(dye PRIVATE particle_array mpipp mpi_particle picgen)

if (${ENABLE_TESTS})
  test_gen(dye ensemble mpipp) # DONE
  test_gen(dye dynamic_balance mpipp) # TODO
endif()
