add_library(dye SHARED "ensemble.cpp" "dynamic_balance.cpp")
target_include_directories(dye PRIVATE ${KERNEL_DIR})
target_link_libraries(dye PRIVATE particle_array mpi mpi_particle picgen)

if (${ENABLE_TESTS})
  test_gen(dye ensemble mpi) # DONE
  test_gen(dye dynamic_balance mpi) # TODO
endif()
