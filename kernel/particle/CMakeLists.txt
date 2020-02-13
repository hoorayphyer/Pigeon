add_library( particle_scattering SHARED scattering.cpp )
target_include_directories(particle_scattering PUBLIC ${KERNEL_DIR})
target_link_libraries(particle_scattering PRIVATE pic_rc)

add_library( particle_forces SHARED forces.cpp )
target_include_directories(particle_forces PUBLIC ${KERNEL_DIR})
target_link_libraries(particle_forces PRIVATE pic_rc)

add_library( particle_array SHARED array.cpp )
target_include_directories(particle_array PUBLIC ${KERNEL_DIR})
target_link_libraries(particle_array PRIVATE pic_rc)

add_library(particle_updater SHARED "updater.cpp")
target_link_libraries(particle_updater PRIVATE particle_forces particle_scattering particle_array msh_current_deposition msh_mesh_shape_interplay pic_rc)

add_library(mpi_particle SHARED "mpi_particle.cpp")
target_include_directories(mpi_particle PUBLIC ${KERNEL_DIR})
target_link_libraries(mpi_particle PRIVATE mpi pic_rc)

add_library( particle_migration SHARED migration.cpp )
target_include_directories(particle_migration PUBLIC ${KERNEL_DIR} PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(particle_migration PRIVATE mpi mpi_particle pic_rc)

if (${ENABLE_TESTS})
  test_gen(particle state)
  # test_gen(particle particle) # TODO some tests on this are in particle_array
  # test_gen(particle virtual_particle) # TODO
  test_gen(particle array)
  test_gen(particle map)
  test_gen(particle forces)
  test_gen(particle sorter)
  test_gen(particle migration mpi) # TODO add test on buffer has size 0 and nonzero capacity

  msperf_gen( particle updater )
endif()

