add_library( particle_scattering SHARED scattering.cpp )
target_include_directories(particle_scattering PUBLIC ${KERNEL_DIR})
target_link_libraries(particle_scattering PRIVATE picgen)

add_library( particle_forces SHARED forces.cpp )
target_include_directories(particle_forces PUBLIC ${KERNEL_DIR})
target_link_libraries(particle_forces PRIVATE picgen)

add_library( particle_array SHARED array.cpp )
target_include_directories(particle_array PUBLIC ${KERNEL_DIR})
target_link_libraries(particle_array PRIVATE picgen)

add_library(particle_updater SHARED "updater.cpp")
target_link_libraries(particle_updater PRIVATE particle_forces particle_scattering particle_array msh_current_deposition msh_mesh_shape_interplay picgen)

add_library( particle_migration SHARED migration.cpp )
target_include_directories(particle_migration PUBLIC ${KERNEL_DIR} PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(particle_migration PRIVATE mpipp picgen)

if (${ENABLE_TESTS})
  test_gen(particle state) # DONE
  # test_gen(particle particle) # TODO some tests on this are in particle_array
  # test_gen(particle virtual_particle) # TODO
  # test_gen(particle c_particle) # TODO
  test_gen(particle array) # DONE
  test_gen(particle forces) # DONE
  test_gen(particle migration mpipp) # TODO add test on buffer has size 0 and nonzero capacity
  target_include_directories( test_particle_migration PRIVATE ${CMAKE_SOURCE_DIR}) # TODOL this is to include pic.hpp for Specs

  msperf_gen( particle updater )
endif()

