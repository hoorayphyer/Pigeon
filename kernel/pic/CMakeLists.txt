add_executable(pic "pic.cpp")
target_precompile_headers(pic PRIVATE "${CMAKE_CURRENT_LIST_DIR}/simulator.hpp" "${CMAKE_CURRENT_LIST_DIR}/vitals.hpp" "${CMAKE_CURRENT_LIST_DIR}/modules.hpp" "${CMAKE_CURRENT_LIST_DIR}/argparser.hpp")
target_include_directories(pic PRIVATE ${KERNEL_DIR} ${CMAKE_CURRENT_LIST_DIR}/../../.subtrees/toml/include)
target_link_libraries(pic PRIVATE old_field_solver field_sync particle_updater particle_migration particle_forces particle_scattering particle_annihilation dye mpi mpi_particle pic_rc io filesys logger ckpt particle_array timer msh_mesh_shape_interplay)
