add_executable(pic "pic.cpp")
target_include_directories(pic PRIVATE ${KERNEL_DIR} )
target_link_libraries(pic PRIVATE haugbolle_solver old_field_solver field_sync particle_updater particle_migration particle_forces particle_scattering dye mpipp mpi_particle picgen io filesys logger ckpt particle_array timer msh_mesh_shape_interplay)
