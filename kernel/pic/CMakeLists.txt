add_executable(pic "pic.cpp")
target_include_directories(pic PRIVATE ${KERNEL_DIR} )
target_link_libraries(pic PRIVATE haugbolle_solver old_field_solver field_sync particle_updater particle_migration particle_forces particle_scattering dye mpipp mpi_particle picgen io filesys logger ckpt particle_array timer msh_mesh_shape_interplay)

# FIXME not sure why the RPATH info in mpipp is not transitive here. Build? Install?
find_package(MPI)

target_compile_options(pic PUBLIC "${MPI_CXX_COMPILE_OPTIONS}")
separate_arguments(MPI_CXX_LINK_OPTIONS NATIVE_COMMAND ${MPI_CXX_LINK_FLAGS}) # a workaround since findMPI doesn't provide this yet
target_link_options(pic PUBLIC ${MPI_CXX_LINK_OPTIONS})
target_link_libraries(pic PRIVATE ${MPI_CXX_LIBRARIES})
