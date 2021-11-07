add_library(confFile SHARED "conf_file.cpp")
target_include_directories(confFile BEFORE PRIVATE ${KERNEL_DIR} ${CMAKE_CURRENT_LIST_DIR}/../../.subtrees/toml ${CMAKE_CURRENT_LIST_DIR}/../../.subtrees/fmt/include)

add_library(simulationBuilder SHARED "builder.cpp")
target_include_directories(simulationBuilder BEFORE PRIVATE ${KERNEL_DIR})
# TODO pic_rc for pic.hpp
target_link_libraries(simulationBuilder PRIVATE mpi io logger ckpt pic_rc particle_array)


add_library(simulatior SHARED "simulator.cpp")
target_include_directories(simulatior BEFORE PRIVATE ${KERNEL_DIR})
# TODO pic_rc for pic.hpp
target_link_libraries(simulatior PRIVATE field_sync dye mpi io filesys logger ckpt pic_rc particle_array timer)

# add_library(simulator SHARED "simulator.cpp")
# target_include_directories(simulator BEFORE PRIVATE ${KERNEL_DIR} ${CMAKE_CURRENT_LIST_DIR}/../../.subtrees/toml/include)
# target_link_libraries(simulator PRIVATE log_spherical_solver field_sync particle_updater particle_migration particle_forces particle_scattering particle_annihilation dye mpi mpi_particle pic_rc io filesys logger ckpt particle_array timer msh_mesh_shape_interplay tracing silopp)
