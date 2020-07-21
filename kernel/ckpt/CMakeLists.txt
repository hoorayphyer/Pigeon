add_library(ckpt SHARED "checkpoint.cpp")
target_include_directories(ckpt PRIVATE ${KERNEL_DIR} )
target_link_libraries( ckpt PRIVATE silopp dye pic_rc filesys mpi particle_array)

if (${ENABLE_TESTS})
  test_gen(ckpt checkpoint mpi filesys silopp dye stdc++fs ) #TODO this is only for gnu. In GCC 9.1 this is not necessary
endif()

add_executable(ckpt_reader "ckpt_reader.cpp")
target_link_directories(ckpt_reader PRIVATE ${KERNEL_DIR}/../infrastructure/lib)
target_link_libraries(ckpt_reader PRIVATE silopp filesys)
