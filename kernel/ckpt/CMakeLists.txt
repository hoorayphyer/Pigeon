add_library(ckpt SHARED "checkpoint.cpp")
target_include_directories(ckpt PRIVATE ${KERNEL_DIR} )
target_link_libraries( ckpt PRIVATE silopp dye picgen filesys )

# test_gen(io io mpipp filesys silopp)
# target_link_libraries( test_io_io PRIVATE stdc++fs) #TODO this is only for gnu. In GCC 9.1 this is not necessary
