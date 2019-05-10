add_library(io SHARED "io.cpp")
target_include_directories(io PRIVATE ${KERNEL_DIR} )
target_compile_definitions(io PRIVATE APPARENT_DATA_DIR="${CMAKE_SOURCE_DIR}/Data")
target_link_libraries( io PRIVATE silopp dye picgen filesys field_communication )

test_gen(io io mpipp filesys silopp) # DONE
target_link_libraries( test_io_io PRIVATE stdc++fs) #TODO this is only for gnu. In GCC 9.1 this is not necessary
