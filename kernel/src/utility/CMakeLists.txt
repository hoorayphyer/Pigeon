add_library(filesys SHARED "filesys.cpp")
target_include_directories(filesys PRIVATE ${KERNEL_DIR}/include )
target_link_libraries(filesys PRIVATE stdc++fs) #TODO this is only for gnu
