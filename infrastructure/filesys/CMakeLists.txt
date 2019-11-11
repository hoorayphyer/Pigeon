add_library(filesys SHARED "filesys.cpp")
target_link_libraries(filesys PRIVATE stdc++fs) #FIXME this is only for gnu. Not needed in GCC 9

if (${ENABLE_TESTS})
  test_gen(filesys filesys filesys)
endif()
