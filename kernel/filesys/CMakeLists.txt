add_library(filesys SHARED "filesys.cpp")
target_link_libraries(filesys PRIVATE stdc++fs) #TODO this is only for gnu

if (${ENABLE_TESTS})
  test_gen(filesys filesys filesys)
endif()
