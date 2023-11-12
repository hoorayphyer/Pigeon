add_library(filesys SHARED "filesys.cpp")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(filesys PRIVATE stdc++fs)
endif()

if (${ENABLE_TESTS})
  test_gen(filesys filesys filesys)
endif()
