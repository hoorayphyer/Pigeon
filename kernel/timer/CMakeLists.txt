add_library(timer SHARED "timer.cpp")
target_include_directories(timer PRIVATE ${KERNEL_DIR})

if (${ENABLE_TESTS})
  test_gen(timer timer timer)
endif()
