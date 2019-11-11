add_library(timer SHARED "timer.cpp")

if (${ENABLE_TESTS})
  test_gen(timer timer timer)
endif()
