add_library(timer SHARED "timer.cpp")
target_include_directories(timer PRIVATE ${KERNEL_DIR})

test_gen(timer timer timer)
