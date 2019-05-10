set(CMAKE_CXX_FLAGS "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC") # -Wmissing-braces # TODO add this to target public properties
add_library(testfw SHARED "catch_main.cpp")
target_include_directories(testfw PUBLIC ${KERNEL_DIR}/apt ${KERNEL_DIR}/../.subtrees/Catch2/single_include)

add_library(testfw_mpi SHARED "catch_main_mpi++.cpp")
target_include_directories(testfw_mpi PUBLIC ${KERNEL_DIR} ${KERNEL_DIR}/../.subtrees/Catch2/single_include )
target_link_libraries(testfw_mpi PRIVATE mpipp)
