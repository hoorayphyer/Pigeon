cmake_minimum_required(VERSION 3.16)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project("UNDETERMINED")

option(ENABLE_TESTS "Enable tests" OFF)
option(ENABLE_DEBUG "Enable debug" OFF)

include(../examples/CMakePreamble.txt)

add_library(pic_rc INTERFACE)
target_include_directories(pic_rc INTERFACE ${CMAKE_CURRENT_LIST_DIR})

# set(SHADOW_DIR ${CMAKE_CURRENT_LIST_DIR}/shadows) # uncomment to enable shadowing
add_subdirectory(../kernel ${CMAKE_BINARY_DIR}/kernel)

include(../examples/CMakePostamble.txt)
