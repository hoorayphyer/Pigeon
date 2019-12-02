cmake_minimum_required(VERSION 3.14)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project (InfraStructure)
option(ENABLE_TESTS "Enable tests" OFF) #FIXME

# message ("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
if (${ENABLE_TESTS})
  enable_testing()
endif()

# Set CXX flags
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-std=c++17 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC ${CXX_EXTRA_FLAGS}")
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE_INIT "-xHost -no-prec-div")
else()
  string(APPEND CMAKE_CXX_FLAGS_RELEASE_INIT "-march=native")
endif()

set(CMAKE_CXX_FLAGS_DEBUG_INIT "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE_INIT})
# set(CMAKE_CXX_COMPILER "g++")

# enable colored cmake message.
# Copied from https://stackoverflow.com/questions/18968979/how-to-get-colorized-output-with-cmake
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(Yellow      "${Esc}[33m")
endif()

set(KERNEL_DIR ${CMAKE_CURRENT_LIST_DIR})

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
endif()

function(test_gen category name)
  set(test_target "test_${category}_${name}")
  add_executable( ${test_target} "tests/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE ${PROJECT_SOURCE_DIR})
  target_compile_options( ${test_target} PRIVATE -O0 ) # -Wmissing-braces

  set(link_testfw testfw)
  foreach( lib IN LISTS ARGN)
    target_link_libraries(${test_target} PRIVATE ${lib})
    if( ${lib} STREQUAL "mpipp" )
      set(link_testfw testfw_mpi)
    endif()
  endforeach()
  target_link_libraries( ${test_target} PRIVATE ${link_testfw} )
  add_test( NAME TEST_${category}_${name} COMMAND ${test_target} )
endfunction()

add_subdirectory(logger)
add_subdirectory(apt)
add_subdirectory(timer)
add_subdirectory(filesys)
add_subdirectory(mpipp)
add_subdirectory(silopp)
add_subdirectory(random)
add_subdirectory(testfw)
