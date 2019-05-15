cmake_minimum_required(VERSION 3.13)
project () # specify your project name here
enable_testing()
message ("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

# Set CXX flags
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC -xHost -no-prec-div ${CXX_EXTRA_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC -march=native ${CXX_EXTRA_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
# set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_FLAGS_DEBUG "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
endif()

function(add_custom_test name)
  set(test_target "test_custom_${name}")
  add_executable( ${test_target} "./tests/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE "../kernel" ${CMAKE_CURRENT_LIST_DIR} )

  set(link_testfw testfw)
  foreach( lib IN LISTS ARGN)
    target_link_libraries(${test_target} PRIVATE ${lib})
    if( ${lib} STREQUAL "mpipp" )
      set(link_testfw testfw_mpi)
    endif()
  endforeach()
  target_link_libraries( ${test_target} PRIVATE ${link_testfw} )
  add_test( NAME TEST_CUSTOM_${name} COMMAND ${test_target} )
endfunction()

add_library(picgen INTERFACE)
target_include_directories(picgen INTERFACE ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_BINARY_DIR ${CMAKE_CURRENT_LIST_DIR}/.build)

add_subdirectory(../kernel ${CMAKE_BINARY_DIR}/kernel)

# add_custom_test()
