cmake_minimum_required(VERSION 3.13)
project () # specify your project name here
option(ENABLE_TESTS "Enable tests" ON)
option(ENABLE_DEBUG "Enable debug" OFF)

# message ("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")
if (${ENABLE_TESTS})
    enable_testing()
endif()

if (${ENABLE_DEBUG})
  add_compile_definitions(PIC_DEBUG)
endif()

# Set CXX flags
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC ${CXX_EXTRA_FLAGS}")
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  string(APPEND CMAKE_CXX_FLAGS_RELEASE "-xHost -no-prec-div")
else()
  string(APPEND CMAKE_CXX_FLAGS_RELEASE "-march=native")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
# set(CMAKE_CXX_COMPILER "g++")

function(add_custom_test name)
  set(test_target "test_custom_${name}")
  add_executable( ${test_target} "./tests/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE "../kernel" ${CMAKE_CURRENT_LIST_DIR} )
  target_compile_options( ${test_target} PIRVATE -O0 ) # -Wmissing-braces
  set_target_properties( ${test_target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

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

if (${ENABLE_TESTS})
  # add_custom_test()
endif()
