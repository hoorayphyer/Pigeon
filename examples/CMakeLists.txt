cmake_minimum_required(VERSION 3.13)
project (Pulsar)
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


# if(LAB_ROOT_DIR)
#   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LAB_ROOT_DIR}/lib)
#   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${LAB_ROOT_DIR})
# else()
#   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
#   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
# endif()

set(CMAKE_CXX_FLAGS_DEBUG "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")

add_library(picgen INTERFACE)
target_include_directories(picgen INTERFACE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_SOURCE_DIR}/include)
set(CMAKE_BINARY_DIR ${CMAKE_CURRENT_LIST_DIR}/.build)

add_subdirectory(../kernel ${CMAKE_BINARY_DIR}/kernel)
