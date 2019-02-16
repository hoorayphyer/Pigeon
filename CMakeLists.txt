cmake_minimum_required(VERSION 2.8)
project (Aperture)
message ("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

# enable colored cmake message.
# Copied from https://stackoverflow.com/questions/18968979/how-to-get-colorized-output-with-cmake
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(Yellow      "${Esc}[33m")
endif()

# Set CXX flags
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC -xHost -no-prec-div ${CXX_EXTRA_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -pthread -O3 -g -DNDEBUG -Wall -Wextra -fPIC -march=native ${CXX_EXTRA_FLAGS}")
endif()
# set(CMAKE_CXX_FLAGS_RELEASE "-pthread -O3 -DNDEBUG -Wall -Wno-unused-variable -std=gnu++14 -fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")
#set(CMAKE_CXX_FLAGS "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++14 -fPIC")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
# set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
# set(Project_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

# Enable testing
enable_testing(true)

# if(LAB_ROOT_DIR)
#   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LAB_ROOT_DIR}/lib)
#   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${LAB_ROOT_DIR})
# else()
#   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
#   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
# endif()

find_package(MPI)
# find_package(Silo REQUIRED)
set(CMAKE_CXX_FLAGS_DEBUG "-pthread -g3 -O0 -Wall -Wno-unused-variable -std=gnu++17 -fPIC")

add_subdirectory(src)
# add_subdirectory(tests)
# add_subdirectory(tools)
