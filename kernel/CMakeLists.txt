# enable colored cmake message.
# Copied from https://stackoverflow.com/questions/18968979/how-to-get-colorized-output-with-cmake
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(Yellow      "${Esc}[33m")
endif()

set(KERNEL_DIR ${SHADOW_DIR} ${CMAKE_CURRENT_LIST_DIR})

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
endif()

# FIXME use imported library
include_directories(${CMAKE_CURRENT_LIST_DIR}/../infrastructure)

add_library(mpi INTERFACE)
# FIXME not sure why the RPATH info in mpipp is not transitive here. Build? Install?
find_package(MPI)
target_compile_options(mpi INTERFACE "${MPI_CXX_COMPILE_OPTIONS}")
target_include_directories(mpi INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../infrastructure ${MPI_CXX_INCLUDE_DIRS})
# NOTE somehow the following two lines should be commented out for things to work
# separate_arguments(MPI_CXX_LINK_OPTIONS NATIVE_COMMAND ${MPI_CXX_LINK_FLAGS}) # a workaround since findMPI doesn't provide this yet
# target_link_options(mpi INTERFACE ${MPI_CXX_LINK_OPTIONS})
target_link_directories(mpi INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../infrastructure/lib)
target_link_libraries(mpi INTERFACE ${MPI_CXX_LIBRARIES})
target_link_libraries(mpi INTERFACE -lmpipp)

add_library(testfw INTERFACE)
target_include_directories(testfw INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../.subtrees/Catch2/single_include)
target_link_directories(testfw INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../infrastructure/lib)
target_link_libraries(testfw INTERFACE -ltestfw)

add_library(testfw_mpi INTERFACE)
target_include_directories(testfw_mpi INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../.subtrees/Catch2/single_include)
target_link_directories(testfw_mpi INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../infrastructure/lib)
target_link_libraries(testfw_mpi INTERFACE -ltestfw_mpi)

add_library(logger INTERFACE)
target_link_directories(logger INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../infrastructure/lib)
target_link_libraries(logger INTERFACE -llogger)

function(test_gen category name)
  set(test_target "test_${category}_${name}")
  add_executable( ${test_target} "tests/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE ${KERNEL_DIR})
  target_compile_options( ${test_target} PRIVATE -O0 ) # -Wmissing-braces

  set(link_testfw testfw)
  foreach( lib IN LISTS ARGN)
    target_link_libraries(${test_target} PRIVATE ${lib})
    if( ${lib} STREQUAL "mpi" )
      set(link_testfw testfw_mpi mpi)
    endif()
  endforeach()
  target_link_libraries( ${test_target} PRIVATE ${link_testfw} )
  add_test( NAME TEST_${category}_${name} COMMAND ${test_target} )
endfunction()

function(msperf_gen category name)
  set(test_target "msperf_${category}_${name}")
  add_executable( ${test_target} "tests/msperf_${name}.cpp")
  target_include_directories(${test_target} PRIVATE ${KERNEL_DIR} ${CMAKE_SOURCE_DIR})
  target_compile_options( ${test_target} PRIVATE -O3 )

  set(link_testfw testfw)
  foreach( lib IN LISTS ARGN)
    target_link_libraries(${test_target} PRIVATE ${lib})
    if( ${lib} STREQUAL "mpi" )
      set(link_testfw testfw_mpi mpi)
    endif()
  endforeach()
  target_link_libraries( ${test_target} PRIVATE ${link_testfw} )
  target_link_libraries( ${test_target} PRIVATE timer logger )
endfunction()

if (${ENABLE_DEBUG})
  add_subdirectory(debug)
  link_libraries(debugger logger)
endif()

add_subdirectory(field)
add_subdirectory(msh)
add_subdirectory(particle)

add_subdirectory(dye)
add_subdirectory(io)
add_subdirectory(ckpt)

add_subdirectory(simulator)
