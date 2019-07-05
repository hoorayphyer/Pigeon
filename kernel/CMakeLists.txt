# enable colored cmake message.
# Copied from https://stackoverflow.com/questions/18968979/how-to-get-colorized-output-with-cmake
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(Yellow      "${Esc}[33m")
endif()

set(KERNEL_DIR ${CMAKE_CURRENT_LIST_DIR})

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
endif()

function(test_gen category name)
  set(test_target "test_${category}_${name}")
  add_executable( ${test_target} "tests/test_${name}.cpp")
  target_include_directories(${test_target} PRIVATE ${KERNEL_DIR})
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

function(msperf_gen category name)
  set(test_target "msperf_${category}_${name}")
  add_executable( ${test_target} "tests/msperf_${name}.cpp")
  target_include_directories(${test_target} PRIVATE ${KERNEL_DIR} ${CMAKE_SOURCE_DIR})
  target_compile_options( ${test_target} PRIVATE -O3 )

  set(link_testfw testfw)
  foreach( lib IN LISTS ARGN)
    target_link_libraries(${test_target} PRIVATE ${lib})
    if( ${lib} STREQUAL "mpipp" )
      set(link_testfw testfw_mpi)
    endif()
  endforeach()
  target_link_libraries( ${test_target} PRIVATE ${link_testfw} )
  target_link_libraries( ${test_target} PRIVATE timer logger )
endfunction()

if (${ENABLE_DEBUG})
  add_subdirectory(debug)
  link_libraries(debugger)
endif()

add_subdirectory(apt)
add_subdirectory(timer)

add_subdirectory(filesys)
add_subdirectory(mpipp)
add_subdirectory(silopp)
add_subdirectory(logger)
add_subdirectory(random)

add_subdirectory(manifold)

add_subdirectory(field)
add_subdirectory(msh)
add_subdirectory(particle)

add_subdirectory(dye)
add_subdirectory(io)
add_subdirectory(ckpt)

add_executable(pic "simulator.cpp")
target_include_directories(pic PRIVATE ${KERNEL_DIR} )
target_link_libraries(pic PRIVATE field_updater field_sync particle_updater particle_migration particle_forces particle_scattering dye mpipp picgen io filesys logger ckpt particle_array timer)

if(${ENABLE_TESTS})
  add_subdirectory(testfw)
endif()
