# enable colored cmake message.
# Copied from https://stackoverflow.com/questions/18968979/how-to-get-colorized-output-with-cmake
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset "${Esc}[m")
  set(Yellow      "${Esc}[33m")
endif()

set(KERNEL_DIR ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(src)
add_subdirectory(dye)
add_subdirectory(picgen)
add_subdirectory(tests)
