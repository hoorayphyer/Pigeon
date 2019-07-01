add_library( debugger SHARED "debugger.cpp" )
target_include_directories( debugger PRIVATE ${KERNEL_DIR} )
