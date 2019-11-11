add_library( logger SHARED "logger.cpp" )
target_include_directories( logger PRIVATE ${PROJECT_SOURCE_DIR} )
