add_library(old_field_solver SHARED "updater.cpp" "Fields.cpp" "FieldUpdater.cpp" "FiniteDiff.cpp")
target_include_directories(old_field_solver PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_SOURCE_DIR} ${KERNEL_DIR})
target_link_libraries(old_field_solver PRIVATE field_sync)
