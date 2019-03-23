add_library(old_field_solver
  STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/old_field_solver_adapter.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Fields.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/FieldUpdater.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/FiniteDiff.cpp")
target_include_directories(old_field_solver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(old_field_solver PRIVATE mpipp_lib)
