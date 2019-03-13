add_library(old_field_solver
  STATIC
  "${CMAKE_SOURCE_DIR}/old_field_solver/old_field_solver_adapter.cpp"
  "${CMAKE_SOURCE_DIR}/old_field_solver/Fields.cpp"
  "${CMAKE_SOURCE_DIR}/old_field_solver/FieldUpdater.cpp"
  "${CMAKE_SOURCE_DIR}/old_field_solver/FiniteDiff.cpp")
target_include_directories(old_field_solver PRIVATE "${CMAKE_SOURCE_DIR}/old_field_solver")
link_mpipp(old_field_solver)
