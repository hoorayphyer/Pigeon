add_library(old_field_solver STATIC)
target_include_directories(old_field_solver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src")
target_sources(old_field_solver PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/old_field_solver_adapter.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/Fields.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/FieldUpdater.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/FiniteDiff.cpp"
  "${CMAKE_SOURCE_DIR}/picgen/field_communication.cpp")
target_link_libraries(old_field_solver PRIVATE mpipp picgen)
