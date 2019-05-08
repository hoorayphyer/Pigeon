add_library(old_field_solver STATIC)
target_include_directories(old_field_solver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${KERNEL_DIR}/src")
target_sources(old_field_solver PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/old_field_solver_adapter.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/Fields.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/FieldUpdater.cpp"
  "${CMAKE_CURRENT_LIST_DIR}/FiniteDiff.cpp"
  "${KERNEL_DIR}/simulator/field_communication.cpp")
target_link_libraries(old_field_solver PRIVATE mpipp picgen)
