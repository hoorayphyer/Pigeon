add_library(haugbolle_solver SHARED "updater.cpp")
target_link_libraries(haugbolle_solver PRIVATE field_sync pic_rc)
