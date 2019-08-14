add_library(cartesian_field_updater SHARED "updater.cpp")
target_link_libraries(cartesian_field_updater PRIVATE field_sync picgen)
