add_library(log_spherical_solver SHARED "updater.cpp")
target_link_libraries(log_spherical_solver PRIVATE field_sync pic_rc)

if (${ENABLE_TESTS})
  test_gen(field log_spherical_solver)
endif()

test_gen(field log_spherical_solver) # FIXME
