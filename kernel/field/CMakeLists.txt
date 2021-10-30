add_library( field_sync SHARED sync.cpp )
target_include_directories( field_sync PUBLIC ${KERNEL_DIR} ${CMAKE_SOURCE_DIR} )
target_link_libraries( field_sync PRIVATE pic_rc mpi )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/log_spherical_solver)

if (${ENABLE_TESTS})
  test_gen(field field) # DONE
  test_gen(field mesh) # DONE
  test_gen(field sync mpi) # DONE
  # FIXME
  # msperf_gen(field updater mpi old_field_solver field_sync) # TODO old field solver can only work with Scalar with is float now
endif()
