add_library( field_sync SHARED sync.cpp )
target_include_directories( field_sync PUBLIC ${KERNEL_DIR} ${CMAKE_SOURCE_DIR} )
target_link_libraries( field_sync PRIVATE picgen mpipp )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/old_field_solver)

add_library(field_updater SHARED "updater.cpp")
target_link_libraries(field_updater PRIVATE field_sync picgen old_field_solver)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cartesian_solver)

if (${ENABLE_TESTS})
  test_gen(field field) # DONE
  test_gen(field mesh) # DONE
  test_gen(field sync mpipp) # DONE
  msperf_gen(field updater mpipp old_field_solver field_sync) # TODO old field solver can only work with Scalar with is float now
endif()
