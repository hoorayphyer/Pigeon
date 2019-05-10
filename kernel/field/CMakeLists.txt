add_library( field_communication SHARED communication.cpp )
target_include_directories( field_communication PUBLIC ${KERNEL_DIR} ${CMAKE_SOURCE_DIR} )
target_link_libraries( field_communication PRIVATE picgen mpipp )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/old_field_solver)

add_library(field_updater SHARED "updater.cpp")
target_link_libraries(field_updater PRIVATE field_communication picgen old_field_solver)

test_gen(field field) # DONE
test_gen(field mesh) # DONE
test_gen(field communication mpipp) # DONE
