add_library(data_exporter SHARED "data_export.cpp")
target_include_directories(data_exporter PRIVATE ${KERNEL_DIR} ${KERNEL_DIR}/include )
target_compile_definitions(data_exporter PRIVATE APPARENT_DATA_DIR="${CMAKE_SOURCE_DIR}/Data")
target_link_libraries( data_exporter PRIVATE silopp dye picgen filesys field_communication )

add_executable( test_data_export test_data_export.cpp )
target_include_directories( test_data_export PRIVATE ${KERNEL_DIR} ${KERNEL_DIR}/tests )
target_link_libraries( test_data_export PRIVATE silopp catch_main_mpipp mpipp filesys )
target_link_libraries( test_data_export PRIVATE stdc++fs) #TODO this is only for gnu
set_target_properties( test_data_export
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin"
  )

add_test( NAME TEST_DATA_EXPORT COMMAND test_data_export )
