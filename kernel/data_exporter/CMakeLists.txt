add_library(data_exporter SHARED "data_export.cpp")
target_include_directories(data_exporter PRIVATE ${KERNEL_DIR} ${KERNEL_DIR}/include )
target_link_libraries( data_exporter PRIVATE silopp dye picgen filesys field_communication )
