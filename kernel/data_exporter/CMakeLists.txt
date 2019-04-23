add_library(dataexport SHARED "data_export.cpp")
target_include_directories(dataexport PRIVATE ${KERNEL_DIR} ${KERNEL_DIR}/include )
target_link_libraries( dataexport PRIVATE silopp dye picgen )
