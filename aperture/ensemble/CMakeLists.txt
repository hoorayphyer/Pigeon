add_library(aperture_dynamic_ensemble STATIC "ensemble.cpp" "dynamic_balance.cpp")
target_include_directories(aperture_dynamic_ensemble PRIVATE ${CMAKE_SOURCE_DIR}/aperture)
target_link_libraries(aperture_dynamic_ensemble PRIVATE mpipp_lib)
