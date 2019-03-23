add_library(aperture_ensemble STATIC "ensemble.cpp")
target_include_directories(aperture_ensemble PRIVATE ${CMAKE_SOURCE_DIR}/aperture)
target_link_libraries(aperture_ensemble PRIVATE mpipp_lib)

add_library(aperture_dynamic_balance STATIC "dynamic_balance.cpp")
target_include_directories(aperture_dynamic_balance PRIVATE ${CMAKE_SOURCE_DIR}/aperture)
target_link_libraries(aperture_dynamic_balance PRIVATE mpipp_lib aperture_ensemble)
