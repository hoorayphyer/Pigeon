# TODO use INTERFACE library
add_library(dye STATIC "ensemble.cpp" "dynamic_balance.cpp")
target_include_directories(dye PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(dye PRIVATE mpipp)
