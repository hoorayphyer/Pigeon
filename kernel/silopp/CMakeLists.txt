list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} )
find_package(Silo REQUIRED)
add_library(silopp SHARED "silo++.cpp" "silo_optlist.cpp" "silo_operations.cpp")
target_include_directories(silopp PUBLIC ${Silo_INCLUDE_DIR} ${KERNEL_DIR} )
target_link_libraries( silopp PUBLIC ${Silo_LIBRARIES} )

test_gen(silopp silo++ mpipp silopp filesys)
