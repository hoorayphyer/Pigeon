list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} )
find_package(Silo REQUIRED)
add_library(silopp SHARED "silo++.cpp" "silo_optlist.cpp" "silo_operations.cpp" "silo_navigation.cpp" "silo_reader.cpp" "silo_toc.cpp")
target_include_directories(silopp PUBLIC ${Silo_INCLUDE_DIR} ${KERNEL_DIR} )
target_link_libraries( silopp PUBLIC ${Silo_LIBRARIES} )

add_executable(dbtoc "dbtoc.cpp")
target_include_directories(dbtoc PRIVATE ${Silo_INCLUDE_DIR} ${KERNEL_DIR} )
target_link_libraries( dbtoc PRIVATE silopp )

if (${ENABLE_TESTS})
  test_gen(silopp silo++ mpipp silopp filesys)
endif()
