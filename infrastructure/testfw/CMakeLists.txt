add_library(testfw SHARED "catch_main.cpp")
target_include_directories(testfw PUBLIC ${PROJECT_SOURCE_DIR}/../.subtrees/Catch2/single_include)

add_library(testfw_mpi SHARED "catch_main_mpi++.cpp")
target_include_directories(testfw_mpi PUBLIC ${PROJECT_SOURCE_DIR}/../.subtrees/Catch2/single_include )
target_link_libraries(testfw_mpi PUBLIC mpipp)
