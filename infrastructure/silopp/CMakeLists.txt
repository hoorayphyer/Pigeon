find_package(HDF5 REQUIRED)

list(GET HDF5_INCLUDE_DIRS 0 HDF5_INC) # in case more than one dir is found

set(SILO_ROOT "${CMAKE_CURRENT_LIST_DIR}/silo-llnl")
include(ExternalProject)
ExternalProject_Add(silo-llnl
  PREFIX "${SILO_ROOT}"
  INSTALL_DIR "${SILO_ROOT}"
  URL https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-shared=yes --enable-fortran=no --enable-static=no --enable-optimization=yes --enable-browser=no --with-hdf5=${HDF5_INC},${HDF5_INC}/../lib
  STEP_TARGETS install
  )

add_library(silopp SHARED "silo++.cpp" "silo_optlist.cpp" "silo_operations.cpp" "silo_navigation.cpp" "silo_reader.cpp" "silo_toc.cpp")
add_dependencies(silopp silo-llnl-install)
target_include_directories(silopp PUBLIC "${SILO_ROOT}/include" ${PROJECT_SOURCE_DIR})
target_link_directories( silopp PUBLIC "${SILO_ROOT}/lib" )
target_link_libraries( silopp PUBLIC siloh5 )

add_executable(dbtoc "dbtoc.cpp")
target_include_directories(dbtoc PRIVATE "${SILO_ROOT}/include" ${PROJECT_SOURCE_DIR})
target_link_libraries( dbtoc PRIVATE silopp )

if (${ENABLE_TESTS})
  test_gen(silopp silo++ mpipp silopp filesys)
endif()
