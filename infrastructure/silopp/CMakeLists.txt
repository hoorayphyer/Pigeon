set(SILO_ROOT "${CMAKE_CURRENT_LIST_DIR}/silo-llnl")
include(ExternalProject)

find_package(HDF5)

set(build_own_hdf5 False)
if ( ${HDF5_FOUND} AND (NOT ${HDF5_IS_PARALLEL}) )
  list(GET HDF5_INCLUDE_DIRS 0 HDF5_INC) # in case more than one dir is found
else()
  # build the serial hdf5 when the system doesn't provide one
  set(pkgver "1.10.6")
  set(hdf5_ser_dir "hdf5_serial")

  ExternalProject_Add(hdf5-serial
    PREFIX "${SILO_ROOT}/${hdf5_ser_dir}"
    INSTALL_DIR "${SILO_ROOT}/${hdf5_ser_dir}"
    URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${pkgver}/src/hdf5-${pkgver}.tar.bz2
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-static --disable-sharedlib-rpath --enable-build-mode=production --enable-hl --enable-cxx --with-pic --with-zlib --with-szlib
    STEP_TARGETS install
    )
  set(HDF5_INC "${SILO_ROOT}/${hdf5_ser_dir}/include")
  set(build_own_hdf5 True)
endif()

ExternalProject_Add(silo-llnl
  PREFIX "${SILO_ROOT}"
  INSTALL_DIR "${SILO_ROOT}"
  URL https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-shared=yes --enable-fortran=no --enable-static=no --enable-optimization=yes --enable-browser=no --with-hdf5=${HDF5_INC},${HDF5_INC}/../lib
  STEP_TARGETS install
  )

if ( ${build_own_hdf5})
  ExternalProject_Add_StepDependencies(silo-llnl configure hdf5-serial-install)
endif()


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
