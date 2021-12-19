set(SILO_ROOT "${CMAKE_CURRENT_LIST_DIR}/silo-llnl")
include(ExternalProject)

find_package(HDF5)

# silo (as latest as silo-4.11, no patch number) doesn't work with HDF5 higher than 1.10.7.
set(hdf5_max_ver "1.10.7")

set(build_own_hdf5 False)
if (${HDF5_FOUND} AND (NOT ${HDF5_IS_PARALLEL}) AND (${HDF5_VERSION} VERSION_LESS_EQUAL ${hdf5_max_ver}) )
  list(GET HDF5_INCLUDE_DIRS 0 HDF5_INC) # in case more than one dir is found
else()
  # build the serial hdf5 when the system doesn't provide one
  set(hdf5_ser_dir "hdf5_serial")

  # extract the major and minor parts of the version
  string(REPLACE "." ";" version_list ${hdf5_max_ver})
  list(GET version_list 0 hdf5_max_ver_major)
  list(GET version_list 1 hdf5_max_ver_minor)

  ExternalProject_Add(hdf5-serial
    PREFIX "${SILO_ROOT}/${hdf5_ser_dir}"
    INSTALL_DIR "${SILO_ROOT}/${hdf5_ser_dir}"
    URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_max_ver_major}.${hdf5_max_ver_minor}/hdf5-${hdf5_max_ver}/src/hdf5-${hdf5_max_ver}.tar.bz2
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-static --disable-sharedlib-rpath --enable-build-mode=production --enable-hl --enable-cxx --with-pic --with-zlib --with-szlib
    STEP_TARGETS install
    )
  set(HDF5_INC "${SILO_ROOT}/${hdf5_ser_dir}/include")
  set(build_own_hdf5 True)
endif()

ExternalProject_Add(silo-llnl
  PREFIX "${SILO_ROOT}"
  INSTALL_DIR "${SILO_ROOT}"
  URL https://wci.llnl.gov/sites/wci/files/2021-09/silo-4.11.tgz
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

## FIXME may have problem with parallel hdf5?
add_executable(combine_silos "combine_silos.cpp")
target_include_directories(combine_silos PRIVATE "${SILO_ROOT}/include" ${PROJECT_SOURCE_DIR}) ## FIXME silo not needed, use correct hdf5
target_link_libraries( combine_silos PRIVATE filesys ${HDF5_LIBRARIES})

if (${ENABLE_TESTS})
  test_gen(silopp silo++ mpipp silopp filesys)
endif()
