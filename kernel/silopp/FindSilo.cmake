# (Slightly adapted from S. Johnson)
# - Find LLNL's Silo library
# This module defines
# Silo_INCLUDE_DIR, where to find blitz/blitz.h, etc.
# Silo_LIBRARIES, libraries to link against to use Silo.
# Silo_FOUND, If false, do not try to use Silo.
# also defined, but not for general use are
# Silo_LIBRARY, where to find the Silo library.
# The user should specify the head Silo director, Silo_DIR,
# or Silo_INC and Silo_LIB.

#if (Silo_DIR)
# message( "-- Using Silo_DIR: ${Silo_DIR} " )
#endif()
#if (Silo_LIB)
# message( "-- Using Silo_LIB: ${Silo_LIB} " )
#endif()
#if (Silo_INC)
# message( "-- Using Silo_INC: ${Silo_INC} " )
#endif()

# Set the cmake library path. This is because a clean
# build of cmake 2.8.4 seems to exclude this path. It might
# be a Ubuntu-specific quirk.
set(CMAKE_LIBRARY_PATH "/usr/lib/x86_64-linux-gnu")

find_path(Silo_INCLUDE_DIR
    NAMES silo.h
    PATHS ${Silo_DIR}/include
          ${Silo_INC}
)

find_library(Silo_LIBRARY
    NAMES siloh5 silo siloxx
    PATHS ${Silo_DIR}/lib
          ${Silo_LIB}
)

if (Silo_LIBRARY MATCHES "siloh5")
    MESSAGE( "-- Note: This Silo installation needs HDF5." )
    FIND_PACKAGE( HDF5 REQUIRED )
endif()

# Requires ZLib
find_package(ZLIB REQUIRED)

# handle the QUIETLY and REQUIRED arguments and set Silo_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Silo
    "Silo could not be found. Try setting Silo_DIR or Silo_LIB and Silo_INC."
    Silo_LIBRARY Silo_INCLUDE_DIR ZLIB_LIBRARIES HDF5_LIBRARIES)

if(SILO_FOUND)
    mark_as_advanced(Silo_INCLUDE_DIR Silo_LIBRARY HDF5_LIBRARIES ZLIB_LIBRARY)
    set( Silo_LIBRARIES ${Silo_LIBRARY} ${HDF5_LIBRARIES} ${ZLIB_LIBRARIES})
endif(SILO_FOUND)
