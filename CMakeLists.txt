# ======================================================================
#  larcorealg main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larcorealg/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.12...3.14 FATAL_ERROR)

project(larcorealg LANGUAGES CXX)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( canvas )
find_ups_product( cetbuildtools )
find_ups_product( larcoreobj )
find_ups_root()

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(CetMake)
include(BasicPlugin)

# source
add_subdirectory(larcorealg)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
