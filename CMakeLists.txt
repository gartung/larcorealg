# ======================================================================
#  larcorealg main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larcorealg/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================
# use cmake 3.14 or later
cmake_minimum_required(VERSION 3.18.2 FATAL_ERROR)
cmake_policy(VERSION 3.14)
# Policy CMP0096 must be set to NEW to preserve leading zeros in the
# VERSION set by project().
cmake_policy(SET CMP0096 NEW)



find_package(cetmodules 2.05.03 REQUIRED)
project(larcorealg VERSION 09.00.01 LANGUAGES CXX)

# cetbuildtools contains our cmake modules ### MIGRATE-ACTION-REQUIRED: remove
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} ) ### MIGRATE-ACTION-REQUIRED: remove
IF (NOT CETBUILDTOOLS_VERSION) ### MIGRATE-ACTION-REQUIRED: remove
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules") ### MIGRATE-ACTION-REQUIRED: remove
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules ### MIGRATE-ACTION-REQUIRED: remove
		      $ENV{CETBUILDTOOLS_DIR}/Modules ### MIGRATE-ACTION-REQUIRED: remove
		      ${CMAKE_MODULE_PATH}) ### MIGRATE-ACTION-REQUIRED: remove

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

# these are minimum required versions, not the actual product versions
find_ups_product( canvas ) ### MIGRATE-ACTION-RECOMMENDED: use cet_find_package() with target semantics for linking
find_ups_product( larcoreobj ) ### MIGRATE-ACTION-RECOMMENDED: use cet_find_package() with target semantics for linking
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

# packaging utility
cet_cmake_config(COMPATIBILITY SameMajorVersion)
