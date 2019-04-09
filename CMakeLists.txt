# ======================================================================
#  larcorealg main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larcorealg/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 3.10 or later
cmake_minimum_required (VERSION 3.10)

project(larcorealg VERSION 08.05.00 LANGUAGES CXX)
set(fcl_dir job)
set(gdml_dir gdml)
find_package(cetmodules REQUIRED)
include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_package(Boost REQUIRED COMPONENTS date_time filesystem graph program_options
  regex system unit_test_framework)
find_package( larcoreobj REQUIRED)

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(CetMake)
include(BasicPlugin)

# source
add_subdirectory(larcorealg)

# tests 
add_subdirectory(test)
add_subdirectory(cmakeConfig)
