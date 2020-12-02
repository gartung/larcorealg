# ======================================================================
#
# Testing
#
# ======================================================================

cet_test( printBacktrace_test
  LIBRARIES
    larcorealg_CoreUtils
    cetlib_except::cetlib_except
  )

cet_test(UncopiableAndUnmovableClass_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT
  )


cet_test( ProviderPack_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( RealComparisons_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( NumericUtils_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( DereferenceIterator_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( DumpUtils_test
  LIBRARIES
    cetlib::cetlib
    ROOT::Core
    ROOT::Physics
  USE_BOOST_UNIT )

cet_test( SortByPointers_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( span_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( MetaUtils_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( ContainerMeta_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( DebugUtils_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( counter_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( zip_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( enumerate_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( values_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( get_elements_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

cet_test( operations_test
  LIBRARIES
    cetlib::cetlib
  USE_BOOST_UNIT )

# test quiet_Math_Functor.h
get_property(IncludeDirectories DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)

cet_test(MathFunctor_workaround HANDBUILT
  TEST_EXEC ./quiet_Math_Functor_test.sh
  TEST_ARGS ${IncludeDirectories}
  DATAFILES quiet_Math_Functor_test.sh
  )
