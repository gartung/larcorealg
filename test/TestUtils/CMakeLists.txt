
cet_test(NameSelector_test
  SOURCE NameSelector_test.cxx
  LIBRARIES
    larcorealg_TestUtils
    canvas::canvas
  USE_BOOST_UNIT
)

cet_test(ProviderList_test
  LIBRARIES
    canvas::canvas
  USE_BOOST_UNIT
)

cet_test(StopWatch_test)
