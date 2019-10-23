/**
 * @file   larcorealg/test/CoreUtils/operations_test.cc
 * @brief  Unit test for `larcorealg/CoreUtils/operations.h`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 5, 2019
 * @see    `larcorealg/CoreUtils/operations.h`
 */


// Boost libraries
#define BOOST_TEST_MODULE ( operations_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "larcorealg/CoreUtils/operations.h"
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"

// C/C++ standard libraries
#include <vector>
#include <algorithm>
#include <memory> // std::unique_ptr<>
#include <numeric> // std::iota()



//------------------------------------------------------------------------------
void test_AddressTaker_documentation() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int*> ptrs(data.size());
   * std::transform
   *   (data.begin(), data.end(), ptrs.begin(), util::AddressTaker{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  
  std::vector<int*> ptrs(data.size());
  std::transform
    (data.begin(), data.end(), ptrs.begin(), util::AddressTaker{});
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the test
  BOOST_CHECK_EQUAL(ptrs.size(), data.size());
  for (auto&& [ value, ptr ]: util::zip(data, ptrs)) {
    
    BOOST_CHECK(ptr);
    BOOST_CHECK_EQUAL(*ptr, value);
    BOOST_CHECK_EQUAL(ptr, &value);
    
  } // for
  
  
} // test_AddressTaker_documentation()

  
//------------------------------------------------------------------------------
void test_takeAddress_documentation() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int*> ptrs(data.size());
   * std::transform
   *   (data.begin(), data.end(), ptrs.begin(), util::takeAddress());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  
  std::vector<int*> ptrs(data.size());
  std::transform
    (data.begin(), data.end(), ptrs.begin(), util::takeAddress());
  
  // the test
  BOOST_CHECK_EQUAL(ptrs.size(), data.size());
  for (auto&& [ value, ptr ]: util::zip(data, ptrs)) {
    
    BOOST_CHECK(ptr);
    BOOST_CHECK_EQUAL(*ptr, value);
    BOOST_CHECK_EQUAL(ptr, &value);
    
  } // for
  
  
} // test_takeAddress_documentation()

  
//------------------------------------------------------------------------------
void test_takeAddress() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  std::vector<int const*> dataPtr;
  std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr),
    util::takeAddress());
  
  for (auto&& [ value, ptr ]: util::zip(data, dataPtr)) {
    
    BOOST_CHECK(ptr);
    BOOST_CHECK_EQUAL(*ptr, value);
    BOOST_CHECK_EQUAL(ptr, &value);
    
  } // for
  
} // test_takeAddress()


//------------------------------------------------------------------------------
void test_takeAddress_whyBother() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  std::vector<int const*> dataPtr;
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // `std::addressof()` approach
  using addressof_t = int const*(*)(int const&);
  
  std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr),
    ((addressof_t) &std::addressof));
  
  for (auto&& [ value, ptr ]: util::zip(data, dataPtr)) {
    
    BOOST_CHECK(ptr);
    BOOST_CHECK_EQUAL(*ptr, value);
    BOOST_CHECK_EQUAL(ptr, &value);
    
  } // for
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // lambda approach
  dataPtr.clear();
  
  auto takeAddress = [](auto&& ref){ return std::addressof(ref); };
  
  std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtr),
    takeAddress);
  
  for (auto&& [ value, ptr ]: util::zip(data, dataPtr)) {
    
    BOOST_CHECK(ptr);
    BOOST_CHECK_EQUAL(*ptr, value);
    BOOST_CHECK_EQUAL(ptr, &value);
    
  } // for
  
} // test_takeAddress_whyBother()


//------------------------------------------------------------------------------
void test_Dereferencer_documentation() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  std::vector<int*> ptrs(data.size());
  std::transform
    (data.begin(), data.end(), ptrs.begin(), util::takeAddress());
  
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> values(ptrs.size());
   * std::transform
   *   (ptrs.cbegin(), ptrs.cend(), values.begin(), util::Dereferencer{});
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  
  std::vector<int> values(ptrs.size());
  std::transform
    (ptrs.cbegin(), ptrs.cend(), values.begin(), util::Dereferencer{});
  
  // the test
  BOOST_CHECK_EQUAL(values.size(), data.size());
  for (auto&& [ value, orig ]: util::zip(data, values)) {
    
    BOOST_CHECK_EQUAL(value, orig);
    
  } // for
  
} // test_Dereferencer_documentation()

  
//------------------------------------------------------------------------------
void test_dereference_documentation() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  std::vector<int*> ptrs(data.size());
  std::transform
    (data.begin(), data.end(), ptrs.begin(), util::takeAddress());
  
  /*
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> values(ptrs.size());
   * std::transform
   *   (ptrs.cbegin(), ptrs.cend(), values.begin(), util::dereference());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  
  std::vector<int> values(ptrs.size());
  std::transform
    (ptrs.cbegin(), ptrs.cend(), values.begin(), util::dereference());
  
  // the test
  BOOST_CHECK_EQUAL(values.size(), data.size());
  for (auto&& [ value, orig ]: util::zip(data, values)) {
    
    BOOST_CHECK_EQUAL(value, orig);
    
  } // for
  
  
} // test_dereference_documentation()


//------------------------------------------------------------------------------
void test_dereference_C_ptr() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  std::vector<int const*> dataPtrs;
  std::transform(data.cbegin(), data.cend(), std::back_inserter(dataPtrs),
    util::takeAddress());
  
  std::vector<int> dataAgain;
  std::transform(
    dataPtrs.cbegin(), dataPtrs.cend(),
    std::back_inserter(dataAgain),
    util::dereference()
    );
  
  BOOST_CHECK_EQUAL(dataAgain.size(), data.size());
  for (auto&& [ value, valueAgain ]: util::zip(data, dataAgain)) {
    
    BOOST_CHECK_EQUAL(valueAgain, value);
    BOOST_CHECK_NE(&valueAgain, &value);
    
  } // for
  
} // test_dereference_C_ptr()


//------------------------------------------------------------------------------
void test_dereference_unique_ptr() {
  
  std::vector<int> data(10U);
  std::iota(data.begin(), data.end(), 0U);
  
  std::vector<std::unique_ptr<int>> dataPtrs;
  dataPtrs.reserve(data.size());
  for (auto&& value: data) dataPtrs.push_back(std::make_unique<int>(value)); 
  
  std::vector<int> dataAgain;
  std::transform(
    dataPtrs.cbegin(), dataPtrs.cend(),
    std::back_inserter(dataAgain),
    util::dereference()
    );
  
  BOOST_CHECK_EQUAL(dataAgain.size(), data.size());
  for (auto&& [ value, valueAgain ]: util::zip(data, dataAgain)) {
    
    BOOST_CHECK_EQUAL(valueAgain, value);
    BOOST_CHECK_NE(&valueAgain, &value);
    
  } // for
  
} // test_dereference_unique_ptr()


//------------------------------------------------------------------------------
void test_dereference_uncopiable() {
  
  struct ToughInt: private lar::UncopiableAndUnmovableClass {
    int value = 0;
  }; // ToughInt
  
  ToughInt value;
  ToughInt const* pValue = &value;
  
  BOOST_CHECK_EQUAL(pValue, &(util::dereference()(pValue)));
  
} // test_dereference_uncopiable()


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(takeAddress_testcase) {

  test_takeAddress();
  test_takeAddress_documentation();
  test_AddressTaker_documentation();
  test_takeAddress_whyBother();

} // BOOST_AUTO_TEST_CASE(takeAddress_testcase)


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(dereference_testcase) {

  test_dereference_C_ptr();
  test_dereference_unique_ptr();
  test_dereference_uncopiable();
  test_dereference_documentation();
  test_Dereferencer_documentation();

} // BOOST_AUTO_TEST_CASE(dereference_testcase)


//------------------------------------------------------------------------------
