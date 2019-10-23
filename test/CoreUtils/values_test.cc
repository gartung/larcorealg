/**
 * @file   values_test.cc
 * @brief  Test of `util::values()`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 9, 2019
 * 
 */


// testing library
#include "larcorealg/CoreUtils/values.h"

// Boost libraries
#define BOOST_TEST_MODULE ( values_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// C/C++ libraries
#include <numeric> // std::iota()
#include <vector>
#include <map>
#include <utility> // std::as_const()
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
void test_values() {
  
  // prepare the data set
  constexpr std::size_t N = 7U;
  
  std::vector<float> data(N);
  std::iota(data.begin(), data.end(), 0.0F);
  
  //
  // static checks
  //
  // for a vector, we expect exactly the same vector object as return value:
  static_assert(std::is_same_v<
    decltype(util::values(data)),
    std::add_lvalue_reference_t<std::vector<float>>
    >
    );
  BOOST_CHECK_EQUAL(&util::values(data), &data);
  
  //
  // read-only access
  //
  {
    std::size_t i = 0U;
    for (auto&& value: util::values(std::as_const(data))) {
      BOOST_CHECK_EQUAL(value, data[i]);
      BOOST_CHECK_EQUAL(std::addressof(value), std::addressof(data[i]));
      ++i;
    } // for
    BOOST_CHECK_EQUAL(i, data.size());
  }
  
  //
  // read/write access
  //
  {
    std::size_t i = 0U;
    for (auto&& value: util::values(data)) {
      BOOST_CHECK_EQUAL(value, data[i]);
      BOOST_CHECK_EQUAL(std::addressof(value), std::addressof(data[i]));
      ++i;
    } // for
    BOOST_CHECK_EQUAL(i, data.size());
  }
  
} // test_values()


void test_const_values() {
  
  // prepare the data set
  constexpr std::size_t N = 7U;
  
  std::vector<float> data(N);
  std::iota(data.begin(), data.end(), 0.0F);
  
  //
  // static checks
  //
  // for a vector, we expect exactly the same vector object as return value:
  static_assert(std::is_same_v<
    decltype(util::const_values(data)),
    std::add_lvalue_reference_t<std::vector<float> const>
    >
    );
  BOOST_CHECK_EQUAL(&util::const_values(data), &data);
  
  //
  // read-only access
  //
  std::size_t i = 0U;
  for (auto&& value: util::const_values(data)) {
    static_assert(std::is_const_v<std::remove_reference_t<decltype(value)>>);
    BOOST_CHECK_EQUAL(value, data[i]);
    BOOST_CHECK_EQUAL(std::addressof(value), std::addressof(data[i]));
    ++i;
  } // for
  BOOST_CHECK_EQUAL(i, data.size());
  
} // test_const_values()


// -----------------------------------------------------------------------------
template<template <typename Key, typename Value, typename... Args> typename Map>
void test_values_map() {
  
  // prepare the data set
  constexpr std::size_t N = 7U;
  constexpr double factor = 5.0;
  
  std::vector<int> keys(N);
  std::iota(keys.begin(), keys.end(), 0U);
  Map<int, float> data;
  for (int key: keys) data.emplace(key, key * factor);
  
  std::size_t i = 0U;
  for (auto&& value: util::values(data)) {
    static_assert(std::is_same_v<decltype(value), float&>);
    BOOST_CHECK_EQUAL(value, data[i]);
    BOOST_CHECK_EQUAL(std::addressof(value), std::addressof(data[i]));
    ++i;
  } // for
  BOOST_CHECK_EQUAL(i, data.size());
  
} // test_values_map()


// -----------------------------------------------------------------------------
void test_values_documentation() {
  /*
   * The promise:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::map<int, float> data { { 1, 4.0F }, { 3, 12.0F }, { 2, 8.0F } };
   * std::vector<float> values;
   * 
   * for (float value: util::values(data)) 
   *   values.push_back(value);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will result in `values` vector being of size `3` and with values
   * `{ 4.0F, 8.0F, 12.0F }`
   */
  
  std::map<int, float> data { { 1, 4.0F }, { 3, 12.0F }, { 2, 8.0F } };
  std::vector<float> values;
  
  for (float value: util::values(data)) 
    values.push_back(value);
  
  // the test
  BOOST_CHECK_EQUAL(values.size(), 3U);
  BOOST_CHECK_EQUAL(values[0],  4.0F);
  BOOST_CHECK_EQUAL(values[1],  8.0F);
  BOOST_CHECK_EQUAL(values[2], 12.0F);
  
} // test_values_documentation()


// -----------------------------------------------------------------------------
// BEGIN Test cases  -----------------------------------------------------------
// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(values_testcase) {
  
  test_values();
  test_const_values();
  test_values_map<std::map>();
  
  test_values_documentation();
  
} // BOOST_AUTO_TEST_CASE(values_testcase)


// -----------------------------------------------------------------------------
// END Test cases  -------------------------------------------------------------
// -----------------------------------------------------------------------------
