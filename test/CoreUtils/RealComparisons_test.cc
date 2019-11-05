/**
 * @file    RealComparisons_test.cc
 * @brief   Unit test for RealComparisons class
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    November 4th, 2016
 * @see     RealComparisons.h
 */

// Boost libraries
#define BOOST_TEST_MODULE ( RealComparisons_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"

// C/C++ standard libraries
#include <ostream>
#include <cmath>


BOOST_AUTO_TEST_CASE(test_RealComparisons) {

  // instantiate a RealComparisons with two classes
  lar::util::RealComparisons<float> check(1e-5);

  double const sqrt2 = std::sqrt(2);

  // check zero()
  double const epsilon = 2.0 - (sqrt2 * sqrt2);
  BOOST_CHECK( check.zero(epsilon));
  BOOST_CHECK( check.zero(0.0));
  BOOST_CHECK( check.zero(1e-5));
  BOOST_CHECK( check.zero(-1e-5));
  BOOST_CHECK(!check.zero(1.01e-5));
  BOOST_CHECK(!check.zero(-1.01e-5 * 1.01));

  // check nonZero()
  BOOST_CHECK(!check.nonZero(epsilon));
  BOOST_CHECK(!check.nonZero(0.0));
  BOOST_CHECK(!check.nonZero(1e-5));
  BOOST_CHECK(!check.nonZero(-1e-5));
  BOOST_CHECK( check.nonZero(1.01e-5));
  BOOST_CHECK( check.nonZero(-1.01e-5 * 1.01));

  // check equal()
  BOOST_CHECK(!check.equal(sqrt2, 1.4142));
  BOOST_CHECK( check.nonEqual(sqrt2, 1.4142));
  BOOST_CHECK( check.equal(sqrt2, 1.414213));
  BOOST_CHECK(!check.nonEqual(sqrt2, 1.414213));

  // check strictlyNegative()
  BOOST_CHECK(!check.strictlyNegative(+1e-5 + 1e-7)); // outside tolerance
  BOOST_CHECK(!check.strictlyNegative(+1e-5 - 1e-7)); // within tolerance
  BOOST_CHECK(!check.strictlyNegative(0.0));          // zero
  BOOST_CHECK(!check.strictlyNegative(-1e-5 + 1e-7)); // within tolerance
  BOOST_CHECK( check.strictlyNegative(-1e-5 - 1e-7)); // outside tolerance

  // check strictlyPositive()
  BOOST_CHECK( check.strictlyPositive(+1e-5 + 1e-7)); // outside tolerance
  BOOST_CHECK(!check.strictlyPositive(+1e-5 - 1e-7)); // within tolerance
  BOOST_CHECK(!check.strictlyPositive(0.0));          // zero
  BOOST_CHECK(!check.strictlyPositive(-1e-5 + 1e-7)); // within tolerance
  BOOST_CHECK(!check.strictlyPositive(-1e-5 - 1e-7)); // outside tolerance

  // check nonNegative()
  BOOST_CHECK( check.nonNegative(+1e-5 + 1e-7)); // outside tolerance
  BOOST_CHECK( check.nonNegative(+1e-5 - 1e-7)); // within tolerance
  BOOST_CHECK( check.nonNegative(0.0));          // zero
  BOOST_CHECK( check.nonNegative(-1e-5 + 1e-7)); // within tolerance
  BOOST_CHECK(!check.nonNegative(-1e-5 - 1e-7)); // outside tolerance

  // check nonPositive()
  BOOST_CHECK(!check.nonPositive(+1e-5 + 1e-7)); // outside tolerance
  BOOST_CHECK( check.nonPositive(+1e-5 - 1e-7)); // within tolerance
  BOOST_CHECK( check.nonPositive(0.0));          // zero
  BOOST_CHECK( check.nonPositive(-1e-5 + 1e-7)); // within tolerance
  BOOST_CHECK( check.nonPositive(-1e-5 - 1e-7)); // outside tolerance

  // check strictlySmaller()
  BOOST_CHECK(!check.strictlySmaller(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_CHECK(!check.strictlySmaller(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_CHECK(!check.strictlySmaller(1.0, 1.0));        // equal
  BOOST_CHECK(!check.strictlySmaller(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_CHECK( check.strictlySmaller(1.0, 1.0 + 1e-4)); // outside tolerance

  // check nonSmaller()
  BOOST_CHECK( check.nonSmaller(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_CHECK( check.nonSmaller(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_CHECK( check.nonSmaller(1.0, 1.0));        // equal
  BOOST_CHECK( check.nonSmaller(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_CHECK(!check.nonSmaller(1.0, 1.0 + 1e-4)); // outside tolerance

  // check strictlyGreater()
  BOOST_CHECK( check.strictlyGreater(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_CHECK(!check.strictlyGreater(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_CHECK(!check.strictlyGreater(1.0, 1.0));        // equal
  BOOST_CHECK(!check.strictlyGreater(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_CHECK(!check.strictlyGreater(1.0, 1.0 + 1e-4)); // outside tolerance

  // check nonGreater()
  BOOST_CHECK(!check.nonGreater(1.0, 1.0 - 1e-4)); // outside tolerance
  BOOST_CHECK( check.nonGreater(1.0, 1.0 - 1e-7)); // within tolerance
  BOOST_CHECK( check.nonGreater(1.0, 1.0));        // equal
  BOOST_CHECK( check.nonGreater(1.0, 1.0 + 1e-7)); // within tolerance
  BOOST_CHECK( check.nonGreater(1.0, 1.0 + 1e-4)); // outside tolerance

  // check within()
  BOOST_CHECK(!check.within(sqrt2, 0., 1.41420));
  BOOST_CHECK( check.within(sqrt2, 0., 1.41421));
  BOOST_CHECK( check.within(sqrt2, 1.41422, 2.));
  BOOST_CHECK(!check.within(sqrt2, 1.41423, 2.));

  // check inverted limits
  BOOST_CHECK(!check.within(sqrt2, 1.41421, 0.));
  BOOST_CHECK(check.withinSorted(sqrt2, 1.41421, 0.));

} // BOOST_AUTO_TEST_CASE(test_RealComparisons)
