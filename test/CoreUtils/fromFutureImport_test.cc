/**
 * @file   fromFutureImport_test.cc
 * @brief  Tests some features in `larcorealg/CoreUtils/fromFutureImport_test.h`
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   April 28, 2016
 * @see    `larcorealg/CoreUtils/fromFutureImport_test.h
 */

/*
 * Boost Magic: define the name of the module;
 * and do that before the inclusion of Boost unit test headers
 * because it will change what they provide.
 * Among the those, there is a main() function and some wrapping catching
 * unhandled exceptions and considering them test failures, and probably more.
 */
#define BOOST_TEST_MODULE ( fromFutureImport_test )

// LArSoft libraries
#include "larcorealg/CoreUtils/fromFutureImport.h"

// Boost libraries
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK_THROW()...

// C/C++ standard libraries
#include <string_view>
#include <string>
#include <stdexcept> // std::logic_error
#include <type_traits>


//------------------------------------------------------------------------------
void string_view_concatenation_test() {
  using namespace std::string_literals;
  using namespace std::string_view_literals;
  
  using namespace util::not_std::string_view_ops;
  
  static_assert(std::is_same_v<decltype("<A>"s + "<B>"sv), std::string>);
  static_assert(std::is_same_v<decltype("<A>"sv + "<B>"s), std::string>);
  
  BOOST_CHECK_EQUAL("<A>"s + "<B>"sv, "<A><B>");
  BOOST_CHECK_EQUAL("<A>"sv + "<B>"s, "<A><B>");
} // string_view_concatenation_test()


//------------------------------------------------------------------------------
void checkLength(std::string_view const& s, std::size_t const maxLength) {
  using namespace std::string_literals;
  using namespace util::not_std::string_view_ops;
  
  if (s.length() <= maxLength) return;
  
  throw std::logic_error(
    "String '"s + s + "' is "
    + std::to_string(s.length()) + " characters long!"
    );
  
} // void checkLength()

void string_view_concatenation_documentation_test() {
  using namespace std::string_view_literals;
  
  /*
   * The promise is `checkLength()` above, compiling; we also run it.
   */
  
  BOOST_CHECK_NO_THROW(checkLength("test a string"sv, 30U));
  
  BOOST_CHECK_THROW(checkLength("test a string"sv, 10U), std::logic_error);
  
} // string_view_concatenation_documentation_test()


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(string_view_concatenation_testcase) {
  
  string_view_concatenation_test();
  string_view_concatenation_documentation_test();
  
} // BOOST_AUTO_TEST_CASE(string_view_concatenation_testcase)


//------------------------------------------------------------------------------
