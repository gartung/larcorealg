/**
 * @file   geo_vectors_parse_test.cc
 * @brief  Test of `larcorealg/Geometry/geo_vectors_parse.h` utilities.
 * @author Gianluca Petrillo (petrillo@slac.standford.edu)
 * @date   November 8, 2019
 */


// Boost libraries
#define BOOST_TEST_MODULE ( geo_vectors_parse_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_parse.h"
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h"

// ROOT libraries
// #include "TVector2.h"
#include "TVector3.h"
// #include "TLorentzVector.h"

// C/C++ standard library
#include <string>
#include <string_view>
#include <stdexcept> // std::runtime_error


//------------------------------------------------------------------------------
template <typename Vector, typename String = std::string>
void parse3D_string_test() {
  
  using Vector_t = Vector;
  using String_t = String;
  
  Vector_t const V1 { 2.5, 5.0, 3.0 };
  
  String_t const V1str1 { "( 2.5, 5.0, 3.0 )" };
  auto const V1parsed1 = geo::vect::parse::parse<Vector_t>(V1str1);
  BOOST_CHECK_EQUAL(V1parsed1, V1);
  
  String_t const V1str2 { "( 2.5; 5.0; 3.0 )" };
  auto const V1parsed2 = geo::vect::parse::parse<Vector_t>(V1str2);
  BOOST_CHECK_EQUAL(V1parsed2, V1);
  
  String_t const V1str3 { " { 2.5 ; 5.0 ; 30e-1 } " };
  auto const V1parsed3 = geo::vect::parse::parse<Vector_t>(V1str3);
  BOOST_CHECK_EQUAL(V1parsed3, V1);
  
  String_t const V1str4 { "( 2.5; 5.0; 3.0 )" };
  auto const V1parsed4 = geo::vect::parse::parse<Vector_t>(V1str4);
  BOOST_CHECK_EQUAL(V1parsed4, V1);
  
  //
  // check errors
  // 
  String_t const VerrStr1 { "( 2.5, 5.0, 3.0 " };
  BOOST_CHECK_THROW(geo::vect::parse::parse<Vector_t>(VerrStr1), std::runtime_error);
  
  String_t const VerrStr2 { "( 2.5, 5.0, )" };
  BOOST_CHECK_THROW(geo::vect::parse::parse<Vector_t>(VerrStr2), std::runtime_error);
  
  String_t const VerrStr3 { "( 2.5, 5.0; 3.0 )" };
  BOOST_CHECK_THROW(geo::vect::parse::parse<Vector_t>(VerrStr3), std::runtime_error);
  
  String_t const VerrStr4 { "( 2.5, 5.0  3.0 )" };
  BOOST_CHECK_THROW(geo::vect::parse::parse<Vector_t>(VerrStr4), std::runtime_error);
  
  String_t const VerrStr5 { "( , 5.0,  3.0 )" };
  BOOST_CHECK_THROW(geo::vect::parse::parse<Vector_t>(VerrStr5), std::runtime_error);
  
  String_t const VerrStr6 { " 2.5 , 5.0,  3.0 " };
  BOOST_CHECK_THROW(geo::vect::parse::parse<Vector_t>(VerrStr6), std::runtime_error);
  
  String_t const VerrStr7 { "( 2.5 , 5.0,  3.0 }" };
  BOOST_CHECK_THROW(geo::vect::parse::parse<Vector_t>(VerrStr7), std::runtime_error);
  
  
} // parse3D_string_test()


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(parse3D_string_testcase) {
  
  parse3D_string_test<geo::Vector_t, std::string>();
  parse3D_string_test<geo::Point_t , std::string>();
  parse3D_string_test<TVector3     , std::string>();
  
} // BOOST_AUTO_TEST_CASE(parse3D_string_testcase)

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(parse3D_stringview_testcase) {
  
  parse3D_string_test<geo::Vector_t, std::string_view>();
  parse3D_string_test<geo::Point_t , std::string_view>();
  parse3D_string_test<TVector3     , std::string_view>();
  
} // BOOST_AUTO_TEST_CASE(parse3D_stringview_testcase)

//------------------------------------------------------------------------------


