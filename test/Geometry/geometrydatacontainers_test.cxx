/**
 * @file   geometrydatacontainers_test.cxx
 * @brief  Unit test for GeometryDataContainers.h library.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 2nd, 2018
 *
 *
 */

// LArSoft libraries
#include "larcorealg/Geometry/GeometryDataContainers.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Boost libraries
#include "cetlib/quiet_unit_test.hpp" // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK()
#include <boost/test/floating_point_comparison.hpp> // BOOST_CHECK_CLOSE()


//------------------------------------------------------------------------------
template <typename T>
struct Summer {
  
  T sum = T { 0 };
  
  void operator() (T v) { sum += v; }
  
  T get() const { return sum; }
  void reset() { sum = T{0}; }
  
}; // struct Summer

//------------------------------------------------------------------------------
void TPCDataContainerTest() {
  
  static constexpr std::size_t NCryostats = 2U;
  static constexpr std::size_t NTPCs      = 3U;
  static constexpr std::size_t N          = NCryostats * NTPCs;
  
  geo::TPCDataContainer<int> data(NCryostats, NTPCs);

  BOOST_CHECK(!data.empty());
  BOOST_CHECK_EQUAL(data.size(), N);
  BOOST_CHECK_GE(data.capacity(), N);
  
  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      BOOST_CHECK_EQUAL((data[{ c, t }]), 0);
  
  BOOST_CHECK_EQUAL(data.firstID(), geo::TPCID(0, 0));
  BOOST_CHECK_EQUAL(data.lastID(), geo::TPCID(1, 2));

  BOOST_CHECK( data.hasTPC({ 0,  0}));
  BOOST_CHECK( data.hasTPC({ 0,  1}));
  BOOST_CHECK( data.hasTPC({ 0,  2}));
  BOOST_CHECK(!data.hasTPC({ 0,  3}));
  BOOST_CHECK(!data.hasTPC({ 0,  4}));
  BOOST_CHECK( data.hasTPC({ 1,  0}));
  BOOST_CHECK( data.hasTPC({ 1,  1}));
  BOOST_CHECK( data.hasTPC({ 1,  2}));
  BOOST_CHECK(!data.hasTPC({ 1,  3}));
  BOOST_CHECK(!data.hasTPC({ 1,  4}));
  BOOST_CHECK(!data.hasTPC({ 2,  0}));
  BOOST_CHECK(!data.hasTPC({ 2,  1}));
  BOOST_CHECK(!data.hasTPC({ 2,  2}));
  BOOST_CHECK(!data.hasTPC({ 2,  3}));
  BOOST_CHECK(!data.hasTPC({ 2,  4}));

  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 0,  0}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 0,  1}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 0,  2}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 0,  3}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 0,  4}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 1,  0}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 1,  1}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 1,  2}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 1,  3}));
  BOOST_CHECK( data.hasCryostat(geo::TPCID{ 1,  4}));
  BOOST_CHECK(!data.hasCryostat(geo::TPCID{ 2,  0}));
  BOOST_CHECK(!data.hasCryostat(geo::TPCID{ 2,  1}));
  BOOST_CHECK(!data.hasCryostat(geo::TPCID{ 2,  2}));
  BOOST_CHECK(!data.hasCryostat(geo::TPCID{ 2,  3}));
  BOOST_CHECK(!data.hasCryostat(geo::TPCID{ 2,  4}));

  data[{0, 0}] = 4;
  BOOST_CHECK_EQUAL((data[{0, 0}]), 4);
  BOOST_CHECK_EQUAL(data.at({0, 0}), 4);
  data[{0, 0}] = 5;
  BOOST_CHECK_EQUAL((data[{0, 0}]), 5);
  BOOST_CHECK_EQUAL(data.at({0, 0}), 5);

  data[{0, 1}] = 6;
  BOOST_CHECK_EQUAL((data[{0, 1}]), 6);
  BOOST_CHECK_EQUAL(data.at({0, 1}), 6);

  BOOST_CHECK_EQUAL((data[{0, 0}]),  5);

  data[{0, 2}] = 7;
  BOOST_CHECK_EQUAL((data[{0, 2}]), 7);
  BOOST_CHECK_EQUAL(data.at({0, 2}), 7);

  BOOST_CHECK_EQUAL((data[{0, 0}]),  5);
  BOOST_CHECK_EQUAL((data[{0, 1}]),  6);

  data[{1, 0}] = 15;
  BOOST_CHECK_EQUAL((data[{1, 0}]), 15);
  BOOST_CHECK_EQUAL(data.at({1, 0}), 15);

  BOOST_CHECK_EQUAL((data[{0, 0}]),  5);
  BOOST_CHECK_EQUAL((data[{0, 1}]),  6);
  BOOST_CHECK_EQUAL((data[{0, 2}]),  7);

  data[{1, 1}] = 16;
  BOOST_CHECK_EQUAL((data[{1, 1}]), 16);
  BOOST_CHECK_EQUAL(data.at({1, 1}), 16);

  BOOST_CHECK_EQUAL((data[{0, 0}]),  5);
  BOOST_CHECK_EQUAL((data[{0, 1}]),  6);
  BOOST_CHECK_EQUAL((data[{0, 2}]),  7);
  BOOST_CHECK_EQUAL((data[{1, 0}]), 15);

  data[{1, 2}] = 17;
  BOOST_CHECK_EQUAL((data[{1, 2}]), 17);
  BOOST_CHECK_EQUAL(data.at({1, 2}), 17);

  BOOST_CHECK_EQUAL((data[{0, 0}]),  5);
  BOOST_CHECK_EQUAL((data[{0, 1}]),  6);
  BOOST_CHECK_EQUAL((data[{0, 2}]),  7);
  BOOST_CHECK_EQUAL((data[{1, 0}]), 15);
  BOOST_CHECK_EQUAL((data[{1, 1}]), 16);

  BOOST_CHECK_THROW(data.at({0, 3}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 4}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 4}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 4}), std::out_of_range);

  BOOST_CHECK_EQUAL(data.first(), 5);
  data.first() = -5;
  BOOST_CHECK_EQUAL((data[{0, 0}]), -5);
  BOOST_CHECK_EQUAL(data.first(), -5);
  data.first() =  5;

  BOOST_CHECK_EQUAL(data.last(), 17);
  data.last() = -17;
  BOOST_CHECK_EQUAL((data[{1U, 2U}]), -17);
  BOOST_CHECK_EQUAL(data.last(), -17);
  data.last() =  17;
  
  auto const& constData = data;

  BOOST_CHECK_EQUAL
    (std::addressof(constData.first()), std::addressof(data.first()));
  BOOST_CHECK_EQUAL
    (std::addressof(constData.last()), std::addressof(data.last()));

  BOOST_CHECK_EQUAL((constData[{0, 0}]), (data[{0, 0}]));
  BOOST_CHECK_EQUAL((constData[{0, 1}]), (data[{0, 1}]));
  BOOST_CHECK_EQUAL((constData[{0, 2}]), (data[{0, 2}]));
  BOOST_CHECK_EQUAL((constData[{1, 0}]), (data[{1, 0}]));
  BOOST_CHECK_EQUAL((constData[{1, 1}]), (data[{1, 1}]));
  BOOST_CHECK_EQUAL((constData[{1, 2}]), (data[{1, 2}]));
  BOOST_CHECK_EQUAL(constData.at({0, 0}), data.at({0, 0}));
  BOOST_CHECK_EQUAL(constData.at({0, 1}), data.at({0, 1}));
  BOOST_CHECK_EQUAL(constData.at({0, 2}), data.at({0, 2}));
  BOOST_CHECK_EQUAL(constData.at({1, 0}), data.at({1, 0}));
  BOOST_CHECK_EQUAL(constData.at({1, 1}), data.at({1, 1}));
  BOOST_CHECK_EQUAL(constData.at({1, 2}), data.at({1, 2}));

  BOOST_CHECK_THROW(constData.at({0, 3}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 4}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 4}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 4}), std::out_of_range);
  
  
  data.fill(14);
  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      BOOST_CHECK_EQUAL((data[{ c, t }]), 14);
  
  data.apply([](int& v){ v *= 2; });
  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      BOOST_CHECK_EQUAL((data[{ c, t }]), 28);
  
  Summer<int> summer;
  static_assert(std::is_same_v<decltype(data.apply(summer)), Summer<int>&>);
  data.apply(summer);
  BOOST_CHECK_EQUAL(summer.get(), N * 28);
  
  summer.reset();
  static_assert
    (std::is_same_v<decltype(constData.apply(summer)), Summer<int>&>);
  constData.apply(summer);
  BOOST_CHECK_EQUAL(summer.get(), N * 28);
  
  auto summer1 = data.apply(Summer<int>{});
  BOOST_CHECK_EQUAL(summer1.get(), N * 28);
  
  auto summer2 = constData.apply(Summer<int>{});
  BOOST_CHECK_EQUAL(summer2.get(), N * 28);
  
  data.clear();
  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      BOOST_CHECK_EQUAL((data[{ c, t }]), 0);
  
  
} // TPCDataContainerTest()


//------------------------------------------------------------------------------
void PlaneDataContainerTest() {

  static constexpr std::size_t NCryostats = 2U;
  static constexpr std::size_t NTPCs      = 3U;
  static constexpr std::size_t NPlanes    = 2U;
  static constexpr std::size_t N          = NCryostats * NTPCs * NPlanes;
  geo::PlaneDataContainer<int> data(NCryostats, NTPCs, NPlanes);

  BOOST_CHECK(!data.empty());
  BOOST_CHECK_EQUAL(data.size(), N);
  BOOST_CHECK_GE(data.capacity(), N);

  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      for (auto p: util::counter<unsigned int>(NPlanes)) 
        BOOST_CHECK_EQUAL((data[{ c, t, p }]), 0);
  
  BOOST_CHECK_EQUAL(data.firstID(), geo::PlaneID(0, 0, 0));
  BOOST_CHECK_EQUAL(data.lastID(), geo::PlaneID(1, 2, 1));

  BOOST_CHECK( data.hasPlane({ 0, 0, 0}));
  BOOST_CHECK( data.hasPlane({ 0, 0, 1}));
  BOOST_CHECK(!data.hasPlane({ 0, 0, 2}));
  BOOST_CHECK( data.hasPlane({ 0, 1, 0}));
  BOOST_CHECK( data.hasPlane({ 0, 1, 1}));
  BOOST_CHECK(!data.hasPlane({ 0, 1, 2}));
  BOOST_CHECK( data.hasPlane({ 0, 2, 0}));
  BOOST_CHECK( data.hasPlane({ 0, 2, 1}));
  BOOST_CHECK(!data.hasPlane({ 0, 2, 2}));
  BOOST_CHECK(!data.hasPlane({ 0, 3, 0}));
  BOOST_CHECK(!data.hasPlane({ 0, 3, 1}));
  BOOST_CHECK(!data.hasPlane({ 0, 3, 2}));
  BOOST_CHECK(!data.hasPlane({ 0, 4, 0}));
  BOOST_CHECK(!data.hasPlane({ 0, 4, 1}));
  BOOST_CHECK(!data.hasPlane({ 0, 4, 2}));
  BOOST_CHECK( data.hasPlane({ 1, 0, 0}));
  BOOST_CHECK( data.hasPlane({ 1, 0, 1}));
  BOOST_CHECK(!data.hasPlane({ 1, 0, 2}));
  BOOST_CHECK( data.hasPlane({ 1, 1, 0}));
  BOOST_CHECK( data.hasPlane({ 1, 1, 1}));
  BOOST_CHECK(!data.hasPlane({ 1, 1, 2}));
  BOOST_CHECK( data.hasPlane({ 1, 2, 0}));
  BOOST_CHECK( data.hasPlane({ 1, 2, 1}));
  BOOST_CHECK(!data.hasPlane({ 1, 2, 2}));
  BOOST_CHECK(!data.hasPlane({ 1, 3, 0}));
  BOOST_CHECK(!data.hasPlane({ 1, 3, 1}));
  BOOST_CHECK(!data.hasPlane({ 1, 3, 2}));
  BOOST_CHECK(!data.hasPlane({ 1, 4, 0}));
  BOOST_CHECK(!data.hasPlane({ 1, 4, 1}));
  BOOST_CHECK(!data.hasPlane({ 1, 4, 2}));
  BOOST_CHECK(!data.hasPlane({ 2, 0, 0}));
  BOOST_CHECK(!data.hasPlane({ 2, 0, 1}));
  BOOST_CHECK(!data.hasPlane({ 2, 0, 2}));
  BOOST_CHECK(!data.hasPlane({ 2, 1, 0}));
  BOOST_CHECK(!data.hasPlane({ 2, 1, 1}));
  BOOST_CHECK(!data.hasPlane({ 2, 1, 2}));
  BOOST_CHECK(!data.hasPlane({ 2, 2, 0}));
  BOOST_CHECK(!data.hasPlane({ 2, 2, 1}));
  BOOST_CHECK(!data.hasPlane({ 2, 2, 2}));
  BOOST_CHECK(!data.hasPlane({ 2, 3, 0}));
  BOOST_CHECK(!data.hasPlane({ 2, 3, 1}));
  BOOST_CHECK(!data.hasPlane({ 2, 3, 2}));
  BOOST_CHECK(!data.hasPlane({ 2, 4, 0}));
  BOOST_CHECK(!data.hasPlane({ 2, 4, 1}));
  BOOST_CHECK(!data.hasPlane({ 2, 4, 2}));

  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 0, 0}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 0, 1}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 0, 2}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 1, 0}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 1, 1}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 1, 2}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 2, 0}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 2, 1}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 0, 2, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 0, 3, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 0, 3, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 0, 3, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 0, 4, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 0, 4, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 0, 4, 2}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 0, 0}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 0, 1}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 0, 2}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 1, 0}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 1, 1}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 1, 2}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 2, 0}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 2, 1}));
  BOOST_CHECK( data.hasTPC(geo::PlaneID{ 1, 2, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 1, 3, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 1, 3, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 1, 3, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 1, 4, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 1, 4, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 1, 4, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 0, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 0, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 0, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 1, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 1, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 1, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 2, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 2, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 2, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 3, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 3, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 3, 2}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 4, 0}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 4, 1}));
  BOOST_CHECK(!data.hasTPC(geo::PlaneID{ 2, 4, 2}));

  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 0, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 0, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 0, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 1, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 1, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 1, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 2, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 2, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 2, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 3, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 3, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 3, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 4, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 4, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 0, 4, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 0, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 0, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 0, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 1, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 1, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 1, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 2, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 2, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 2, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 3, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 3, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 3, 2}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 4, 0}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 4, 1}));
  BOOST_CHECK( data.hasCryostat(geo::PlaneID{ 1, 4, 2}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 0, 0}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 0, 1}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 0, 2}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 1, 0}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 1, 1}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 1, 2}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 2, 0}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 2, 1}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 2, 2}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 3, 0}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 3, 1}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 3, 2}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 4, 0}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 4, 1}));
  BOOST_CHECK(!data.hasCryostat(geo::PlaneID{ 2, 4, 2}));


  data[{0, 0, 0}] = 4;
  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   4);
  BOOST_CHECK_EQUAL(data.at({0, 0, 0}),    4);
  data[{0, 0, 0}] = 5;
  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(data.at({0, 0, 0}),    5);

  data[{0, 0, 1}] = 6;
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(data.at({0, 0, 1}),    6);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);

  data[{0, 1, 0}] = 15;
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(data.at({0, 1, 0}),   15);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);

  data[{0, 1, 1}] = 16;
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(data.at({0, 1, 1}),   16);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);

  data[{0, 2, 0}] = 25;
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);
  BOOST_CHECK_EQUAL(data.at({0, 2, 0}),   25);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);

  data[{0, 2, 1}] = 26;
  BOOST_CHECK_EQUAL(  (data[{0, 2, 1}]),  26);
  BOOST_CHECK_EQUAL(data.at({0, 2, 1}),   26);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);

  data[{1, 0, 0}] = 105;
  BOOST_CHECK_EQUAL(  (data[{1, 0, 0}]), 105);
  BOOST_CHECK_EQUAL(data.at({1, 0, 0}),  105);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 1}]),  26);

  data[{1, 0, 1}] = 106;
  BOOST_CHECK_EQUAL(  (data[{1, 0, 1}]), 106);
  BOOST_CHECK_EQUAL(data.at({1, 0, 1}),  106);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 1}]),  26);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 0}]), 105);

  data[{1, 1, 0}] = 115;
  BOOST_CHECK_EQUAL(  (data[{1, 1, 0}]), 115);
  BOOST_CHECK_EQUAL(data.at({1, 1, 0}),  115);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 1}]),  26);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 0}]), 105);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 1}]), 106);

  data[{1, 1, 1}] = 116;
  BOOST_CHECK_EQUAL(  (data[{1, 1, 1}]), 116);
  BOOST_CHECK_EQUAL(data.at({1, 1, 1}),  116);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 1}]),  26);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 0}]), 105);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 1}]), 106);
  BOOST_CHECK_EQUAL(  (data[{1, 1, 0}]), 115);

  data[{1, 2, 0}] = 125;
  BOOST_CHECK_EQUAL(  (data[{1, 2, 0}]), 125);
  BOOST_CHECK_EQUAL(data.at({1, 2, 0}),  125);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 1}]),  26);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 0}]), 105);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 1}]), 106);
  BOOST_CHECK_EQUAL(  (data[{1, 1, 0}]), 115);
  BOOST_CHECK_EQUAL(  (data[{1, 1, 1}]), 116);

  data[{1, 2, 1}] = 126;
  BOOST_CHECK_EQUAL(  (data[{1, 2, 1}]), 126);
  BOOST_CHECK_EQUAL(data.at({1, 2, 1}),  126);

  BOOST_CHECK_EQUAL(  (data[{0, 0, 0}]),   5);
  BOOST_CHECK_EQUAL(  (data[{0, 0, 1}]),   6);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 0}]),  15);
  BOOST_CHECK_EQUAL(  (data[{0, 1, 1}]),  16);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 0}]),  25);
  BOOST_CHECK_EQUAL(  (data[{0, 2, 1}]),  26);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 0}]), 105);
  BOOST_CHECK_EQUAL(  (data[{1, 0, 1}]), 106);
  BOOST_CHECK_EQUAL(  (data[{1, 1, 0}]), 115);
  BOOST_CHECK_EQUAL(  (data[{1, 1, 1}]), 116);
  BOOST_CHECK_EQUAL(  (data[{1, 2, 0}]), 125);


  BOOST_CHECK_THROW(data.at({0, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({0, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({1, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(data.at({2, 3, 2}), std::out_of_range);

  BOOST_CHECK_EQUAL(data.first(), 5);
  data.first() = -5;
  BOOST_CHECK_EQUAL((data[{0, 0, 0}]), -5);
  BOOST_CHECK_EQUAL(data.first(), -5);
  data.first() =  5;

  BOOST_CHECK_EQUAL(data.last(), 126);
  data.last() = -126;
  BOOST_CHECK_EQUAL((data[{1U, 2U, 1U}]), -126);
  BOOST_CHECK_EQUAL(data.last(), -126);
  data.last() =  126;

  auto const& constData = data;

  BOOST_CHECK_EQUAL
    (std::addressof(constData.first()), std::addressof(data.first()));
  BOOST_CHECK_EQUAL
    (std::addressof(constData.last()), std::addressof(data.last()));

  BOOST_CHECK_EQUAL((constData[{0, 0, 0}]), (data[{0, 0, 0}]));
  BOOST_CHECK_EQUAL((constData[{0, 0, 1}]), (data[{0, 0, 1}]));
  BOOST_CHECK_EQUAL((constData[{0, 1, 0}]), (data[{0, 1, 0}]));
  BOOST_CHECK_EQUAL((constData[{0, 1, 1}]), (data[{0, 1, 1}]));
  BOOST_CHECK_EQUAL((constData[{0, 2, 0}]), (data[{0, 2, 0}]));
  BOOST_CHECK_EQUAL((constData[{0, 2, 1}]), (data[{0, 2, 1}]));
  BOOST_CHECK_EQUAL((constData[{1, 0, 0}]), (data[{1, 0, 0}]));
  BOOST_CHECK_EQUAL((constData[{1, 0, 1}]), (data[{1, 0, 1}]));
  BOOST_CHECK_EQUAL((constData[{1, 1, 0}]), (data[{1, 1, 0}]));
  BOOST_CHECK_EQUAL((constData[{1, 1, 1}]), (data[{1, 1, 1}]));
  BOOST_CHECK_EQUAL((constData[{1, 2, 0}]), (data[{1, 2, 0}]));
  BOOST_CHECK_EQUAL((constData[{1, 2, 1}]), (data[{1, 2, 1}]));
  BOOST_CHECK_EQUAL(constData.at({0, 0, 0}), data.at({0, 0, 0}));
  BOOST_CHECK_EQUAL(constData.at({0, 0, 1}), data.at({0, 0, 1}));
  BOOST_CHECK_EQUAL(constData.at({0, 1, 0}), data.at({0, 1, 0}));
  BOOST_CHECK_EQUAL(constData.at({0, 1, 1}), data.at({0, 1, 1}));
  BOOST_CHECK_EQUAL(constData.at({0, 2, 0}), data.at({0, 2, 0}));
  BOOST_CHECK_EQUAL(constData.at({0, 2, 1}), data.at({0, 2, 1}));
  BOOST_CHECK_EQUAL(constData.at({1, 0, 0}), data.at({1, 0, 0}));
  BOOST_CHECK_EQUAL(constData.at({1, 0, 1}), data.at({1, 0, 1}));
  BOOST_CHECK_EQUAL(constData.at({1, 1, 0}), data.at({1, 1, 0}));
  BOOST_CHECK_EQUAL(constData.at({1, 1, 1}), data.at({1, 1, 1}));
  BOOST_CHECK_EQUAL(constData.at({1, 2, 0}), data.at({1, 2, 0}));
  BOOST_CHECK_EQUAL(constData.at({1, 2, 1}), data.at({1, 2, 1}));

  BOOST_CHECK_THROW(constData.at({0, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 4, 0}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 4, 1}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({0, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({1, 3, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 0, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 1, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 2, 2}), std::out_of_range);
  BOOST_CHECK_THROW(constData.at({2, 3, 2}), std::out_of_range);


  data.fill(14);
  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      for (auto p: util::counter<unsigned int>(NPlanes)) 
        BOOST_CHECK_EQUAL((data[{ c, t, p }]), 14);
  
  data.apply([](int& v){ v *= 2; });
  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      for (auto p: util::counter<unsigned int>(NPlanes)) 
        BOOST_CHECK_EQUAL((data[{ c, t, p }]), 28);
  
  Summer<int> summer;
  static_assert(std::is_same_v<decltype(data.apply(summer)), Summer<int>&>);
  data.apply(summer);
  BOOST_CHECK_EQUAL(summer.get(), N * 28);
  
  summer.reset();
  static_assert
    (std::is_same_v<decltype(constData.apply(summer)), Summer<int>&>);
  constData.apply(summer);
  BOOST_CHECK_EQUAL(summer.get(), N * 28);
  
  auto summer1 = data.apply(Summer<int>{});
  BOOST_CHECK_EQUAL(summer1.get(), N * 28);
  
  auto summer2 = constData.apply(Summer<int>{});
  BOOST_CHECK_EQUAL(summer2.get(), N * 28);
  
  data.clear();
  for (auto c: util::counter<unsigned int>(NCryostats)) 
    for (auto t: util::counter<unsigned int>(NTPCs)) 
      for (auto p: util::counter<unsigned int>(NPlanes)) 
        BOOST_CHECK_EQUAL((data[{ c, t, p }]), 0);
  
  
} // PlaneDataContainerTest()


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(TPCDataContainerTestCase) {
  TPCDataContainerTest();
} // TPCDataContainerTestCase


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(PlaneDataContainerTestCase) {
  PlaneDataContainerTest();
} // PlaneDataContainerTestCase


//------------------------------------------------------------------------------
