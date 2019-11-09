/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeo.cxx
 * @brief   Representation of a readout plane with pixels as sensitive elements.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 23, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`
 * @ingroup Geometry
 */

// class header
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Zaxis()
#if 0
#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "larcorealg/Geometry/SimpleGeo.h" // lar::util::simple_geo::Rectangle
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::rounded01()
#include "larcorealg/CoreUtils/RealComparisons.h" // makeVector3DComparison()
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#endif // 0

// C/C++ standard library
#include <utility> // std::move()
#if 0
#include <sstream> // std::ostringstream
#include <algorithm> // std::clamp(), std::swap()
#include <array>
#include <limits> // std::numeric_limits<>
#include <cmath> // std::cos(), std::abs(), ...
#include <cassert>

#endif // 0


// -----------------------------------------------------------------------------
// --- geo::PixelPlaneGeo
// -----------------------------------------------------------------------------
std::regex const geo::PixelPlaneGeo::DefaultPixelPattern {
  ".*pixel.*",
  std::regex::basic | std::regex::icase | std::regex::optimize
  };


// -----------------------------------------------------------------------------
geo::PixelPlaneGeo::PixelPlaneGeo(
  TGeoNode const& node, geo::TransformationMatrix&& trans,
  std::regex const& pixelNamePattern /* = DefaultPixelPattern */
  )
  : geo::PixelPlaneGeoBase(node, std::move(trans))
{

  /*
   * REMINDER: do not rely on virtual methods from derived classes here, as
   *           they might not be available yet (the derived class constructor
   *           hasn't been run yet at this point)
   */
  
  RectPixelGeometry_t const pixelGeometry
    = ExtractPixelGeometry(pixelNamePattern);
  
  // use the base class standard initialization for pixel geometry:
  geo::PixelPlaneGeoBase::InitializePixelGeometry(pixelGeometry);
  
} // geo::PixelPlaneGeo::PixelPlaneGeo()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::ExtractPixelGeometry
  (std::regex const& pixelNamePattern) const
  -> RectPixelGeometry_t
{
  
  return {
    { // axis specifications:
      RectPixelGeometry_t::AxisInfo_t {
        geo::Yaxis<LocalVector_t>(), // pixels along local y axis
        Depth(),                     // full `Depth()` length
        std::nullopt,                // number of pixels not specified
        3.0                          // 3 cm pitch along local y axis
      },
      RectPixelGeometry_t::AxisInfo_t {
        geo::Zaxis<LocalVector_t>(), // pixels along local z axis
        Width(),                     // full `Width()` length
        std::nullopt,                // number of pixels not specified
        4.0                          // 4 cm pitch along local y axis
      }
    },
    std::nullopt // origin of the pixels not specified
  };
  
} // geo::PixelPlaneGeo::ExtractPixelGeometry()


// -----------------------------------------------------------------------------
