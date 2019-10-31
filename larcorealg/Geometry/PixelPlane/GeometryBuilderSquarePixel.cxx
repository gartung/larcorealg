/**
 * @file   larcorealg/Geometry/PixelPlane/GeometryBuilderSquarePixel.cxx
 * @brief  Implementation of pixel-based geometry extractor (impl. file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 23, 2019
 * @see    `larcorealg/Geometry/PixelPlane/GeometryBuilderSquarePixel.h`
 */

// library header
#include "larcorealg/Geometry/PixelPlane/GeometryBuilderSquarePixel.h"

// LArSoft libraries
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"


//------------------------------------------------------------------------------
geo::GeometryBuilderSquarePixel::GeometryBuilderSquarePixel
  (Config const& config)
  : geo::GeometryBuilderStandard(config)
{
  MF_LOG_TRACE("GeometryBuilder")
    << "Loading geometry builder: GeometryBuilderSquarePixel";
} // geo::GeometryBuilderSquarePixel::GeometryBuilderSquarePixel()


//------------------------------------------------------------------------------
auto geo::GeometryBuilderSquarePixel::doMakePlane(Path_t& path) -> PlanePtr_t {
  return std::make_unique<geo::PixelPlaneGeo>(
    path.current(), path.currentTransformation<geo::TransformationMatrix>()
    );
} // geo::GeometryBuilderSquarePixel::doMakePlane()


//------------------------------------------------------------------------------

