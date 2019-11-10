/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h
 * @brief   Readout plane for a pixel detector.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    November 5, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.cxx`
 * @ingroup GeometryPixel
 */

#ifndef LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H
#define LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeoBase.h"

// C/C++ standard libraries
#include <regex>


// -----------------------------------------------------------------------------
namespace geo { class PixelPlaneGeo; }

/**
 * @brief Geometry description for a pixel readout plane.
 * 
 * This class implements the `geo:PlaneGeo` interface via
 * `geo::PixelPlaneGeoBase` base class, which should be referred to for any
 * detail on the class behaviour.
 * 
 * This class is an implementation of `geo::PixelPlaneGeoBase` which only adds
 * a specific initialization pattern to it.
 * 
 * 
 * Initialization
 * ===============
 * 
 * This class uses the standard pixel geometry initialization algorithm
 * `geo::PixelPlaneGeoBase::InitializePixelGeometry()`.
 * The initialization of the pixel geometry happens immediately at construction,
 * just after the base class initialization.
 * 
 * The pixel geometry information needed by the algorithm is obtained from the
 * GDML/ROOT description of the plane geometry.
 * First, direct information is sought as metadata in the geometry description
 * (see `geo::PixelPlaneGeoBase::ReadPixelGeometryFromMetadata()`), then the
 * rest of the information is filled from pixel volume information in the same
 * geometry (see `geo::PixelPlaneGeoBase::ExtractPixelGeometry()`).
 * 
 */
class geo::PixelPlaneGeo: public geo::PixelPlaneGeoBase {
  
    public:
  
  /// Default pattern for recognizing by name a sensitive element geometry node.
  static std::regex const DefaultPixelPattern;
  
  /**
   * @brief Constructor: extracts pixel information from geometry description.
   * @param node GDML/ROOT object describing the plane itself
   * @param trans transformation describing position and orientation of the
   *              plane in the world (based on local-to-world transformation)
   * @param pixelNamePattern name of pixel nodes in the geometry description
   * 
   * The pattern is a regular expression: every volume matching it will be
   * considered as a sensitive element (i.e. a pixel).
   */
  PixelPlaneGeo(
    TGeoNode const& node, geo::TransformationMatrix&& trans,
    std::regex const& pixelNamePattern = DefaultPixelPattern
    );
  
  
}; // class geo::PixelPlaneGeo


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H
