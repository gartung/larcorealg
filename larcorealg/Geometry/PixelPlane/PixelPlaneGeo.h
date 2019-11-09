/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h
 * @brief   Representation of a readout plane with pixels as sensitive elements.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 23, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.cxx`
 * @ingroup Geometry
 */

#ifndef LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H
#define LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeoBase.h"

// C/C++ standard libraries
#include <regex>


// -----------------------------------------------------------------------------
namespace geo {
  class PixelPlaneGeo;
} // namespace geo


/**
 * @brief Geometry information for a single pixel plane.
 * @ingroup Geometry
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
 * The details are documented in `ExtractPixelGeometry()` method.
 * 
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
  
  
    private:
  
  // --- BEGIN --- Initialization procedures -----------------------------------
  
  /**
   * @brief Extracts in some way the pixel information.
   * @param pixelNamePattern pattern used to recognize pixel nodes
   * @return the pixel geometry information as complete as possible
   * 
   * This method produces information about the placement and geometry of the
   * pixels on the plane, to be taken by
   * `geo::PixelPlaneGeoBase::InitializePixelGeometry()` and turned into actual
   * plane set up.
   * 
   * @todo Document the extraction algorithm, expectations and abilities.
   */
  RectPixelGeometry_t ExtractPixelGeometry
    (std::regex const& pixelNamePattern) const;
  
  // --- END --- Initialization procedures -------------------------------------
  
  
}; // class PixelPlaneGeo


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H
