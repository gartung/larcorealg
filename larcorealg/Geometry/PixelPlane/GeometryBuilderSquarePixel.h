/**
 * @file   larcorealg/Geometry/PixelPlane/GeometryBuilderSquarePixel.h
 * @brief  Implementation of pixel-based detector geometry extractor.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 23, 2019
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/PixelPlane/GeometryBuilderSquarePixel.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_PIXELPLANE_GEOMETRYBUILDERSQUAREPIXEL_H
#define LARCOREALG_GEOMETRY_PIXELPLANE_GEOMETRYBUILDERSQUAREPIXEL_H

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"



namespace geo {
  
  /**
   * @brief Geometry builder which uses `geo::PixelPlane` for sensitive planes.
   * 
   * This builder works like `geo::GeometryBuilderWireless`, with the exception
   * that it creates `geo::PixelPlane` instead of `geo::WirePlaneGeo` objects.
   * The content of the plane volumes in the geometry is ignored.
   */
  class GeometryBuilderSquarePixel: public geo::GeometryBuilderStandard {
    
      public:
    
    struct Config: public geo::GeometryBuilderStandard::Config {
      
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      // so far, only the standard configuration is needed
      
    }; // struct Config
    
    
    /// Constructor: uses the specified configuration.
    GeometryBuilderSquarePixel(Config const& config);
    
    //
    // we don't expand the public interface here
    //
    
      protected:
    
    // --- BEGIN Plane information ---------------------------------------------
    /// @name Pixel plane information
    /// @{

    using PlanePtr_t = std::unique_ptr<geo::PlaneGeo>;
    using Planes_t = std::vector<PlanePtr_t>;

    /// Create a `geo::PixelPlane` for each plane volume in the geometry.
    virtual PlanePtr_t doMakePlane(Path_t& path) override;

    /// @}
    // --- END Plane information -----------------------------------------------


    // --- BEGIN Wire information ----------------------------------------------
    /// @name Wire information
    /// @{
    
    /// Core implementation of `extractWires()`: no wires returned whatsoever.
    virtual Wires_t doExtractWires(Path_t&) { return {}; }
    
    /// @}
    // --- END Wire information ------------------------------------------------
    
    
  }; // class GeometryBuilderSquarePixel
  
} // namespace geo


#endif // LARCOREALG_GEOMETRY_PIXELPLANE_GEOMETRYBUILDERSQUAREPIXEL_H
