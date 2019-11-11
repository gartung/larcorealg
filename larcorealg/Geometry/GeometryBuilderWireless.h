/**
 * @file   larcorealg/Geometry/GeometryBuilderWireless.h
 * @brief  Implementation of wireless geometry extractor.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 26, 2019
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/GeometryBuilderStandard.cxx`
 *         `larcorealg/Geometry/GeometryBuilderWireless.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYBUILDERWIRELESS_H
#define LARCOREALG_GEOMETRY_GEOMETRYBUILDERWIRELESS_H

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"



namespace geo {
  
  /**
   * @brief Geometry builder which ignores wires on wire planes.
   * 
   * This builder works like `geo::GeometryBuilderStandard`, with the exception
   * that it does not consider the wires on the wire plane objects: wires may
   * or may not exist.
   * 
   */
  class GeometryBuilderWireless: public geo::GeometryBuilderStandard {
    
      public:
    
    /// Constructor: uses the specified configuration.
    GeometryBuilderWireless(Config const& config);
    
    //
    // we don't expand the public interface here
    //
    
      protected:
    
    // --- BEGIN Wire information ----------------------------------------------
    /// @name Wire information
    /// @{
    
    /// Core implementation of `extractWires()`: no wires returned whatsoever.
    virtual Wires_t doExtractWires(Path_t&) override { return {}; }
    
    /// @}
    // --- END Wire information ------------------------------------------------
    
    
  }; // class GeometryBuilderWireless
  
} // namespace geo


#endif // LARCOREALG_GEOMETRY_GEOMETRYBUILDERWIRELESS_H
