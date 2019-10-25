/**
 * @file   larcorealg/Geometry/GeometryData.h
 * @brief  Simple data structure holding the data of the geometry.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @see    `larcorealg/Geometry/GeometryCore.h`
 * @ingroup Geometry
 * 
 * This is a header-only library.
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYDATA_H
#define LARCOREALG_GEOMETRY_GEOMETRYDATA_H

// LArSoft libraries
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/AuxDetGeo.h"

// C/C++ standard libraries
#include <vector>


namespace geo {
  
  /*
   * To decouple the channel mapping algorithm from `geo::GeometryCore`,
   * this data structure, which is in fact the core data of `geo::GeometryCore`,
   * is detached and included by both.
   */
  
  /// Data in the geometry description.
  struct GeometryData_t {

    /// Type of list of cryostats.
    using CryostatList_t = std::vector<geo::CryostatGeo>;
    /// Type of list of auxiliary detectors
    using AuxDetList_t = std::vector<geo::AuxDetGeo*>;

    CryostatList_t cryostats; ///< The detector cryostats.
    AuxDetList_t   auxDets;   ///< The auxiliary detectors.

  }; // GeometryData_t

} // namespace geo

#endif // LARCOREALG_GEOMETRY_GEOMETRYDATA_H
