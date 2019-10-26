/**
 * @file   larcorealg/Geometry/GeoObjectSorterNull.h
 * @brief  A geometry sorter that does not sort anything.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 25, 2019
 * 
 * There is no implementation file for this class.
 */

#ifndef LARCOREALG_GEOMETRY_GEOOBJECTSORTERNULL_H
#define LARCOREALG_GEOMETRY_GEOOBJECTSORTERNULL_H

// LArSoft libraries
#include "larcorealg/Geometry/GeoObjectSorter.h"


// -----------------------------------------------------------------------------
namespace geo { class GeoObjectSorterNull; }

/**
 * @brief A geometry sorter that does not sort anything.
 * @see `geo::GeoObjectSorter`
 * @ingroup Geometry
 */
class geo::GeoObjectSorterNull: public geo::GeoObjectSorter {

    public:

  virtual void SortAuxDets
    (std::vector<geo::AuxDetGeo*>&) const override {}
  
  virtual void SortAuxDetSensitive
    (std::vector<geo::AuxDetSensitiveGeo*>&) const override {}
  
  virtual void SortCryostats
    (std::vector<geo::CryostatGeo*>&) const override {}

  virtual void SortTPCs
    (std::vector<geo::TPCGeo*>&) const override {}

  virtual void SortPlanes
    (std::vector<geo::PlaneGeo*>&, geo::DriftDirection_t const&) const override
    {}

  virtual void SortWires
    (std::vector<geo::WireGeo*>&) const override {}

  virtual void SortOpDets
    (std::vector<geo::OpDetGeo*>&) const override {}

}; // class geo::GeoObjectSorterNull


// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEOOBJECTSORTERNULL_H
