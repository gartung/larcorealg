////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTER_H
#define GEO_GEOOBJECTSORTER_H

#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace geo{

  class AuxDetGeo;
  class AuxDetSensitiveGeo;
  class CryostatGeo;
  class TPCGeo;
  class PlaneGeo;
  class WireGeo;
  class OpDetGeo;
  /// @ingroup Geometry
  class GeoObjectSorter {

  public:

    virtual ~GeoObjectSorter() = default;

    virtual void SortAuxDets        (std::vector<geo::AuxDetGeo*>          & adgeo)   const = 0;
    virtual void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo)  const = 0;
    virtual void SortCryostats      (std::vector<geo::CryostatGeo*>        & cgeo)    const = 0;
    virtual void SortTPCs     	    (std::vector<geo::TPCGeo*>      	  & tgeo)     const = 0;
    virtual void SortPlanes   	    (std::vector<geo::PlaneGeo*>       	  & pgeo,
			      	     geo::DriftDirection_t     	    const & driftDir) const = 0;
    virtual void SortWires    	    (std::vector<geo::WireGeo*>     	  & wgeo)     const = 0;
    virtual void SortOpDets        (std::vector<geo::OpDetGeo*>          & opdet) const;
  private:

  };

}

#endif // GEO_GEOOBJECTSORTER_H
