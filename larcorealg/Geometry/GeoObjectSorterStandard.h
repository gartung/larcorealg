////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterStandard.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERSTANDARD_H
#define GEO_GEOOBJECTSORTERSTANDARD_H

#include <vector>

#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"  // for DriftDirec...

namespace fhicl { class ParameterSet; }

namespace geo{

  class AuxDetGeo;
  class AuxDetSensitiveGeo;
  class CryostatGeo;
  class PlaneGeo;
  class TPCGeo;
  class WireGeo;

  /// @ingroup Geometry
  class GeoObjectSorterStandard : public GeoObjectSorter {
  public:

    GeoObjectSorterStandard(fhicl::ParameterSet const& p);
    ~GeoObjectSorterStandard();

    void SortAuxDets        (std::vector<geo::AuxDetGeo*>          & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo)   const;
    void SortCryostats      (std::vector<geo::CryostatGeo*>        & cgeo)     const;
    void SortTPCs     	    (std::vector<geo::TPCGeo*>      	   & tgeo)     const;
    void SortPlanes   	    (std::vector<geo::PlaneGeo*>    	   & pgeo,
		      	     geo::DriftDirection_t     	     const & driftDir) const;
    void SortWires    	    (std::vector<geo::WireGeo*>     	   & wgeo)     const;

  private:

  };

}

#endif // GEO_GEOOBJECTSORTERSTANDARD_H
