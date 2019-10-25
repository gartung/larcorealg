/**
 * @file   larcorealg/Geometry/details/extractMaxGeometryElements.h
 * @brief  Algorithm discovering the number of elements in the geometry.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 19, 2019
 * 
 * This is a header only library.
 */

#ifndef LARCOREALG_GEOMETRY_DETAILS_EXTRACTMAXGEOMETRYELEMENTS_H
#define LARCOREALG_GEOMETRY_DETAILS_EXTRACTMAXGEOMETRYELEMENTS_H


// LArSoft libraries
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryData.h"

// C/C++ standard libraries
#include <array>
#include <cstddef> // std::size_t


namespace geo::details {
  
  /**
   * @brief Extracts the maximum number of elements per type.
   * @tparam Levels the number of detector elements to discover
   * @param Cryostats the sorted list of cryostats in the detector
   * @return an array with maximum number of cryostats, TPCs, planes and wires
   * 
   * The returned array includes:
   * 
   * * index `[0]`: number of cryostats
   * * index `[1]`: maximum number of TPCs in any of the cryostats
   *     (enabled only if `Levels` is `2` or higher)
   * * index `[2]`: maximum number of wire planes in any of the TPCs
   *     (enabled only if `Levels` is `3` or higher)
   * * index `[3]`: maximum number of wires in any of the wire planes
   *     (enabled only if `Levels` is `4`)
   * 
   */
  template <std::size_t Levels = 4U>
  static std::array<unsigned int, Levels> extractMaxGeometryElements
    (geo::GeometryData_t::CryostatList_t const& Cryostats);

} // namespace geo::details


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
template <std::size_t Levels /* = 4U */>
std::array<unsigned int, Levels>
geo::details::extractMaxGeometryElements
  (geo::GeometryData_t::CryostatList_t const& Cryostats)
{
  static_assert(Levels > 0U);
  static_assert(Levels <= 4U);
  
  std::array<unsigned int, Levels> maxElements;
  maxElements.fill(0U);
  
  auto setMax = [&maxElements](std::size_t index, unsigned int value)
    { if (maxElements[index] < value) maxElements[index] = value; };
  
  setMax(0U, Cryostats.size());
  if constexpr(Levels > 1U) {
    for (geo::CryostatGeo const& cryo: Cryostats) {
      setMax(1U, cryo.NTPC());
      if constexpr(Levels > 2U) {
        for (geo::TPCGeo const& TPC: cryo.IterateTPCs()) {
          setMax(2U, TPC.Nplanes());
          if constexpr(Levels > 3U) {
            for (geo::PlaneGeo const& plane: TPC.IteratePlanes()) {
              setMax(3U, plane.Nwires());
            } // for planes
          } // if do wires
        } // for TPCs
      } // if do planes
    } // for cryostats
  } // if do TPCs
  
  return maxElements;
} // geo::details::extractMaxGeometryElements()


// -----------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_DETAILS_EXTRACTMAXGEOMETRYELEMENTS_H
