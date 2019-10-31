/**
 * @file   larcorealg/Geometry/PixelPlane/ChannelMapPixelAlg.cxx
 * @brief  Simple channel mapping for pixel readout (imlementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 24, 2019
 * @see    `larcorealg/Geometry/PixelPlane/ChannelMapPixelAlg.h`
 */

// library header
#include "larcorealg/Geometry/PixelPlane/ChannelMapPixelAlg.h"

// LArSoft libraries
#include "larcorealg/Geometry/details/extractMaxGeometryElements.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug::printBacktrace()
#include "larcorealg/CoreUtils/counter.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"


// -----------------------------------------------------------------------------
geo::ChannelMapPixelAlg::ChannelMapPixelAlg(Config const& config) {
  
  MF_LOG_DEBUG("ChannelMapAlg")
    << "Loading channel mapping algorithm: ChannelMapPixelAlg";
  
} // geo::ChannelMapPixelAlg::ChannelMapPixelAlg()


// -----------------------------------------------------------------------------
void geo::ChannelMapPixelAlg::Initialize(geo::GeometryData_t const& geodata) {
  
  fillFirstChannelInfo(geodata.cryostats);
  
  fillGeometrySizes(geodata.cryostats);
  
  makePlaneChannelMap();
  
} // geo::ChannelMapPixelAlg::Initialize()


// -----------------------------------------------------------------------------
void geo::ChannelMapPixelAlg::Uninitialize() {
  
  fFirstChannelInfo.clear();
  fGeometrySizes.reset();
  fPlaneInfo.reset();
  
} // geo::ChannelMapPixelAlg::Uninitialize()


// -----------------------------------------------------------------------------
std::vector<geo::WireID> geo::ChannelMapPixelAlg::ChannelToWire
  (raw::ChannelID_t channel) const
{
  assert(!fFirstChannelInfo.empty());
  auto const* channelInfo = fFirstChannelInfo.find(channel);
  if (!channelInfo) {
    throw cet::exception("Geometry")
      << "geo::ChannelMapPixelAlg::ChannelToWire(" << channel
      << "): invalid channel requested (must be lower than "
      << Nchannels() << ")\n";
  }
  return {
    geo::WireID{ channelInfo->data().plane, (channel - channelInfo->channel()) }
    };
} // geo::ChannelMapPixelAlg::ChannelToWire()


// -----------------------------------------------------------------------------
unsigned int geo::ChannelMapPixelAlg::Nchannels() const {
  assert(!fFirstChannelInfo.empty());
  return static_cast<unsigned int>(fFirstChannelInfo.endChannel());
} // geo::ChannelMapPixelAlg::Nchannels()


// -----------------------------------------------------------------------------
unsigned int geo::ChannelMapPixelAlg::Nchannels
  (readout::ROPID const& ropid) const
  { return NChannels(ropid); }


// -----------------------------------------------------------------------------
double geo::ChannelMapPixelAlg::WireCoordinate
  (double YPos, double ZPos, geo::PlaneID const& planeID) const
{
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ChannelMapPixelAlg");
  e << "ChannelMapPixelAlg does not support `WireCoordinate()` call."
    "\nPlease update calling software to use geo::PlaneGeo::WireCoordinate()`:"
    "\n\n";
  
  lar::debug::printBacktrace(e, 2U);
  
  throw e;
} // geo::ChannelMapPixelAlg::WireCoordinate()


// -----------------------------------------------------------------------------
geo::WireID geo::ChannelMapPixelAlg::NearestWireID
  (const TVector3& worldPos, geo::PlaneID const& planeID) const
{
  
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ChannelMapPixelAlg");
  e << "ChannelMapPixelAlg does not support `NearestWireID()` call."
    "\nPlease update calling software to use geo::PlaneGeo::NearestWireID()`:"
    "\n\n";
  
  lar::debug::printBacktrace(e, 2U);
  
  throw e;
} // geo::ChannelMapPixelAlg::NearestWireID()


// -----------------------------------------------------------------------------
raw::ChannelID_t geo::ChannelMapPixelAlg::PlaneWireToChannel
  (geo::WireID const& wireID) const
{
  return planeInfo(wireID).firstChannel + wireID.Wire;
} // geo::ChannelMapPixelAlg::PlaneWireToChannel()


//------------------------------------------------------------------------------
std::set<geo::PlaneID> const& geo::ChannelMapPixelAlg::PlaneIDs() const {
  
  /*
   * this should NOT be called... it shouldn't be here at all!
   */
  
  cet::exception e("ChannelMapPixelAlg");
  e << "geo::ChannelMapPixelAlg does not support `PlaneIDs()` call."
    "\nPlease update calling software to use geo::GeometryCore::IteratePlanes()`"
    "\n\n";
  
  lar::debug::printBacktrace(e, 2U);
  
  throw e;
  
} // geo::ChannelMapPixelAlg::PlaneIDs()


// -----------------------------------------------------------------------------
unsigned int geo::ChannelMapPixelAlg::NTPCsets
  (readout::CryostatID const& cryoid) const
{
  return hasCryostat(cryoid)? NTPCs(cryoid): 0U;
} // geo::ChannelMapPixelAlg::NTPCsets()


// -----------------------------------------------------------------------------
unsigned int geo::ChannelMapPixelAlg::MaxTPCsets() const {
  return MaxTPCs();
} // geo::ChannelMapPixelAlg::MaxTPCsets()


// -----------------------------------------------------------------------------
bool geo::ChannelMapPixelAlg::HasTPCset
  (readout::TPCsetID const& tpcsetid) const
{
  return tpcsetid && hasTPCset(tpcsetid);
} // geo::ChannelMapPixelAlg::HasTPCset()


// -----------------------------------------------------------------------------
readout::TPCsetID geo::ChannelMapPixelAlg::TPCtoTPCset
  (geo::TPCID const& tpcid) const
{
  return convertTPCtoTPCset(tpcid);
} // geo::ChannelMapPixelAlg::TPCtoTPCset()


// -----------------------------------------------------------------------------
std::vector<geo::TPCID> geo::ChannelMapPixelAlg::TPCsetToTPCs
  (readout::TPCsetID const& tpcsetid) const
{
  std::vector<geo::TPCID> TPCs;
  if (tpcsetid) TPCs.push_back(convertTPCsetToTPC(tpcsetid));
  return TPCs;
} // geo::ChannelMapPixelAlg::TPCsetToTPCs()


// -----------------------------------------------------------------------------
geo::TPCID geo::ChannelMapPixelAlg::FirstTPCinTPCset
  (readout::TPCsetID const& tpcsetid) const
{
  return tpcsetid? convertTPCsetToTPC(tpcsetid): geo::TPCID{};
} // geo::ChannelMapPixelAlg::FirstTPCinTPCset()


// -----------------------------------------------------------------------------
unsigned int geo::ChannelMapPixelAlg::NROPs
  (readout::TPCsetID const& tpcsetid) const
{
  return hasTPCset(tpcsetid)? NPlanes(convertTPCsetToTPC(tpcsetid)): 0U;
} // geo::ChannelMapPixelAlg::NROPs()


// -----------------------------------------------------------------------------
unsigned int geo::ChannelMapPixelAlg::MaxROPs() const
  { return MaxPlanes(); }


// -----------------------------------------------------------------------------
bool geo::ChannelMapPixelAlg::HasROP(readout::ROPID const& ropid) const {
  return ropid && hasROP(ropid);
} // geo::ChannelMapPixelAlg::HasROP()


// -----------------------------------------------------------------------------
readout::ROPID geo::ChannelMapPixelAlg::WirePlaneToROP
  (geo::PlaneID const& planeid) const
{
  return planeid? convertPlaneToROP(planeid): readout::ROPID{};
} // geo::ChannelMapPixelAlg::WirePlaneToROP


// -----------------------------------------------------------------------------
std::vector<geo::PlaneID> geo::ChannelMapPixelAlg::ROPtoWirePlanes
  (readout::ROPID const& ropid) const
{
  std::vector<geo::PlaneID> planes;
  if (ropid) planes.push_back(convertROPtoPlane(ropid));
  return planes;
} // geo::ChannelMapPixelAlg::ROPtoWirePlanes()


// -----------------------------------------------------------------------------
std::vector<geo::TPCID> geo::ChannelMapPixelAlg::ROPtoTPCs
  (readout::ROPID const& ropid) const
{
  std::vector<geo::TPCID> TPCs;
  if (ropid) TPCs.push_back(convertROPtoPlane(ropid));
  return TPCs;
} // geo::ChannelMapPixelAlg::ROPtoTPCs()


// -----------------------------------------------------------------------------
    /// Returns the ID of the ROP the channel belongs to (invalid if none)
readout::ROPID geo::ChannelMapPixelAlg::ChannelToROP
  (raw::ChannelID_t channel) const
{
  assert(!fFirstChannelInfo.empty());
  auto const* channelInfo = fFirstChannelInfo.find(channel);
  return
    channelInfo? convertPlaneToROP(channelInfo->data().plane): readout::ROPID{};
} // geo::ChannelMapPixelAlg::ChannelToROP()


// -----------------------------------------------------------------------------
raw::ChannelID_t geo::ChannelMapPixelAlg::FirstChannelInROP
  (readout::ROPID const& ropid) const
{
  return ropid
    ? planeInfo(convertROPtoPlane(ropid)).firstChannel: raw::InvalidChannelID;
} // geo::ChannelMapPixelAlg::FirstChannelInROP()


// -----------------------------------------------------------------------------
geo::PlaneID geo::ChannelMapPixelAlg::FirstWirePlaneInROP
  (readout::ROPID const& ropid) const
{
  return ropid? convertROPtoPlane(ropid): geo::PlaneID{};
} // geo::ChannelMapPixelAlg::FirstWirePlaneInROP()


// -----------------------------------------------------------------------------
geo::GeoObjectSorter const& geo::ChannelMapPixelAlg::Sorter() const {
  throw geo::ChannelMapAlg::NoSorter("ChannelMapPixelAlg")
    << "ChannelMapPixelAlg does not provide any sorter."
    << "\nPlease configure one explicitly in the geometry service."
    << "\n";
} // geo::ChannelMapPixelAlg::Sorter()


// -----------------------------------------------------------------------------
geo::SigType_t geo::ChannelMapPixelAlg::SignalTypeForChannelImpl
  (raw::ChannelID_t const channel) const
{
  return raw::isValidChannelID(channel)? geo::kCollection: geo::kMysteryType;
} // geo::ChannelMapPixelAlg::SignalTypeForChannelImpl()


// -----------------------------------------------------------------------------
void geo::ChannelMapPixelAlg::fillFirstChannelInfo
  (geo::GeometryData_t::CryostatList_t const& cryostats)
{
  
  // we do not have geo::GeometryCore, so we need to iterate planes the long way
  assert(fFirstChannelInfo.empty());
  
  for (geo::CryostatGeo const& cryo: cryostats) {
    
    for (geo::TPCGeo const& tpc: cryo.IterateTPCs()) {
      
      for (geo::PlaneGeo const& plane: tpc.IteratePlanes()) {
        
        fFirstChannelInfo.append(
          {
            plane.ID() // plane
          },
          plane.Nwires()
          );
        
      } // for planes
      
    } // for TPCs
    
  } // for cryostats
  
  
} // geo::ChannelMapPixelAlg::fillFirstChannelInfo()


// -----------------------------------------------------------------------------
void geo::ChannelMapPixelAlg::fillGeometrySizes
  (geo::GeometryData_t::CryostatList_t const& cryostats)
{
  assert(!fGeometrySizes.has_value());
  
  // we do not have geo::GeometryCore, so we need to iterate planes the long way
  
  auto const [ NCryostats, MaxTPCs, MaxPlanes ]
    = geo::details::extractMaxGeometryElements<3U>(cryostats);
  
  fGeometrySizes.emplace(NCryostats, MaxTPCs, MaxPlanes);
  
  fGeometrySizes->nCryostats = cryostats.size();
  for (geo::CryostatGeo const& cryo: cryostats) {
    fGeometrySizes->nTPCs[cryo.ID().Cryostat] = cryo.NElements();
    
    for (geo::TPCGeo const& tpc: cryo.IterateTPCs()) {
      fGeometrySizes->nPlanes[tpc.ID()] = tpc.NElements();
      
      for (geo::PlaneGeo const& plane: tpc.IteratePlanes()) {
        fGeometrySizes->nWires[plane.ID()] = plane.NElements();
        
      } // for planes
      
    } // for TPCs
    
  } // for cryostats
  
  
} // geo::ChannelMapPixelAlg::fillGeometrySizes()


// -----------------------------------------------------------------------------
void geo::ChannelMapPixelAlg::makePlaneChannelMap() {
  
  /*
   * Inverting the existing map is slow, so we create one anew.
   * Assertion will ensure that the map is consistent.
   */
  assert(!fPlaneInfo.has_value());
  
  // FIXME use resizeAs()
  fPlaneInfo.emplace(NCryostats(), MaxTPCs(), MaxPlanes());
  
  raw::ChannelID_t firstChannel = 0U;
  for (auto c: util::counter<geo::CryostatID::CryostatID_t>(NCryostats())) {
    geo::CryostatID const cid { c };
    for (auto t: util::counter<geo::TPCID::TPCID_t>(NTPCs(cid))) {
      geo::TPCID const tid { cid, t };
      for (auto p: util::counter<geo::PlaneID::PlaneID_t>(NPlanes(tid))) {
        geo::PlaneID const pid { tid, p };
        
        assert(fFirstChannelInfo.findData(firstChannel)->plane == pid);
        
        fPlaneInfo.value()[pid].firstChannel = firstChannel;
        firstChannel += NWires(pid);
        
      } // for planes in TPC
      
    } // for TPCs in cryostat
    
  } // for cryostats
  assert((firstChannel - raw::ChannelID_t{0U}) == NChannels());
  
} // geo::ChannelMapPixelAlg::makePlaneChannelMap()


// -----------------------------------------------------------------------------
