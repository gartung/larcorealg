/**
 * @file   larcorealg/Geometry/GeometryBuilderStandard.cxx
 * @brief  Standard implementation of geometry extractor (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/Geometry/GeometryBuilderStandard.h`
 */

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderStandard.h"
#include "larcorealg/Geometry/WirePlaneGeo.h"
#include "larcorealg/Geometry/SenseWireGeo.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT libraries

// C++ standard library
#include <algorithm> // std::move()
#include <string_view>
#include <type_traits> // std::is_convertible_v


namespace {

  //----------------------------------------------------------------------------
  template <typename Dest, typename Src>
  Dest& extendCollection(Dest& dest, Src&& src) {
    std::move(src.begin(), src.end(), std::back_inserter(dest));
    return dest;
  } // extend()


  //----------------------------------------------------------------------------
  template <typename Dest, typename Src>
  std::unique_ptr<Dest> into_unique_ptr(Src&& src) {
    using DestPtr_t = std::unique_ptr<Dest>;
    
    if constexpr (std::is_convertible_v<std::decay_t<decltype(src)>, DestPtr_t>)
      return std::move(src);
    else
      return std::make_unique<Dest>(std::move(src));
  } // into_unique_ptr()
  
  
  //----------------------------------------------------------------------------

} // local namespace



//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::GeometryBuilderStandard(Config const& config)
  : fMaxDepth(config.maxDepth())
  , fCryostatPattern(config.cryostatPattern()) // default flags
  , fTPCPattern(config.TPCPattern()) // default flags
  , fPlanePattern(config.planePattern()) // default flags
  , fSensElemPattern(config.sensElemPattern()) // default flags
  , fAuxDetPattern(config.auxDetPattern()) // default flags
  , fAuxDetSensPattern(config.auxDetSensPattern()) // default flags
  , fOpDetPattern(config.opDetPattern()) // default flags
{
  MF_LOG_DEBUG("GeometryBuilder")
    << "Loading geometry builder: GeometryBuilderStandard"
    << "\nGeometry element search patterns: "
    << "\n - aux.det.: '" << config.auxDetPattern() << "'"
    << "\n -   sens. : '" << config.auxDetSensPattern() << "'"
    << "\n - cryostat: '" << config.cryostatPattern() << "'"
    << "\n - op.det.:  '" << config.opDetPattern() << "'"
    << "\n - TPC:      '" << config.TPCPattern() << "'"
    << "\n - plane:    '" << config.planePattern() << "'"
    << "\n - sens.el.: '" << config.sensElemPattern() << "'"
    ;
} // geo::GeometryBuilderStandard::GeometryBuilderStandard()


//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::AuxDets_t
geo::GeometryBuilderStandard::doExtractAuxiliaryDetectors(Path_t& path) {

  return doExtractGeometryObjects<
    geo::AuxDetGeo,
    &geo::GeometryBuilderStandard::isAuxDetNode,
    &geo::GeometryBuilderStandard::doMakeAuxDet
    >
    (path, "auxiliary detector");

} // geo::GeometryBuilderStandard::doExtractAuxiliaryDetectors()


//------------------------------------------------------------------------------
geo::AuxDetGeo geo::GeometryBuilderStandard::doMakeAuxDet(Path_t& path) {

  return geo::AuxDetGeo(
    path.current(), path.currentTransformation<geo::TransformationMatrix>(),
    geo::GeometryBuilder::moveToColl(extractAuxDetSensitive(path))
    );

} // geo::GeometryBuilderStandard::doMakeAuxDet()


//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::AuxDetSensitive_t
geo::GeometryBuilderStandard::doExtractAuxDetSensitive(Path_t& path) {
  return doExtractGeometryObjects<
    geo::AuxDetSensitiveGeo,
    &geo::GeometryBuilderStandard::isAuxDetSensitiveNode,
    &geo::GeometryBuilderStandard::makeAuxDetSensitive
    >
    (path, "sensitive auxiliary detector");
} // geo::GeometryBuilderStandard::doExtractAuxDetSensitive()


//------------------------------------------------------------------------------
geo::AuxDetSensitiveGeo geo::GeometryBuilderStandard::doMakeAuxDetSensitive
  (Path_t& path)
{
  return geo::AuxDetSensitiveGeo
    (path.current(), path.currentTransformation<geo::TransformationMatrix>());
} // geo::GeometryBuilderStandard::doMakeAuxDetSensitive()


//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::Cryostats_t
geo::GeometryBuilderStandard::doExtractCryostats(Path_t& path) {

  return doExtractGeometryObjects<
    geo::CryostatGeo,
    &geo::GeometryBuilderStandard::isCryostatNode,
    &geo::GeometryBuilderStandard::makeCryostat
    >
    (path, "cryostat", 1U); // require at least 1

} // geo::GeometryBuilderStandard::doExtractCryostats()


//------------------------------------------------------------------------------
geo::CryostatGeo geo::GeometryBuilderStandard::doMakeCryostat(Path_t& path) {

  return geo::CryostatGeo(
    path.current(), path.currentTransformation<geo::TransformationMatrix>(),
    geo::GeometryBuilder::moveToColl(extractTPCs(path)),
    geo::GeometryBuilder::moveToColl(extractOpDets(path))
    );

} // geo::GeometryBuilderStandard::doMakeCryostat()


//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::OpDets_t
geo::GeometryBuilderStandard::doExtractOpDets(Path_t& path) {
  return doExtractGeometryObjects<
    geo::OpDetGeo,
    &geo::GeometryBuilderStandard::isOpDetNode,
    &geo::GeometryBuilderStandard::makeOpDet
    >
    (path, "optical detector");
} // geo::GeometryBuilderStandard::doExtractOpDets()


//------------------------------------------------------------------------------
geo::OpDetGeo geo::GeometryBuilderStandard::doMakeOpDet(Path_t& path) {
  return geo::OpDetGeo
    (path.current(), path.currentTransformation<geo::TransformationMatrix>());
} // geo::GeometryBuilderStandard::doMakeOpDet()


//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::TPCs_t geo::GeometryBuilderStandard::doExtractTPCs
  (Path_t& path)
{
  return doExtractGeometryObjects<
    geo::TPCGeo,
    &geo::GeometryBuilderStandard::isTPCNode,
    &geo::GeometryBuilderStandard::makeTPC
    >
    (path, "TPC", 1U); // require at least 1

} // geo::GeometryBuilderStandard::doExtractTPCs()


//------------------------------------------------------------------------------
geo::TPCGeo geo::GeometryBuilderStandard::doMakeTPC(Path_t& path) {
  return geo::TPCGeo(
    path.current(), path.currentTransformation<geo::TransformationMatrix>(),
    extractPlanes(path)
    );
} // geo::GeometryBuilderStandard::doMakeTPC()


//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::Planes_t
geo::GeometryBuilderStandard::doExtractPlanes(Path_t& path)
{
  return doExtractGeometryObjects<
    geo::PlaneGeo, PlanePtr_t,
    &geo::GeometryBuilderStandard::isPlaneNode,
    &geo::GeometryBuilderStandard::makePlane
    >
    (path, "anode plane", 1U); // require at least 1

} // geo::GeometryBuilderStandard::doExtractPlanes()


//------------------------------------------------------------------------------
auto geo::GeometryBuilderStandard::doMakePlane(Path_t& path) -> PlanePtr_t {
  return std::make_unique<geo::WirePlaneGeo>(
    path.current(), path.currentTransformation<geo::TransformationMatrix>(),
    extractWires(path)
    );
} // geo::GeometryBuilderStandard::doMakePlane()


//------------------------------------------------------------------------------
geo::GeometryBuilderStandard::Wires_t
geo::GeometryBuilderStandard::doExtractWires(Path_t& path)
{
  return doExtractGeometryObjects<
    geo::SenseWireGeo,
    &geo::GeometryBuilderStandard::isWireNode,
    &geo::GeometryBuilderStandard::makeWire
    >
    (path, "sensitive element", 2U); // require at least 2
    
} // geo::GeometryBuilderStandard::doExtractWires()


//------------------------------------------------------------------------------
geo::SenseWireGeo geo::GeometryBuilderStandard::doMakeWire(Path_t& path) {

  return
    { path.current(), path.currentTransformation<geo::TransformationMatrix>() };

} // geo::GeometryBuilderStandard::doMakeWire()


//------------------------------------------------------------------------------
bool geo::GeometryBuilderStandard::isAuxDetNode(TGeoNode const& node) const {
  return std::regex_search(node.GetName(), fAuxDetPattern);
} // geo::GeometryBuilderStandard::isAuxDetNode()


//------------------------------------------------------------------------------
bool geo::GeometryBuilderStandard::isAuxDetSensitiveNode
  (TGeoNode const& node) const
{
  return std::regex_search(node.GetName(), fAuxDetSensPattern);
} // geo::GeometryBuilderStandard::isAuxDetSensitiveNode()


//------------------------------------------------------------------------------
bool geo::GeometryBuilderStandard::isCryostatNode(TGeoNode const& node) const {
  return std::regex_search(node.GetName(), fCryostatPattern);
} // geo::GeometryBuilderStandard::isCryostatNode()


//------------------------------------------------------------------------------
bool geo::GeometryBuilderStandard::isOpDetNode(TGeoNode const& node) const {
  return std::regex_search(node.GetName(), fOpDetPattern);
} // geo::GeometryBuilderStandard::isOpDetNode()



//------------------------------------------------------------------------------
bool geo::GeometryBuilderStandard::isTPCNode(TGeoNode const& node) const {
  return std::regex_search(node.GetName(), fTPCPattern);
} // geo::GeometryBuilderStandard::isTPCNode()


//------------------------------------------------------------------------------
bool geo::GeometryBuilderStandard::isPlaneNode(TGeoNode const& node) const {
  return std::regex_search(node.GetName(), fPlanePattern);
} // geo::GeometryBuilderStandard::isPlaneNode()


//------------------------------------------------------------------------------
bool geo::GeometryBuilderStandard::isWireNode(TGeoNode const& node) const {
  return std::regex_search(node.GetName(), fSensElemPattern);
} // geo::GeometryBuilderStandard::isWireNode()


//------------------------------------------------------------------------------
template <
  typename ObjGeoIF, typename ObjGeo,
  bool (geo::GeometryBuilderStandard::*IsObj)(TGeoNode const&) const,
  ObjGeo (geo::GeometryBuilderStandard::*MakeObj)(geo::GeometryBuilder::Path_t&)
  >
geo::GeometryBuilder::GeoPtrColl_t<ObjGeoIF>
geo::GeometryBuilderStandard::doExtractGeometryObjects(
  Path_t& path,
  std::string const& objName, /* = "" */
  unsigned int min /* = 0U */
) {

  using ObjColl_t = geo::GeometryBuilder::GeoPtrColl_t<ObjGeoIF>;
  
  ObjColl_t objs;
  
  // for debugging purposes, uncomment this line:
//  MF_LOG_TRACE("GeometryBuilder")
//    << "Looking for " << objName << " in " << std::string(path);

  //
  // if this is a target object, we are set
  //
  if ((this->*IsObj)(path.current())) {
    // trying to be very very smart here:
    //  * make a unique pointer of the made object,
    //  * but only if it is not a unique pointer already
    objs.emplace_back(::into_unique_ptr<ObjGeoIF>((this->*MakeObj)(path)));
    return objs;
  }

  //
  // descend into the next layer down, concatenate the results and return them
  //
  if (path.depth() >= fMaxDepth) return objs; // yep, this is empty

  TGeoVolume const& volume = *(path.current().GetVolume());
  int const n = volume.GetNdaughters();
  for (int i = 0; i < n; ++i) {
    path.append(*(volume.GetNode(i)));
    extendCollection(
      objs,
      doExtractGeometryObjects<ObjGeoIF, ObjGeo, IsObj, MakeObj>
        (path, objName, 0U)
      );
    path.pop();
  } // for
  
  if (objs.size() < min) {
    throw cet::exception("GeometryBuilder")
      << "Collected " << objs.size() << " " << objName
      << " object nodes under path " << std::string(path)
      << " (at least " << min << " required)\n";
  }
  else if (!objs.empty()) { // feels still too noisy
    MF_LOG_TRACE("GeometryBuilder")
      << "Found " << objs.size() << " " << objName << " object nodes under "
      << std::string(path);
  }
  return objs;

} // geo::GeometryBuilderStandard::doExtractGeometryObjects()


//------------------------------------------------------------------------------

