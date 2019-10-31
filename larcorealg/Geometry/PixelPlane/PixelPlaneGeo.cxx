/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeo.cxx
 * @brief   Representation of a readout plane with pixels as sensitive elements.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 23, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`
 * @ingroup Geometry
 */

// class header
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h"

// LArSoft includes
#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "larcorealg/Geometry/SimpleGeo.h" // lar::util::simple_geo::Rectangle
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::rounded01()
#include "larcorealg/CoreUtils/RealComparisons.h" // makeVector3DComparison()
#include "larcorealg/CoreUtils/zip.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Zaxis()

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard library
#include <sstream> // std::ostringstream
#include <utility> // std::move()
#include <algorithm> // std::clamp(), std::swap()
#include <array>
#include <limits> // std::numeric_limits<>
#include <cmath> // std::cos(), std::abs(), ...
#include <cassert>


// -----------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  template <typename T>
  inline bool closeToUnity(T const value, T const tol = 1e-5) {
    
    return std::abs(value - T{ 1 }) <= tol;
    
  } // closeToUnity()
  
  
  // ---------------------------------------------------------------------------
  
} // local namespace


// -----------------------------------------------------------------------------
// --- geo::PixelPlaneGeo
// -----------------------------------------------------------------------------
std::string const geo::PixelPlaneGeo::SensitivePixelPlaneVolumeName {}; // TODO


// -----------------------------------------------------------------------------
geo::PixelPlaneGeo::PixelPlaneGeo(
  TGeoNode const& node,
  geo::TransformationMatrix&& trans
  )
  : geo::PlaneGeo(node, std::move(trans))
{
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Plane extends " << Width() << " cm in " << WidthDir<geo::Vector_t>()
    << " and " << Depth() << " cm in " << DepthDir<geo::Vector_t>();
  
  SetView(geo::k3D); // view is this simple
  
  // TODO some updates should be moved here
  
} // geo::PixelPlaneGeo::PixelPlaneGeo()


// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeo::PixelInfo(
  PixelCoordID_t const coords,
  std::string indent /* = "" */, unsigned int verbosity /* = 1U */
  ) const
{
  std::ostringstream sstr;
  PrintPixelInfo(sstr, coords, indent, verbosity);
  return sstr.str();
} // geo::PixelPlaneGeo::PixelInfo()


// -----------------------------------------------------------------------------
// ---  interface implementation: anode plane
// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doThetaZ() const {
  return fThetaZ;
} // geo::PixelPlaneGeo::doThetaZ()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doPhiZ() const {
  return fPhiZ;
} // geo::PixelPlaneGeo::doPhiZ()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doSinPhiZ() const {
  // if these need to be cached... can do
  return std::sin(fPhiZ);
} // geo::PixelPlaneGeo::doSinPhiZ()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doCosPhiZ() const {
  // if these need to be cached... can do
  return std::cos(fPhiZ);
} // geo::PixelPlaneGeo::doCosPhiZ()


// -----------------------------------------------------------------------------
geo::WireGeo geo::PixelPlaneGeo::doWire(unsigned int iWire) const {
  return getWire(iWire);
} // geo::PixelPlaneGeo::doWire()


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeo::doNwires() const {
  
  return getNsensElem();
  
} // geo::PixelPlaneGeo::doNwires()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWirePitch() const
  { return getSensElemPitch(geo::pixel::ixWireC); }


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::doWireIDincreasesWithZ() const {
  assert(getSensElemDir(geo::pixel::ixWireC).Dot(geo::Zaxis()) > 0.0);
  return true; // pretty much by the definition of the wire ID's
} // geo::PixelPlaneGeo::doWireIDincreasesWithZ()


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::doGetIncreasingWireDirection() const
  { return getSensElemDir(geo::pixel::ixWireC); }


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::doGetWireDirection() const
  { return getSensElemDir(geo::pixel::ixWireD); }


// -----------------------------------------------------------------------------
geo::WireID geo::PixelPlaneGeo::doNearestWireID(geo::Point_t const& pos) const {
  
  using namespace geo::pixel;
  
  PixelCoordID_t nearest = coordsOf(getPlaneCoordinates(pos));
  if (isPixelOnPlane(nearest)) return { ID(), indexOf(nearest) };
  
  // of, something went bad, and we have to implement a crazy protocol.
  PixelCoordID_t const onSpot = nearest;
  
  nearest.main()
    = std::clamp<PixelCoordIndex_t>(nearest.main(), 0, getNsensElem(ixMain));
  nearest.secondary() =
    std::clamp<PixelCoordIndex_t>(nearest.secondary(), 0, getNsensElem(ixSec));
  
  PixelIndex_t const nearestIndex = indexOf(nearest);
  throw InvalidWireError("Geometry", ID(), nearestIndex, InvalidPixelIndex)
    << "Nearest pixel for position " << pos << " in plane " << std::string(ID())
    << " is approx number #" << nearestIndex << " ("
    << nearest.main() << "; " << nearest.secondary() << "), capped from ("
    << onSpot.main() << "; " << onSpot.secondary()  << ")\n";
  
} // geo::PixelPlaneGeo::doNearestWireID()


// -----------------------------------------------------------------------------
geo::WireGeo geo::PixelPlaneGeo::doNearestWire(geo::Point_t const& pos) const {
  
  return getWire(geo::PixelPlaneGeo::doNearestWireID(pos).Wire);
  
} // geo::PixelPlaneGeo::doNearestWire()


// -----------------------------------------------------------------------------
geo::WireID geo::PixelPlaneGeo::doClosestWireID(geo::WireID::WireID_t) const {
  //
  // The problem is that a `geo::WireID` which is invalid does not correspond to
  // a well defined non-existing pixel, because the mapping 2D pixel coordinate
  // to wire number is defined only in the validity range and this definition
  // can't be univocally extended beyond it (well, it can, but it requires
  // conventions that are not satisfactory).
  // 
  // Ok, what now? The calling code might be rewritten to use concepts better
  // suited for a 2D plane. Chances are that this feature that code requires
  // needs to be implemented and the `geo::PlaneGeo` interface extended.
  // Open a ticket!
  //
  
  NotImplemented("This concept is not appropriate to pixel planes.");
  
} // geo::PixelPlaneGeo::doClosestWireID()


// -----------------------------------------------------------------------------
lar::util::simple_geo::Volume<> geo::PixelPlaneGeo::doCoverage() const {
  
  using namespace geo::pixel;
  
  std::array<double, 3U> A, B;
  geo::vect::fillCoords(A,
    fCenter
      - getNsensElem(ixMain) * getSensElemHalfStepDir(ixMain)
      - getNsensElem(ixSec) * getSensElemHalfStepDir(ixSec)
    );
  geo::vect::fillCoords(B,
    fCenter
      + getNsensElem(ixMain) * getSensElemHalfStepDir(ixMain)
      + getNsensElem(ixSec) * getSensElemHalfStepDir(ixSec)
    );
  
  // not sure whether it is a problem if the two points are not sorted...
  for (auto&& [ a, b ]: util::zip(A, B))
    if (a > b) std::swap(a, b);
  
  return { A.data(), B.data() };
} // geo::PixelPlaneGeo::WirePlaneGeo::doCoverage()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doPlaneCoordinateFrom
  (geo::Point_t const& point, geo::WireGeo const& refPixel) const
{
  return getPlaneCoordinateFrom(point, refPixel, geo::pixel::ixWireC);
} // geo::PixelPlaneGeo::doPlaneCoordinateFrom()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doPlaneCoordinate(geo::Point_t const& point) const {
  return getPlaneCoordinate(point, geo::pixel::ixWireC);
} // geo::PixelPlaneGeo::doPlaneCoordinate()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireCoordinate(geo::Point_t const& point) const {
  return getWireCoordinate(point, geo::pixel::ixWireC);
} // geo::PixelPlaneGeo::doWireCoordinate()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::doDecomposePoint(geo::Point_t const& point) const
  -> WireDecomposedVector_t
{
  return fDecompPixel.DecomposePoint(point);
} // geo::PixelPlaneGeo::doDecomposePoint()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doProjectionReferencePoint() const {
  return fDecompPixel.ReferencePoint();
} // geo::PixelPlaneGeo::doProjectionReferencePoint()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::doProjection(geo::Point_t const& point) const
  -> WireCoordProjection_t
{
  return getPlaneCoordinates(point);
} // geo::PixelPlaneGeo::doProjection()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::doProjection(geo::Vector_t const& v) const
  -> WireCoordProjection_t
{
  return fDecompPixel.ProjectVectorOnPlane(v);
} // geo::PixelPlaneGeo::doProjection()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doComposePoint
  (WireDecomposedVector_t const& decomp) const
{
  return fDecompPixel.ComposePoint(decomp); 
} // geo::PixelPlaneGeo::doComposePoint(WireDecomposedVector_t)


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doComposePoint
  (double distance, WireCoordProjection_t const& proj) const
{
  return fDecompPixel.ComposePoint(distance, proj); 
} // geo::PixelPlaneGeo::doComposePoint(double, WireCoordProjection_t)


// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeo::doPlaneInfo
  (std::string indent /* = "" */, unsigned int verbosity /* = 1U */) const
{
  std::ostringstream out;
  PrintPixelPlaneInfo(out, indent, verbosity);
  return out.str();
} // geo::PixelPlaneGeo::doPlaneInfo()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::doUpdateAfterSorting
  (geo::BoxBoundedGeo const& /* TPCbox */)
{
  
  UpdateDecompPixel();
  
  UpdateNpixelsAndPitches();
  
  UpdatePixelDirs();
  
  UpdateDecompPixelOrigin();
  
  UpdatePlaneCenter();
  
  UpdateDirections();
  
  UpdateActiveArea();
  
} // geo::PixelPlaneGeo::doUpdateAfterSorting()


// -----------------------------------------------------------------------------
// --- Polymorphic implementation: wire abstraction
// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireRMax(WireLocator const&) const {
  using namespace geo::pixel;
  return std::max(getSensElemPitch(ixMain), getSensElemPitch(ixSec)) / 2.0;
} // geo::PixelPlaneGeo::doWireRMax()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireRMin(WireLocator const&) const {
  using namespace geo::pixel;
  return std::min(getSensElemPitch(ixMain), getSensElemPitch(ixSec)) / 2.0;
} // geo::PixelPlaneGeo::doWireRMin()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireHalfL(WireLocator const& wloc) const {
  return getPixelHalfL(wloc);
} // geo::PixelPlaneGeo::doWireHalfL()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::doWireFillCenterXYZ
  (WireLocator const& wloc, double* xyz, double localz /* = 0.0 */) const
{
  geo::vect::fillCoords(xyz, doWireGetPositionFromCenter(wloc, localz));
} // geo::PixelPlaneGeo::doWireFillCenterXYZ()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::doWireFillStartXYZ
  (WireLocator const& wloc, double* xyz) const
{
  geo::vect::fillCoords(xyz, doWireGetStart(wloc));
} // geo::PixelPlaneGeo::doWireFillStartXYZ()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::doWireFillEndXYZ
  (WireLocator const& wloc, double* xyz) const
{
  geo::vect::fillCoords(xyz, doWireGetEnd(wloc));
} // geo::PixelPlaneGeo::doWireFillEndXYZ()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doWireGetPositionFromCenter
  (WireLocator const& wloc, double localz) const
{
  return
    doWireGetPositionFromCenterUnbounded(wloc, std::clamp(localz, -1.0, +1.0));
} // geo::PixelPlaneGeo::doWireGetPositionFromCenter()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doWireGetPositionFromCenterUnbounded
  (WireLocator const& wloc, double localz) const
{
  return getSensElemCenter(coordsOf(wloc))
    + localz * getSensElemHalfStepDir(geo::pixel::ixWireD);
} // geo::PixelPlaneGeo::doWireGetPositionFromCenterUnbounded()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doWireGetCenter(WireLocator const& wloc) const
{
  return getSensElemCenter(coordsOf(wloc));
} // geo::PixelPlaneGeo::doWireGetCenter()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doWireGetStart(WireLocator const& wloc) const
{
  using namespace geo::pixel;
  return getSensElemCenter(coordsOf(wloc))
    - getSensElemHalfStepDir(ixWireC)
    - getSensElemHalfStepDir(ixWireD)
    ;
} // geo::PixelPlaneGeo::doWireGetStart()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::doWireGetEnd(WireLocator const& wloc) const
{
  using namespace geo::pixel;
  return getSensElemCenter(coordsOf(wloc))
    + getSensElemHalfStepDir(ixWireC)
    + getSensElemHalfStepDir(ixWireD)
    ;
} // geo::PixelPlaneGeo::doWireGetEnd()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireLength(WireLocator const& wloc) const {
  return 2.0 * getSensElemPitch(geo::pixel::ixWireD);
} // geo::PixelPlaneGeo::doWireLength()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireThetaZ(WireLocator const&) const {
  return ThetaZ();
} // geo::PixelPlaneGeo::doWireThetaZ()


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::doWireIsParallelTo
  (WireLocator const&, geo::WireGeo const& wire) const
{
  return closeToUnity
    (getSensElemDir(geo::pixel::ixWireD).Dot(wire.Direction<geo::Vector_t>()));
} // geo::PixelPlaneGeo::doWireIsParallelTo()


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::doWireDirection(WireLocator const&) const {
  return getSensElemDir(geo::pixel::ixWireD);
} // geo::PixelPlaneGeo::doWireDirection()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireDistanceFrom
  (WireLocator const& wloc, geo::WireGeo const& wire) const
{
  // can't use getSecPlaneCoordinate() because it uses pixel #0 reference
  return fDecompPixel.VectorSecondaryComponent
    (getSensElemCenter(coordsOf(wloc)) - wire.GetCenter<geo::Point_t>());
} // geo::PixelPlaneGeo::doWireDistanceFrom()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doWireComputeZatY0(WireLocator const& wloc) const {
  
  //
  // the prescription of this method is complicate...
  // but given that the main direction (the one we extrapolate on) is aligned
  // with depth direction and that the width and secondary directions coincide,
  // this pretty much turns into the width coordinate of the center of the
  // pixel.
  // We could push the assumption further and reduce all to a shift...
  // not bothering to, though.
  //
  
  return fDecompFrame.PointMainComponent(getSensElemCenter(coordsOf(wloc)));
  
} // geo::PixelPlaneGeo::doWireComputeZatY0()


// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeo::doWireInfo(
  WireLocator const& wloc,
  std::string indent /* = "" */, unsigned int verbosity /* = 1 */
  ) const
{
  return PixelInfo(coordsOf(wloc), indent, verbosity);
} // geo::PixelPlaneGeo::doWireInfo()


// -----------------------------------------------------------------------------
// --- Plane coordinate and ID conversions
// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::indexAt(WireCoordProjection_t const& point) const
  -> PixelIndex_t
{
  return indexOf(coordsOf(point));
} // geo::PixelPlaneGeo::indexAt()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::coordsOf(WireCoordProjection_t const& point) const
  -> PixelCoordID_t
{
  using namespace geo::pixel;
  return { roundCoord(point.X(), ixMain), roundCoord(point.Y(), ixSec) };
} // geo::PixelPlaneGeo::coordsOf(WireCoordProjection_t)


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::coordsOf(WireLocator const& wloc) const
  -> PixelCoordID_t
{
  return coordsOf(wireToPixelIndex(wloc));
} // geo::PixelPlaneGeo::coordsOf(WireLocator)


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::coordsOf(PixelIndex_t const index) const
  -> PixelCoordID_t
{
  //
  // indexOf():
  //   coords.secondary() * getNsensElem(ixMain) + coords.main()
  //   
  
  using geo::pixel::ixMain;
  
  return {
    static_cast<PixelCoordIndex_t>(index % getNsensElem(ixMain)), // main
    static_cast<PixelCoordIndex_t>(index / getNsensElem(ixMain))  // secondary
    };
} // geo::PixelPlaneGeo::coordsOf()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::roundCoord(double coord, DirIndex_t const dir) const
  -> PixelCoordIndex_t
{
  return roundPixelCoord(pixelCoord(coord, dir));
} // geo::PixelPlaneGeo::roundCoord()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::pixelCoord(double coord, DirIndex_t const dir) const
{
  assert(isDirIndex(dir));
  return coord / getSensElemPitch(dir);
} // geo::PixelPlaneGeo::pixelCoord()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::roundPixelCoord(double const coord) const
  -> PixelCoordIndex_t
{
  return static_cast<PixelCoordIndex_t>(coord + 0.5);
} // geo::PixelPlaneGeo::roundPixelCoord()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::indexOf(PixelCoordID_t const& coords) const
  -> PixelIndex_t
{
  using geo::pixel::ixMain;
  return coords.secondary() * getNsensElem(ixMain) + coords.main();
} // geo::PixelPlaneGeo::indexOf()


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isOnPlane(WireCoordProjection_t const& point) const {
  using namespace geo::pixel;
  return isOnPlane(point.X(), ixMain) && isOnPlane(point.Y(), ixSec);
} // geo::PixelPlaneGeo::isOnPlane(WireCoordProjection_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isOnPlane
  (double const coord, DirIndex_t const dir) const
{
  // yeah, strictly speaking it should be `<` instead of `<=`. Judge me.
  double const shiftedCoord = coord + getSensElemPitch(dir) / 2.0;
  return (shiftedCoord >= 0.0) && (shiftedCoord <= getSensElemDirSize(dir));
} // geo::PixelPlaneGeo::isOnPlane(double)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isPixelOnPlane(PixelCoords_t const& coords) const {
  using namespace geo::pixel;
  return isPixelOnPlane(coords.main(), ixMain)
    && isPixelOnPlane(coords.secondary(), ixSec);
} // bool geo::PixelPlaneGeo::isOnPlane(PixelCoords_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isPixelOnPlane(PixelCoordID_t const& coords) const {
  using namespace geo::pixel;
  return isPixelOnPlane(coords.main(), ixMain)
    && isPixelOnPlane(coords.secondary(), ixSec);
} // geo::PixelPlaneGeo::isOnPlane(PixelCoordID_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isPixelOnPlane
  (double const coord, DirIndex_t const dir) const
{
  double const shiftedCoord = 0.5 + coord;
  return (shiftedCoord >= 0.0)
    && (shiftedCoord < static_cast<double>(getNsensElem(dir)));
} // geo::PixelPlaneGeo::isPixelOnPlane(double)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isPixelOnPlane
  (PixelCoordIndex_t const coord, DirIndex_t const dir) const
{
  return
    (coord >= 0) && (coord < static_cast<PixelCoordIndex_t>(getNsensElem(dir)));
} // geo::PixelPlaneGeo::isPixelOnPlane(PixelCoordIndex_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isPixelOnPlane(PixelIndex_t const index) const {
  return index < getNsensElem();
} // geo::PixelPlaneGeo::isOnPlane(PixelIndex_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isWireIDvalid(geo::WireID const& wireid) const {
  return isWireIDvalid(wireid.Wire);
} // geo::PixelPlaneGeo::isWireIDvalid(WireID)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeo::isWireIDvalid(geo::WireID::WireID_t const wire) const {
  return (wire >= 0) && (wire < getNsensElem());
} // geo::PixelPlaneGeo::isWireIDvalid(WireID_t)


// -----------------------------------------------------------------------------
// --- Candidate extensions to `geo::PlaneGeo` interface
// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::doSensElemDir(DirIndex_t const dir) const
  { return getSensElemDir(dir); }


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::doSensElemMainDir() const
  { return getSensElemDir(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::doSensElemSecondaryDir() const
  { return getSensElemDir(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doSensElemDirSize(DirIndex_t const dir) const
  { return getSensElemDirSize(dir); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doSensElemMainDirSize() const
  { return getSensElemDirSize(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doSensElemSecondaryDirSize() const
  { return getSensElemDirSize(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeo::doNsensElem(DirIndex_t const dir) const
  { return getNsensElem(dir); }


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeo::doNsensElemMain() const
  { return getNsensElem(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeo::doNsensElemSecondary() const
  { return getNsensElem(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doSensElemPitch(DirIndex_t const dir) const
  { return getSensElemPitch(dir); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doSensElemPitchMain() const
  { return getSensElemPitch(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::doSensElemPitchSecondary() const
  { return getSensElemPitch(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
geo::WireGeo geo::PixelPlaneGeo::getWire(PixelIndex_t const iWire) const {
  return { *this, { ID(), static_cast<geo::WireID::WireID_t>(iWire) } };
} // geo::PixelPlaneGeo::getWire()


// -----------------------------------------------------------------------------
// --- Implementation of the candidate interface extensions
// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::getSensElemDir(DirIndex_t const dir) const {
  assert(isDirIndex(dir));
  switch (dir) {
    case geo::pixel::ixMain: return fDecompPixel.MainDir();
    case geo::pixel::ixSec:  return fDecompPixel.SecondaryDir();
    default:     return {};
  } // switch
} // geo::PixelPlaneGeo::getSensElemDir()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getSensElemDirSize(DirIndex_t const dir) const {
  //
  // this code assumes that the main direction matches the frame width direction
  // and the secondary direction msatches the frame depth direction
  //
  assert(isDirIndex(dir));
  switch (dir) {
    case geo::pixel::ixMain:
      assert(closeToUnity(std::abs(DepthDir<geo::Vector_t>().Dot(getSensElemDir(dir)))));
      return fFrameSize.Depth();
    case geo::pixel::ixSec:
      assert(closeToUnity(std::abs(WidthDir<geo::Vector_t>().Dot(getSensElemDir(dir)))));
      return fFrameSize.Width();
    default:
      return 0;
  } // switch
} // geo::PixelPlaneGeo::getSensElemDirSize()


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeo::getSensElemHalfStepDir
  (DirIndex_t const dir) const
{
  assert(isDirIndex(dir));
  return fPixelDirs[dir];
} // geo::PixelPlaneGeo::getSensElemHalfStepDir()


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeo::getNsensElem(DirIndex_t const dir) const {
  assert(isDirIndex(dir));
  return fNPixels[dir];
} // geo::PixelPlaneGeo::getNsensElem(DirIndex_t)


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeo::getNsensElem() const {
  using namespace geo::pixel;
  return getNsensElem(ixMain) * getNsensElem(ixSec);
} // geo::PixelPlaneGeo::getNsensElem()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getSensElemPitch(DirIndex_t const dir) const {
  assert(isDirIndex(dir));
  return fPitches[dir];
} // geo::PixelPlaneGeo::getSensElemPitch()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getPlaneCoordinateFrom
  (geo::Point_t const& point, geo::WireGeo const& ref, DirIndex_t const dir)
  const
{
  return getPlaneCoordinate(point - ref.GetCenter<geo::Vector_t>(), dir);
} // geo::PixelPlaneGeo::getPlaneCoordinateFrom()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeo::getPlaneCoordinates(geo::Point_t const& point) const
  -> WireCoordProjection_t
{
  return fDecompPixel.ProjectPointOnPlane(point);
} // geo::PixelPlaneGeo::getPlaneCoordinates()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getPlaneCoordinate
  (geo::Point_t const& point, DirIndex_t const dir) const
{
  // this is just a dispatcher since the function methods are different
  assert(isDirIndex(dir));
  switch (dir) {
    case geo::pixel::ixMain:
      return getMainPlaneCoordinate(point);
    case geo::pixel::ixSec:
      return getSecPlaneCoordinate(point);
    default: // we claimed it's undefined behavior, but actually it secretly is:
      return std::numeric_limits<double>::signaling_NaN();
  } // switch
  
} // geo::PixelPlaneGeo::getPlaneCoordinate()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getMainPlaneCoordinate
  (geo::Point_t const& point) const
{
  return fDecompPixel.PointMainComponent(point);
} // geo::PixelPlaneGeo::getMainPlaneCoordinateFrom()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getSecPlaneCoordinate
  (geo::Point_t const& point) const
{
  return fDecompPixel.PointSecondaryComponent(point);
} // geo::PixelPlaneGeo::getSecPlaneCoordinate()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getWireCoordinate
  (geo::Point_t const& point, DirIndex_t const dir) const
{
  // no rounding, no shifting; pixel #0 is covering coordinates from -0.5 to 0.5
  // in pitch units
  return getPlaneCoordinate(point, dir) / getSensElemPitch(dir);
} // geo::PixelPlaneGeo::getWireCoordinate()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::getSensElemCenter
  (PixelCoordID_t const& coords) const
{
  using namespace geo::pixel;
  return firstPixelCenter()
    + getSensElemHalfStepDir(ixMain) * (2.0 * coords[ixMain])
    + getSensElemHalfStepDir(ixSec) * (2.0 * coords[ixSec])
    ;
} // geo::PixelPlaneGeo::getSensElemCenter()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::getPixelHalfL(WireLocator const&) const {
  return getSensElemPitch(geo::pixel::ixWireD) / 2.0;
} // geo::PixelPlaneGeo::getPixelHalfL()


// -----------------------------------------------------------------------------
// --- implementation details (private methods)
// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdateNpixelsAndPitches() {
  
  //
  // requirements:
  //   * `Width()` and `Depth()` correct
  //
  
  discoverPitches();
  
  for (DirIndex_t dir = 0; dir < geo::pixel::NCoords; ++dir)
    fNPixels[dir] = getSensElemDirSize(dir) / fPitches[dir];
  
} // geo::PixelPlaneGeo::UpdateNpixelsAndPitches()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdateDecompPixel() {
  
  //
  // requires:
  //  * width and depth directions (including their verse)
  //
  
  //
  // direction measured by the secondary coordinate; it is the width direction
  //
  fDecompPixel.SetSecondaryDir
    (geo::vect::rounded01(WidthDir<geo::Vector_t>(), 1e-4));

  //
  // get the axis perpendicular to it on the wire plane
  // (verse is already set to have the base as positive as the frame base is)
  //
  fDecompPixel.SetMainDir
    (geo::vect::rounded01(-DepthDir<geo::Vector_t>(), 1e-4));
  
  //
  // check that the resulting normal matches the plane one
  //
  assert(lar::util::makeVector3DComparison(1e-5)
    .equal(fDecompPixel.NormalDir(), GetNormalDirection<geo::Vector_t>()));
  
} // geo::PixelPlaneGeo::UpdateDecompPixel()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdatePixelDirs() {
  //
  // requirements:
  // * "wire" frame direction definitions completely finalized
  // * pitch sizes
  //
  using namespace geo::pixel;
  fPixelDirs = {
    getSensElemDir(ixMain) * (getSensElemPitch(ixMain) / 2.0),
    getSensElemDir(ixSec) * (getSensElemPitch(ixSec) / 2.0)
    };
} // geo::PixelPlaneGeo::UpdatePixelDirs()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdateDecompPixelOrigin() {

  //
  // update the origin of the reference frame (the middle of the first pixel)
  // requires:
  // * width and depth size
  // * plane center in width and depth direction (thickness is irrelevant)
  // * pixel plane width and depth directions (center is what we find here)
  // * pixel pitches
  //
  
  //
  // for this time only, we need to find the location of the first pixel
  // from scratch: let's pick it from the most negative corner in pixel
  // coordinates (half width and depth before the plane [box] center).
  // Let's pick it also half a thickness length toward the center of the TPC
  // (assuming the normal to the plane is pointing there, as it should already).
  // 
  
  // translate the center of the box half the pixels back
  // (gets to the outer border of the first pixel)
  // and then add back half pixel (gets to the center of that pixel)
  geo::Point_t center = GetBoxCenter<geo::Point_t>();
  for (DirIndex_t dir = 0; dir < geo::pixel::NCoords; ++dir) {
    center += getSensElemHalfStepDir(dir)
      * (1.0 - static_cast<double>(getNsensElem(dir)));
  }
  
  // thickness now... this is expensive:
  double const thickness = extractPlaneThickness();
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Pixel plane box is " << thickness << " cm thick";
  center += (thickness / 2.0) * fDecompPixel.NormalDir();
  
  // save it now
  fDecompPixel.SetOrigin(center);

} // geo::PixelPlaneGeo::UpdateDecompPixelOrigin()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdatePlaneCenter() {
  //
  // update the center of the plane and the center of the frame base
  //
  
  //
  // The center of the sensitive plane is defined as the center of the plane
  // box, translated to the plane the sensitive elements lie on.
  // This assumes that the thickness direction of the box is aligned with
  // the drift direction, so that the translated point is still in the middle
  // of width and depth dimensions.
  // It is possible to remove that assumption by translating the center of the
  // box along the thickness direction enough to bring it to the target plane.
  // The math is just a bit less straightforward, so we don't bother yet.
  //
  // Requirements:
  //  * the pixel decomposition frame must be set up (at least its origin and
  //    normal direction)
  //

  fCenter = GetBoxCenter<geo::Point_t>();

  DriftPoint(fCenter, fDecompPixel.PointNormalComponent(fCenter));
  
  geo::vect::round0(fCenter, 1e-7); // round dimensions less than 1 nm to 0
  
  fDecompFrame.SetOrigin(fCenter); // equivalent to GetCenter() now

} // geo::PixelPlaneGeo::UpdatePlaneCenter()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdateDirections() {
  
  //
  // computes the angles out of the pixel frame directions
  // requires: pixel frame directions
  //
  
  fThetaZ = std::acos(GetWireDirection<geo::Vector_t>().Dot(geo::Zaxis()));
  fPhiZ
    = std::acos(GetIncreasingWireDirection<geo::Vector_t>().Dot(geo::Zaxis()));
  
} // geo::PixelPlaneGeo::UpdateDirections()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdateActiveArea() {
  
  //
  // The active area is defined in the width/depth space which include
  // approximatively all pixels; if pixels are aligned with width and depth
  // (which they are), the active area is *exactly* the area of the pixels.
  // 
  // Requires: number of pixels and pitch in both directions
  //
  
  using namespace geo::pixel;
  double const halfWidth = getSensElemPitch(ixMain)
    * (static_cast<double>(getNsensElem(ixMain)) / 2.0);
  double const halfDepth = getSensElemPitch(ixSec)
    * (static_cast<double>(getNsensElem(ixSec)) / 2.0);
  fActiveArea = { { -halfWidth, halfWidth  }, { -halfDepth, halfDepth } };

} // geo::PixelPlaneGeo::UpdateActiveArea()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::firstPixelCenter() const {
  return fDecompPixel.ReferencePoint();
} // geo::PixelPlaneGeo::firstPixelCenter()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::discoverPitches() {
  // TODO believe it or not, this is not how we expect to discover pixel size
  fPitches = { 3.0, 4.0 }; // mm
} // geo::PixelPlaneGeo::discoverPitches()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeo::extractPlaneThickness() const {
  
  //
  // finding the plane thickness at this point is hard because we have not
  // stored its value in any form and we have forgotten which rotation turned
  // the original box, that we still know, into the final one;
  // so we are going a bit heuristic here:
  // 1) we can get (at a cost) the size of the three sides, but we don't know
  //    which one is which.
  // 2) we remove all the ones we know... the one left must be the thickness.
  //
  
  // behold the plane coordinate box in world coordinates.
  geo::BoxBoundedGeo const& box = BoundingBox();
  
  std::array<double, 3U> sizes = { box.SizeX(), box.SizeY(), box.SizeZ() };
  
  for (double size: { Width(), Depth() }) {
    for (double& boxSize: sizes) {
      if (std::abs(boxSize - size) > 1e-4) continue;
      boxSize = 0.0; // found it: remove it from the list
      break;
    } // for (inner)
  } // for (outer)
  
  // return the survivor
  for (double boxSize: sizes) if (boxSize != 0.0) return boxSize;
  return 0.0; // huh?
} // geo::PixelPlaneGeo::extractPlaneThickness()


// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeo::getDirectionName(DirIndex_t const dir) {
  using namespace std::string_literals;
  switch (dir) {
    case geo::pixel::ixMain: return "main"s;
    case geo::pixel::ixSec:  return "secondary"s;
    default:                 return "invalid"s;
  } // switch
} // geo::PixelPlaneGeo::getDirectionName()


// -----------------------------------------------------------------------------
