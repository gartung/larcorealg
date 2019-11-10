/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeoInterface.cxx
 * @brief   Representation of a readout plane with pixels as sensitive elements.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 23, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeoInterface.h`
 * @ingroup GeometryPixel
 */

// class header
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeoInterface.h"

// LArSoft includes
#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::rounded01()
#include "larcorealg/CoreUtils/RealComparisons.h" // makeVector3DComparison()
#include "larcorealg/CoreUtils/zip.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Zaxis()

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ standard library
#include <sstream> // std::ostringstream
#include <array>
#include <algorithm> // std::clamp(), std::swap(), ...
#include <utility> // std::move()
#include <limits> // std::numeric_limits<>
#include <cmath> // std::cos(), std::abs(), ...
#include <cassert>


// -----------------------------------------------------------------------------
// --- geo::PixelPlaneGeoInterface
// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeoInterface::PixelInfo(
  PixelCoordID_t const coords,
  std::string indent /* = "" */, unsigned int verbosity /* = 1U */
  ) const
{
  std::ostringstream sstr;
  PrintPixelInfo(sstr, coords, indent, verbosity);
  return sstr.str();
} // geo::PixelPlaneGeoInterface::PixelInfo()


// -----------------------------------------------------------------------------
// ---  interface implementation: anode plane
// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doThetaZ() const {
  return fThetaZ;
} // geo::PixelPlaneGeoInterface::doThetaZ()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doPhiZ() const {
  return fPhiZ;
} // geo::PixelPlaneGeoInterface::doPhiZ()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doSinPhiZ() const {
  // if these need to be cached... can do
  return std::sin(fPhiZ);
} // geo::PixelPlaneGeoInterface::doSinPhiZ()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doCosPhiZ() const {
  // if these need to be cached... can do
  return std::cos(fPhiZ);
} // geo::PixelPlaneGeoInterface::doCosPhiZ()


// -----------------------------------------------------------------------------
geo::WireGeo geo::PixelPlaneGeoInterface::doWire(unsigned int iWire) const {
  return getWire(iWire);
} // geo::PixelPlaneGeoInterface::doWire()


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeoInterface::doNwires() const {
  
  return getNsensElem();
  
} // geo::PixelPlaneGeoInterface::doNwires()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWirePitch() const
  { return getSensElemPitch(geo::pixel::ixWireC); }


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::doWireIDincreasesWithZ() const {
  assert(getSensElemDir(geo::pixel::ixWireC).Dot(geo::Zaxis()) > 0.0);
  return true; // pretty much by the definition of the wire ID's
} // geo::PixelPlaneGeoInterface::doWireIDincreasesWithZ()


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::doGetIncreasingWireDirection() const
  { return getSensElemDir(geo::pixel::ixWireC); }


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::doGetWireDirection() const
  { return getSensElemDir(geo::pixel::ixWireD); }


// -----------------------------------------------------------------------------
geo::WireID geo::PixelPlaneGeoInterface::doNearestWireID
  (geo::Point_t const& pos) const
{
  
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
  
} // geo::PixelPlaneGeoInterface::doNearestWireID()


// -----------------------------------------------------------------------------
geo::WireGeo geo::PixelPlaneGeoInterface::doNearestWire
  (geo::Point_t const& pos) const
{
  
  return getWire(geo::PixelPlaneGeoInterface::doNearestWireID(pos).Wire);
  
} // geo::PixelPlaneGeoInterface::doNearestWire()


// -----------------------------------------------------------------------------
geo::WireID geo::PixelPlaneGeoInterface::doClosestWireID
  (geo::WireID::WireID_t) const
{
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
  
} // geo::PixelPlaneGeoInterface::doClosestWireID()


// -----------------------------------------------------------------------------
lar::util::simple_geo::Volume<> geo::PixelPlaneGeoInterface::doCoverage() const
{
  
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
} // geo::PixelPlaneGeoInterface::WirePlaneGeo::doCoverage()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doPlaneCoordinateFrom
  (geo::Point_t const& point, geo::WireGeo const& refPixel) const
{
  return getPlaneCoordinateFrom(point, refPixel, geo::pixel::ixWireC);
} // geo::PixelPlaneGeoInterface::doPlaneCoordinateFrom()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doPlaneCoordinate
  (geo::Point_t const& point) const
{
  return getPlaneCoordinate(point, geo::pixel::ixWireC);
} // geo::PixelPlaneGeoInterface::doPlaneCoordinate()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireCoordinate
  (geo::Point_t const& point) const
{
  return getWireCoordinate(point, geo::pixel::ixWireC);
} // geo::PixelPlaneGeoInterface::doWireCoordinate()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::doDecomposePoint
  (geo::Point_t const& point) const
  -> WireDecomposedVector_t
{
  return fDecompPixel.DecomposePoint(point);
} // geo::PixelPlaneGeoInterface::doDecomposePoint()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doProjectionReferencePoint() const {
  return fDecompPixel.ReferencePoint();
} // geo::PixelPlaneGeoInterface::doProjectionReferencePoint()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::doProjection(geo::Point_t const& point) const
  -> WireCoordProjection_t
{
  return getPlaneCoordinates(point);
} // geo::PixelPlaneGeoInterface::doProjection()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::doProjection(geo::Vector_t const& v) const
  -> WireCoordProjection_t
{
  return fDecompPixel.ProjectVectorOnPlane(v);
} // geo::PixelPlaneGeoInterface::doProjection()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doComposePoint
  (WireDecomposedVector_t const& decomp) const
{
  return fDecompPixel.ComposePoint(decomp); 
} // geo::PixelPlaneGeoInterface::doComposePoint(WireDecomposedVector_t)


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doComposePoint
  (double distance, WireCoordProjection_t const& proj) const
{
  return fDecompPixel.ComposePoint(distance, proj); 
} // geo::PixelPlaneGeoInterface::doComposePoint(double, WireCoordProjection_t)


// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeoInterface::doPlaneInfo
  (std::string indent /* = "" */, unsigned int verbosity /* = 1U */) const
{
  std::ostringstream out;
  PrintPixelPlaneInfo(out, indent, verbosity);
  return out.str();
} // geo::PixelPlaneGeoInterface::doPlaneInfo()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::doUpdateAfterSorting
  (geo::BoxBoundedGeo const& /* TPCbox */)
{
  
  MF_LOG_DEBUG("PixelPlaneGeoInterface") << "doUpdateAfterSorting() for " << ID();
  
  UpdateDecompPixel();
  
  UpdatePixelDirs();
  
  UpdatePlaneCenter();
  
  UpdateAngles();
  
  UpdateActiveArea();
  
} // geo::PixelPlaneGeoInterface::doUpdateAfterSorting()


// -----------------------------------------------------------------------------
// --- Polymorphic implementation: wire abstraction
// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireRMax(WireLocator const&) const {
  using namespace geo::pixel;
  return std::max(getSensElemPitch(ixMain), getSensElemPitch(ixSec)) / 2.0;
} // geo::PixelPlaneGeoInterface::doWireRMax()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireRMin(WireLocator const&) const {
  using namespace geo::pixel;
  return std::min(getSensElemPitch(ixMain), getSensElemPitch(ixSec)) / 2.0;
} // geo::PixelPlaneGeoInterface::doWireRMin()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireHalfL(WireLocator const& wloc) const {
  return getPixelHalfL(wloc);
} // geo::PixelPlaneGeoInterface::doWireHalfL()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::doWireFillCenterXYZ
  (WireLocator const& wloc, double* xyz, double localz /* = 0.0 */) const
{
  geo::vect::fillCoords(xyz, doWireGetPositionFromCenter(wloc, localz));
} // geo::PixelPlaneGeoInterface::doWireFillCenterXYZ()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::doWireFillStartXYZ
  (WireLocator const& wloc, double* xyz) const
{
  geo::vect::fillCoords(xyz, doWireGetStart(wloc));
} // geo::PixelPlaneGeoInterface::doWireFillStartXYZ()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::doWireFillEndXYZ
  (WireLocator const& wloc, double* xyz) const
{
  geo::vect::fillCoords(xyz, doWireGetEnd(wloc));
} // geo::PixelPlaneGeoInterface::doWireFillEndXYZ()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doWireGetPositionFromCenter
  (WireLocator const& wloc, double localz) const
{
  return
    doWireGetPositionFromCenterUnbounded(wloc, std::clamp(localz, -1.0, +1.0));
} // geo::PixelPlaneGeoInterface::doWireGetPositionFromCenter()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doWireGetPositionFromCenterUnbounded
  (WireLocator const& wloc, double localz) const
{
  return getSensElemCenter(coordsOf(wloc))
    + localz * getSensElemHalfStepDir(geo::pixel::ixWireD);
} // geo::PixelPlaneGeoInterface::doWireGetPositionFromCenterUnbounded()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doWireGetCenter
  (WireLocator const& wloc) const
{
  return getSensElemCenter(coordsOf(wloc));
} // geo::PixelPlaneGeoInterface::doWireGetCenter()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doWireGetStart
  (WireLocator const& wloc) const
{
  using namespace geo::pixel;
  return getSensElemCenter(coordsOf(wloc))
    - getSensElemHalfStepDir(ixWireC)
    - getSensElemHalfStepDir(ixWireD)
    ;
} // geo::PixelPlaneGeoInterface::doWireGetStart()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::doWireGetEnd
  (WireLocator const& wloc) const
{
  using namespace geo::pixel;
  return getSensElemCenter(coordsOf(wloc))
    + getSensElemHalfStepDir(ixWireC)
    + getSensElemHalfStepDir(ixWireD)
    ;
} // geo::PixelPlaneGeoInterface::doWireGetEnd()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireLength(WireLocator const& wloc) const
{
  return 2.0 * getSensElemPitch(geo::pixel::ixWireD);
} // geo::PixelPlaneGeoInterface::doWireLength()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireThetaZ(WireLocator const&) const {
  return ThetaZ();
} // geo::PixelPlaneGeoInterface::doWireThetaZ()


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::doWireIsParallelTo
  (WireLocator const&, geo::WireGeo const& wire) const
{
  return lar::util::Vector3DComparison(1e-4).parallel
    (getSensElemDir(geo::pixel::ixWireD), wire.Direction<geo::Vector_t>());
} // geo::PixelPlaneGeoInterface::doWireIsParallelTo()


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::doWireDirection
  (WireLocator const&) const
{
  return getSensElemDir(geo::pixel::ixWireD);
} // geo::PixelPlaneGeoInterface::doWireDirection()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireDistanceFrom
  (WireLocator const& wloc, geo::WireGeo const& wire) const
{
  // can't use getSecPlaneCoordinate() because it uses pixel #0 reference
  return fDecompPixel.VectorSecondaryComponent
    (getSensElemCenter(coordsOf(wloc)) - wire.GetCenter<geo::Point_t>());
} // geo::PixelPlaneGeoInterface::doWireDistanceFrom()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doWireComputeZatY0
  (WireLocator const& wloc) const
{
  
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
  
} // geo::PixelPlaneGeoInterface::doWireComputeZatY0()


// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeoInterface::doWireInfo(
  WireLocator const& wloc,
  std::string indent /* = "" */, unsigned int verbosity /* = 1 */
  ) const
{
  return PixelInfo(coordsOf(wloc), indent, verbosity);
} // geo::PixelPlaneGeoInterface::doWireInfo()


// -----------------------------------------------------------------------------
// --- Plane coordinate and ID conversions
// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::indexAt
  (WireCoordProjection_t const& point) const
  -> PixelIndex_t
{
  return indexOf(coordsOf(point));
} // geo::PixelPlaneGeoInterface::indexAt()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::coordsOf
  (WireCoordProjection_t const& point) const
  -> PixelCoordID_t
{
  using namespace geo::pixel;
  return { roundCoord(point.X(), ixMain), roundCoord(point.Y(), ixSec) };
} // geo::PixelPlaneGeoInterface::coordsOf(WireCoordProjection_t)


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::coordsOf(WireLocator const& wloc) const
  -> PixelCoordID_t
{
  return coordsOf(wireToPixelIndex(wloc));
} // geo::PixelPlaneGeoInterface::coordsOf(WireLocator)


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::coordsOf(PixelIndex_t const index) const
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
} // geo::PixelPlaneGeoInterface::coordsOf()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::roundCoord
  (double coord, DirIndex_t const dir) const -> PixelCoordIndex_t
{
  return roundPixelCoord(pixelCoord(coord, dir));
} // geo::PixelPlaneGeoInterface::roundCoord()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::pixelCoord
  (double coord, DirIndex_t const dir) const
{
  assert(isDirIndex(dir));
  return coord / getSensElemPitch(dir);
} // geo::PixelPlaneGeoInterface::pixelCoord()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::roundPixelCoord(double const coord) const
  -> PixelCoordIndex_t
{
  return static_cast<PixelCoordIndex_t>(coord + 0.5);
} // geo::PixelPlaneGeoInterface::roundPixelCoord()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::indexOf(PixelCoordID_t const& coords) const
  -> PixelIndex_t
{
  using geo::pixel::ixMain;
  return coords.secondary() * getNsensElem(ixMain) + coords.main();
} // geo::PixelPlaneGeoInterface::indexOf()


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isOnPlane
  (WireCoordProjection_t const& point) const
{
  using namespace geo::pixel;
  return isOnPlane(point.X(), ixMain) && isOnPlane(point.Y(), ixSec);
} // geo::PixelPlaneGeoInterface::isOnPlane(WireCoordProjection_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isOnPlane
  (double const coord, DirIndex_t const dir) const
{
  // yeah, strictly speaking it should be `<` instead of `<=`. Judge me.
  double const shiftedCoord = coord + getSensElemPitch(dir) / 2.0;
  return (shiftedCoord >= 0.0) && (shiftedCoord <= getSensElemDirSize(dir));
} // geo::PixelPlaneGeoInterface::isOnPlane(double)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isPixelOnPlane
  (PixelCoords_t const& coords) const
{
  using namespace geo::pixel;
  return isPixelOnPlane(coords.main(), ixMain)
    && isPixelOnPlane(coords.secondary(), ixSec);
} // bool geo::PixelPlaneGeoInterface::isOnPlane(PixelCoords_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isPixelOnPlane
  (PixelCoordID_t const& coords) const
{
  using namespace geo::pixel;
  return isPixelOnPlane(coords.main(), ixMain)
    && isPixelOnPlane(coords.secondary(), ixSec);
} // geo::PixelPlaneGeoInterface::isOnPlane(PixelCoordID_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isPixelOnPlane
  (double const coord, DirIndex_t const dir) const
{
  double const shiftedCoord = 0.5 + coord;
  return (shiftedCoord >= 0.0)
    && (shiftedCoord < static_cast<double>(getNsensElem(dir)));
} // geo::PixelPlaneGeoInterface::isPixelOnPlane(double)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isPixelOnPlane
  (PixelCoordIndex_t const coord, DirIndex_t const dir) const
{
  return
    (coord >= 0) && (coord < static_cast<PixelCoordIndex_t>(getNsensElem(dir)));
} // geo::PixelPlaneGeoInterface::isPixelOnPlane(PixelCoordIndex_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isPixelOnPlane(PixelIndex_t const index) const
{
  return index < getNsensElem();
} // geo::PixelPlaneGeoInterface::isOnPlane(PixelIndex_t)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isWireIDvalid(geo::WireID const& wireid) const
{
  return isWireIDvalid(wireid.Wire);
} // geo::PixelPlaneGeoInterface::isWireIDvalid(WireID)


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoInterface::isWireIDvalid
  (geo::WireID::WireID_t const wire) const
{
  return (wire >= 0) && (wire < getNsensElem());
} // geo::PixelPlaneGeoInterface::isWireIDvalid(WireID_t)


// -----------------------------------------------------------------------------
// --- Candidate extensions to `geo::PlaneGeo` interface
// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::doSensElemDir
  (DirIndex_t const dir) const
  { return getSensElemDir(dir); }


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::doSensElemMainDir() const
  { return getSensElemDir(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::doSensElemSecondaryDir() const
  { return getSensElemDir(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doSensElemDirSize
  (DirIndex_t const dir) const
  { return getSensElemDirSize(dir); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doSensElemMainDirSize() const
  { return getSensElemDirSize(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doSensElemSecondaryDirSize() const
  { return getSensElemDirSize(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeoInterface::doNsensElem
  (DirIndex_t const dir) const
  { return getNsensElem(dir); }


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeoInterface::doNsensElemMain() const
  { return getNsensElem(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeoInterface::doNsensElemSecondary() const
  { return getNsensElem(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doSensElemPitch(DirIndex_t const dir) const
  { return getSensElemPitch(dir); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doSensElemPitchMain() const
  { return getSensElemPitch(geo::pixel::ixMain); }


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::doSensElemPitchSecondary() const
  { return getSensElemPitch(geo::pixel::ixSec); }


// -----------------------------------------------------------------------------
geo::WireGeo geo::PixelPlaneGeoInterface::getWire
  (PixelIndex_t const iWire) const
{
  return { *this, { ID(), static_cast<geo::WireID::WireID_t>(iWire) } };
} // geo::PixelPlaneGeoInterface::getWire()


// -----------------------------------------------------------------------------
// --- Implementation of the candidate interface extensions
// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::getSensElemDir
  (DirIndex_t const dir) const
{
  assert(isDirIndex(dir));
  switch (dir) {
    case geo::pixel::ixMain: return fDecompPixel.MainDir();
    case geo::pixel::ixSec:  return fDecompPixel.SecondaryDir();
    default:     return {};
  } // switch
} // geo::PixelPlaneGeoInterface::getSensElemDir()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getSensElemDirSize
  (DirIndex_t const dir) const
{
  //
  // this code assumes that the main direction matches the frame width direction
  // and the secondary direction msatches the frame depth direction
  //
  constexpr lar::util::Vector3DComparison vcmp(1e-4);
  assert(isDirIndex(dir));
  switch (dir) {
    case geo::pixel::ixMain:
      assert(vcmp.parallel(DepthDir<geo::Vector_t>(), getSensElemDir(dir)));
      return fFrameSize.Depth();
    case geo::pixel::ixSec:
      assert(vcmp.parallel(WidthDir<geo::Vector_t>(), getSensElemDir(dir)));
      return fFrameSize.Width();
    default:
      return 0;
  } // switch
} // geo::PixelPlaneGeoInterface::getSensElemDirSize()


// -----------------------------------------------------------------------------
geo::Vector_t geo::PixelPlaneGeoInterface::getSensElemHalfStepDir
  (DirIndex_t const dir) const
{
  assert(isDirIndex(dir));
  return fPixelDirs[dir];
} // geo::PixelPlaneGeoInterface::getSensElemHalfStepDir()


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeoInterface::getNsensElem
  (DirIndex_t const dir) const
{
  assert(isDirIndex(dir));
  return fNPixels[dir];
} // geo::PixelPlaneGeoInterface::getNsensElem(DirIndex_t)


// -----------------------------------------------------------------------------
unsigned int geo::PixelPlaneGeoInterface::getNsensElem() const {
  using namespace geo::pixel;
  return getNsensElem(ixMain) * getNsensElem(ixSec);
} // geo::PixelPlaneGeoInterface::getNsensElem()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getSensElemPitch
  (DirIndex_t const dir) const
{
  assert(isDirIndex(dir));
  return fPitches[dir];
} // geo::PixelPlaneGeoInterface::getSensElemPitch()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getPlaneCoordinateFrom
  (geo::Point_t const& point, geo::WireGeo const& ref, DirIndex_t const dir)
  const
{
  return getPlaneCoordinate(
    point - (ref.GetCenter<geo::Point_t>() - fDecompPixel.ReferencePoint()),
    dir
    );
} // geo::PixelPlaneGeoInterface::getPlaneCoordinateFrom()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoInterface::getPlaneCoordinates
  (geo::Point_t const& point) const -> WireCoordProjection_t
{
  return fDecompPixel.ProjectPointOnPlane(point);
} // geo::PixelPlaneGeoInterface::getPlaneCoordinates()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getPlaneCoordinate
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
  
} // geo::PixelPlaneGeoInterface::getPlaneCoordinate()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getMainPlaneCoordinate
  (geo::Point_t const& point) const
{
  return fDecompPixel.PointMainComponent(point);
} // geo::PixelPlaneGeoInterface::getMainPlaneCoordinate()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getSecPlaneCoordinate
  (geo::Point_t const& point) const
{
  return fDecompPixel.PointSecondaryComponent(point);
} // geo::PixelPlaneGeoInterface::getSecPlaneCoordinate()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getWireCoordinate
  (geo::Point_t const& point, DirIndex_t const dir) const
{
  // no rounding, no shifting; pixel #0 is covering coordinates from -0.5 to 0.5
  // in pitch units
  return getPlaneCoordinate(point, dir) / getSensElemPitch(dir);
} // geo::PixelPlaneGeoInterface::getWireCoordinate()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::getSensElemCenter
  (PixelCoordID_t const& coords) const
{
  using namespace geo::pixel;
  return firstPixelCenter()
    + getSensElemHalfStepDir(ixMain) * (2.0 * coords[ixMain])
    + getSensElemHalfStepDir(ixSec) * (2.0 * coords[ixSec])
    ;
} // geo::PixelPlaneGeoInterface::getSensElemCenter()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::getPixelHalfL(WireLocator const&) const {
  return getSensElemPitch(geo::pixel::ixWireD) / 2.0;
} // geo::PixelPlaneGeoInterface::getPixelHalfL()


// -----------------------------------------------------------------------------
// --- initialization customization for derived classes
// -----------------------------------------------------------------------------
geo::PixelPlaneGeoInterface::PixelPlaneGeoInterface(
  TGeoNode const& node,
  geo::TransformationMatrix&& trans
  )
  : geo::PlaneGeo(node, std::move(trans))
  , fDecompPixel()
  , fThetaZ(0.0)
  , fPhiZ(0.0)
{

  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PixelPlaneGeoInterface` class. If changing *what* is being
   *      initialized or its *order*, please also update that documentation at
   *      the top of `geo::PixelPlaneGeoInterface` class Doxygen documentation
   *      (in `larcorealg/Geometry/PixelPlane/PixelPlaneGeoInterface.h`, section
   *      "Initialization steps").
   */
  
  /*
   * REMINDER: do not rely on virtual methods from derived classes here, as
   *           they might not be available yet (the derived class constructor
   *           hasn't been run yet at this point)
   */
  
  fNPixels.fill(0);
  fPitches.fill(0.0);
  fPixelDirs.fill({});
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Plane extends " << Width() << " cm in " << WidthDir<geo::Vector_t>()
    << " and " << Depth() << " cm in " << DepthDir<geo::Vector_t>();
  
  SetView(geo::k3D); // view is this simple
  
} // geo::PixelPlaneGeoInterface::PixelPlaneGeoInterface()


// -----------------------------------------------------------------------------
// --- implementation details (private methods)
// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::UpdateDecompPixel() {
  
  constexpr auto vcmp = lar::util::makeVector3DComparison(1e-4); // 1 um tol.
  
  //
  // requires:
  //  * width and depth directions (including their verse)
  //  * old pixel #0 origin, and old pixel frame directions and pitches
  //
  assert(vcmp.nonZero(WidthDir<geo::Vector_t>()));
  assert(vcmp.nonZero(DepthDir<geo::Vector_t>()));
  assert(vcmp.nonZero(fDecompPixel.MainDir()));
  assert(vcmp.nonZero(fDecompPixel.SecondaryDir()));
  
  // we recover the center of the plane, which is going to be the same
  // after the redefinition of the axes; this is undoing a previous step
  // from the center of the plane, which we don't have available here,
  // to the first pixel which is the origin of the 
  geo::Point_t const pixelPlaneCenter
    = fromFirstPixelToCenter(fDecompPixel.ReferencePoint());
  
  //
  // direction measured by the secondary coordinate; it is the width direction
  //
  fDecompPixel.SetSecondaryDir
    (geo::vect::rounded01(WidthDir<geo::Vector_t>(), 1e-4));
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Pixel frame secondary dir set to: " << fDecompPixel.SecondaryDir()
    << " (width: " << WidthDir<geo::Vector_t>() << ")";
  
  //
  // get the axis perpendicular to it on the wire plane
  // (verse is already set to have the base as positive as the frame base is)
  //
  fDecompPixel.SetMainDir
    (geo::vect::rounded01(-DepthDir<geo::Vector_t>(), 1e-4));
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Pixel frame main dir set to: " << fDecompPixel.MainDir()
    << " (depth: " << DepthDir<geo::Vector_t>() << ")";
  //
  // check that the resulting normal matches the plane one
  //
  assert
    (vcmp.equal(fDecompPixel.NormalDir(), GetNormalDirection<geo::Vector_t>()));
  
  //
  // set the new center; it will may lie on a different corner of the plane
  //
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Pixel area center moved: " << fDecompPixel.ReferencePoint()
    << " => " << pixelPlaneCenter << "...";
  
  fDecompPixel.SetReferencePoint(fromCenterToFirstPixel(pixelPlaneCenter));
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Pixel area center moved: ... " << pixelPlaneCenter
    << " => " << fDecompPixel.ReferencePoint();
  
} // geo::PixelPlaneGeoInterface::UpdateDecompPixel()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::UpdatePixelDirs() {
  //
  // requirements:
  // * "wire" frame direction definitions completely finalized
  // * pitch sizes
  //
  using namespace geo::pixel;
  
  constexpr auto vcmp = lar::util::makeVector3DComparison(1e-4); // 1 um tol.
  
  assert(vcmp.nonZero(getSensElemDir(ixMain)));
  
  fPixelDirs = {
    getSensElemDir(ixMain) * (getSensElemPitch(ixMain) / 2.0),
    getSensElemDir(ixSec) * (getSensElemPitch(ixSec) / 2.0)
    };
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Pixel half steps set to: " << fPixelDirs[ixMain] << " (main), "
    << fPixelDirs[ixSec] << " (secondary)";
  
} // geo::PixelPlaneGeoInterface::UpdatePixelDirs()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::UpdatePlaneCenter() {
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
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Pixel frame origin: " << fDecompPixel.ReferencePoint()
    << " => " << fCenter
    << " (box: " << GetBoxCenter<geo::Point_t>() << ")";
  
  fDecompFrame.SetReferencePoint(fCenter); // equivalent to GetCenter() now

} // geo::PixelPlaneGeoInterface::UpdatePlaneCenter()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::UpdateAngles() {
  
  //
  // computes the angles out of the pixel frame directions
  // requires: pixel frame directions
  //
  
  fThetaZ = std::acos(GetWireDirection<geo::Vector_t>().Dot(geo::Zaxis()));
  fPhiZ
    = std::acos(GetIncreasingWireDirection<geo::Vector_t>().Dot(geo::Zaxis()));
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Angles set to: thetaZ=" << fThetaZ << " rad, phiZ=" << fPhiZ << " rad";
  
} // geo::PixelPlaneGeoInterface::UpdateAngles()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoInterface::UpdateActiveArea() {
  
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

  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "Active area set to:"
    << "\n - width: "
      << fActiveArea.width.lower << " -- " << fActiveArea.width.upper
    << "\n - depth: "
      << fActiveArea.depth.lower << " -- " << fActiveArea.depth.upper
    ;
  
} // geo::PixelPlaneGeoInterface::UpdateActiveArea()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::firstPixelCenter() const {
  return fDecompPixel.ReferencePoint();
} // geo::PixelPlaneGeoInterface::firstPixelCenter()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::fromCenterToFirstPixel
  (geo::Point_t const& pixelPlaneCenter) const
{
  
  /*
   * Converts the center of the sensitive area of the plane (`pixelPlaneCenter`)
   * into the center of the first pixel.
   * Requires:
   *  * number of pixels
   *  * directions of the two axes of the pixel plane
   *  * pitch of the pixels
   * 
   */
  using namespace geo::pixel;
  
  constexpr auto vcmp = lar::util::makeVector3DComparison(1e-4);
  
  assert(getNsensElem(ixMain) > 0);
  assert(getNsensElem(ixSec) > 0);
  assert(getSensElemPitch(ixMain) > 0.0);
  assert(getSensElemPitch(ixSec) > 0.0);
  assert(vcmp.nonZero(getSensElemDir(ixMain)));
  assert(vcmp.nonZero(getSensElemDir(ixSec)));
  
  geo::Point_t center = pixelPlaneCenter;
  
  for (DirIndex_t dir = 0; dir < geo::pixel::NCoords; ++dir) {
    center += getSensElemDir(dir) * (
      (1.0 - static_cast<double>(getNsensElem(dir)))
      * getSensElemPitch(dir) / 2.0
      );
  } // for
  
  return center;
} // geo::PixelPlaneGeoInterface::fromCenterToFirstPixel()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeoInterface::fromFirstPixelToCenter
  (geo::Point_t const& pixelCenter) const
{
  
  /*
   * Converts the center of the sensitive area of the plane (`pixelPlaneCenter`)
   * into the center of the first pixel.
   * Requires:
   *  * number of pixels
   *  * directions of the two axes of the pixel plane
   *  * pitch of the pixels
   * 
   */
  using namespace geo::pixel;
  
  constexpr auto vcmp = lar::util::makeVector3DComparison(1e-4);
  
  assert(getNsensElem(ixMain) > 0);
  assert(getNsensElem(ixSec) > 0);
  assert(getSensElemPitch(ixMain) > 0.0);
  assert(getSensElemPitch(ixSec) > 0.0);
  assert(vcmp.nonZero(getSensElemDir(ixMain)));
  assert(vcmp.nonZero(getSensElemDir(ixSec)));
  
  geo::Point_t center = pixelCenter;
  
  for (DirIndex_t dir = 0; dir < geo::pixel::NCoords; ++dir) {
    center -= getSensElemDir(dir) * (
      (1.0 - static_cast<double>(getNsensElem(dir)))
      * getSensElemPitch(dir) / 2.0
      );
  } // for
  
  return center;
} // geo::PixelPlaneGeoInterface::fromFirstPixelToCenter()


// -----------------------------------------------------------------------------
double geo::PixelPlaneGeoInterface::extractPlaneThickness() const {
  
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
} // geo::PixelPlaneGeoInterface::extractPlaneThickness()


// -----------------------------------------------------------------------------
std::string geo::PixelPlaneGeoInterface::getDirectionName
  (DirIndex_t const dir)
{
  using namespace std::string_literals;
  switch (dir) {
    case geo::pixel::ixMain: return "main"s;
    case geo::pixel::ixSec:  return "secondary"s;
    default:                 return "invalid"s;
  } // switch
} // geo::PixelPlaneGeoInterface::getDirectionName()


// -----------------------------------------------------------------------------
