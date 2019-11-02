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
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Zaxis()

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

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
  
  /*
  // ---------------------------------------------------------------------------
  template <typename T>
  bool closeToUnity(T const value, T const tol = 1e-5) {
    return std::abs(value - T{ 1 }) <= tol;
  } // closeToUnity()
  */
  
  // ---------------------------------------------------------------------------
  /// Returns whether the vector `v` is not null.
  template <typename V, typename T = double>
  bool isNull(V const& v, T const tol = 1e-5) { return v.Mag2() <= (tol*tol); }
  
  
  // ---------------------------------------------------------------------------
  /// Returns whether `a` and `b` are orthogonal.
  template <typename VA, typename VB, typename T = double>
  bool areOrthogonal(VA const& a, VB const& b, T const tol = 1e-5) {
    return std::abs(a.Dot(b)) <= tol;
  } // areOrthogonal()
  
  
  // ---------------------------------------------------------------------------
  /// Returns whether `a` and `b` have the same direction (may be opposite).
  template <typename VA, typename VB, typename T = double>
  bool areParallel(VA const& a, VB const& b, T const tol = 1e-5) {
    return isNull(a.Cross(b), tol);
  } // areParallel()
  
  
  // ---------------------------------------------------------------------------
  /// Returns whether `b` projection on `a` is negative
  /// (@f$ \vec{a}\cdot\vec{b} < 0 @f$).
  template <typename VA, typename VB, typename T = double>
  bool areOpposite(VA const& a, VB const& b, T const tol = 1e-5) {
    return a.Dot(b) < -tol;
  } // areOpposite()
  
  
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
  , fDecompPixel()
  , fThetaZ(0.0)
  , fPhiZ(0.0)
{

  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PixelPlaneGeo` class. If changing *what* is being initialized
   *      or its *order*, please also update that documentation at the top
   *      of `geo::PixelPlaneGeo` class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`, section
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
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Plane extends " << Width() << " cm in " << WidthDir<geo::Vector_t>()
    << " and " << Depth() << " cm in " << DepthDir<geo::Vector_t>();
  
  SetView(geo::k3D); // view is this simple
  
  // TODO some updates should be moved here
  SquarePixelGeometry_t pixelGeometry = completePixelGeometry({});
  
  initializePixelGeometry(pixelGeometry);
  
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
  
  UpdatePixelDirs();
  
  UpdatePlaneCenter();
  
  UpdateAngles();
  
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
  return areParallel
    (getSensElemDir(geo::pixel::ixWireD), wire.Direction<geo::Vector_t>());
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
      assert(areParallel(DepthDir<geo::Vector_t>(), getSensElemDir(dir)));
      return fFrameSize.Depth();
    case geo::pixel::ixSec:
      assert(areParallel(WidthDir<geo::Vector_t>(), getSensElemDir(dir)));
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
auto geo::PixelPlaneGeo::completePixelGeometry
  (SquarePixelGeometry_t const& /* info */) const -> SquarePixelGeometry_t
{
  // TODO dummy stuff, and dumb too; assuming _x_ is the thin direction
  
  SquarePixelGeometry_t newInfo {
    {
      SquarePixelGeometry_t::AxisInfo_t
        { geo::Yaxis<LocalVector_t>() * Depth(), 3.0 }, // 3 cm pitch along local y axis
      SquarePixelGeometry_t::AxisInfo_t
        { geo::Zaxis<LocalVector_t>() * Width(), 4.0 }  // 4 cm pitch along local z axis
    },
    geo::origin<LocalPoint_t>()
  };
  
  return newInfo;
} // geo::PixelPlaneGeo::completePixelGeometry()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::initializePixelGeometry
  (SquarePixelGeometry_t const& pixelGeometry)
{
  
  /*
   * Produces the origin, pitch, number and direction of the pixels and grid.
   * This is a elaboration of the raw information from the geometry.
   * 
   * The origin and directions are only partially respected.
   * The directions are used to associate the pitches to the right axis, but
   * the axis are imposed to be the ones of the plane frame.
   * The origin vector is supposed to be one of the four corners of the plane,
   * the one from which the two input axes depart.
   * But as these axes may be swapped, the origin also can be modified
   * accordingly; it is anyway guaranteed to be at one of the four corners
   * of the pixelated area.
   * 
   * Requires:
   *  * local-world transformations be available
   *  * frame geometry to be set (exact origin depth coordinate does not matter)
   * 
   */
  
  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PixelPlaneGeo` class. If changing *what* is being initialized
   *      or its *order*, please also update that documentation at the top
   *      of `geo::PixelPlaneGeo` class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`, section
   *      "Initialization steps").
   */
  
  /*
   * REMINDER: do not rely on virtual methods from derived classes here, as
   *           they might not be available yet (this method is called in class
   *           constructor, when thee derived class constructor hasn't been run
   *           yet)
   */
  
  assert(!isNull(fDecompFrame.MainDir()));
  assert(!isNull(fDecompFrame.SecondaryDir()));
  
  using namespace geo::pixel;
  
  MF_LOG_TRACE("PixelPlaneGeo") << "Information received from plane geometry:"
    << "\n * center of the pixelized area: " << pixelGeometry.center
    << "\n * sides:"
    << "\n    - direction " << pixelGeometry.sides[0U].dir.Unit()
      << ": " << pixelGeometry.sides[0U].side << " cm pixels for "
      << pixelGeometry.sides[0U].dir.R() << " cm of side"
    << "\n    - direction " << pixelGeometry.sides[1U].dir.Unit()
      << ": " << pixelGeometry.sides[1U].side << " cm pixels for "
      << pixelGeometry.sides[1U].dir.R() << " cm of side"
    ;
  
  //
  // 1. assign the correct axis information to each direction;
  //    no nonsense here: axes must follow the frame one way or the other way
  //    (just two ways are allowed: either axis aligned with `WidthDir()`)
  //
  struct ExtendedAxisInfo_t {
    SquarePixelGeometry_t::AxisInfo_t basicInfo;
    geo::Vector_t axisDir; ///< Direction and size in the world frame.
  }; // struct ExtendedAxisInfo_t
  
  static_assert(pixelGeometry.sides.size() == NCoords);
  
  std::array<ExtendedAxisInfo_t, NCoords> axes;
  for (auto&& [ side, axis ]: util::zip(pixelGeometry.sides, axes))
    axis = { side, toWorldCoords(side.dir) };
  
  //
  // make sure that the data on ixSec is aligned with WidthDir()`:
  //
  if (areParallel(axes[ixSec].axisDir, WidthDir())) {
    // all good already
  }
  else if (areParallel(axes[ixSec].axisDir, DepthDir())) {
    std::swap(axes[ixSec], axes[ixMain]);
  }
  else {
    throw cet::exception("PixelPlaneGeo")
      << "initializePixelGeometry(): pixel axis system { "
      << axes[ixSec].axisDir << " x " << axes[ixMain].axisDir
      << " } is misaligned with the plane frame system { "
      << WidthDir<geo::Vector_t>() << " x " << DepthDir<geo::Vector_t>() << " }"
      "\n";
  }
  
  //
  // copy the information into `fPitches` and `fNPixels`, and set the directions
  // (indices follow `ixMain` and `ixSec`)
  //
  for (auto&& [ axisInfo, pitch, nPixels ]: util::zip(axes, fPitches, fNPixels )) {
    pitch = axisInfo.basicInfo.side;
    nPixels = static_cast<unsigned int>(axisInfo.axisDir.R() / pitch);
  } // for
  
  fDecompPixel.SetMainDir
    (geo::vect::rounded01(axes[ixMain].axisDir.Unit(), 1e-5));
  fDecompPixel.SetSecondaryDir
    (geo::vect::rounded01(axes[ixSec].axisDir.Unit(), 1e-5));
  
  //
  // 2. now set the position of the pixel plane;
  //    that is driven by the origin of the pixel decomposition frame,
  //    which corresponds to the position of the first pixel:
  //    that is what we are pursuing now.
  //    We do not bother with the position along thickness direction,
  //    which stays the same as from the geometry
  //    (which is not that bad a choice);
  //    note that `fromCenterToFirstPixel()` requires some of the quantities
  //    just set (pretty much all of them, in fact)
  // 
  fDecompPixel.SetReferencePoint
    (fromCenterToFirstPixel(toWorldCoords(pixelGeometry.center)));
  
  MF_LOG_TRACE("PixelPlaneGeo") << "Directions set to:"
    << "\n * main: " << fDecompPixel.MainDir() << " after geometry dir "
      << axes[ixMain].axisDir << " with " << fNPixels[ixMain]
      << " pixels with " << fPitches[ixMain] << " cm pitch"
    << "\n * secondary: " << fDecompPixel.SecondaryDir() << " after geometry dir "
      << axes[ixSec].axisDir << " with " << fNPixels[ixSec]
      << " pixels with " << fPitches[ixSec] << " cm pitch"
    << "\n * origin: " << fDecompPixel.ReferencePoint()
      << " (center of pixelized area: " << toWorldCoords(pixelGeometry.center)
      << ")"
    ;
  
  
} // geo::PixelPlaneGeo::initializePixelGeometry()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdateDecompPixel() {
  
  //
  // requires:
  //  * width and depth directions (including their verse)
  //  * old pixel #0 origin, and old pixel frame directions and pitches
  //
  assert(!isNull(WidthDir<geo::Vector_t>()));
  assert(!isNull(DepthDir<geo::Vector_t>()));
  assert(!isNull(fDecompPixel.MainDir()));
  assert(!isNull(fDecompPixel.SecondaryDir()));
  
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
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Pixel frame secondary dir set to: " << fDecompPixel.SecondaryDir()
    << " (width: " << WidthDir<geo::Vector_t>() << ")";
  
  //
  // get the axis perpendicular to it on the wire plane
  // (verse is already set to have the base as positive as the frame base is)
  //
  fDecompPixel.SetMainDir
    (geo::vect::rounded01(-DepthDir<geo::Vector_t>(), 1e-4));
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Pixel frame main dir set to: " << fDecompPixel.MainDir()
    << " (depth: " << DepthDir<geo::Vector_t>() << ")";
  //
  // check that the resulting normal matches the plane one
  //
  assert(lar::util::makeVector3DComparison(1e-5)
    .equal(fDecompPixel.NormalDir(), GetNormalDirection<geo::Vector_t>()));
  
  //
  // set the new center; it will may lie on a different corner of the plane
  //
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Pixel area center moved: " << fDecompPixel.ReferencePoint()
    << " => " << pixelPlaneCenter << "...";
  
  fDecompPixel.SetReferencePoint(fromCenterToFirstPixel(pixelPlaneCenter));
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Pixel area center moved: ... " << pixelPlaneCenter
    << " => " << fDecompPixel.ReferencePoint();
  
} // geo::PixelPlaneGeo::UpdateDecompPixel()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdatePixelDirs() {
  //
  // requirements:
  // * "wire" frame direction definitions completely finalized
  // * pitch sizes
  //
  using namespace geo::pixel;
  
  assert(!isNull(getSensElemDir(ixMain)));
  
  fPixelDirs = {
    getSensElemDir(ixMain) * (getSensElemPitch(ixMain) / 2.0),
    getSensElemDir(ixSec) * (getSensElemPitch(ixSec) / 2.0)
    };
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Pixel steps set to: " << fPixelDirs[ixMain] << " (main), "
    << fPixelDirs[ixSec] << " (secondary)";
  
} // geo::PixelPlaneGeo::UpdatePixelDirs()


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
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Pixel frame origin: " << fDecompPixel.ReferencePoint()
    << " => " << fCenter
    << " (box: " << GetBoxCenter<geo::Point_t>() << ")";
  
  fDecompFrame.SetReferencePoint(fCenter); // equivalent to GetCenter() now

} // geo::PixelPlaneGeo::UpdatePlaneCenter()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeo::UpdateAngles() {
  
  //
  // computes the angles out of the pixel frame directions
  // requires: pixel frame directions
  //
  
  fThetaZ = std::acos(GetWireDirection<geo::Vector_t>().Dot(geo::Zaxis()));
  fPhiZ
    = std::acos(GetIncreasingWireDirection<geo::Vector_t>().Dot(geo::Zaxis()));
  
  MF_LOG_TRACE("PixelPlaneGeo")
    << "Angles set to: thetaZ=" << fThetaZ << " rad, phiZ=" << fPhiZ << " rad";
  
} // geo::PixelPlaneGeo::UpdateAngles()


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

  MF_LOG_TRACE("PixelPlaneGeo")
    << "Active area set to:"
    << "\n - width: " << fActiveArea.width.lower << " -- " << fActiveArea.width.upper
    << "\n - depth: " << fActiveArea.depth.lower << " -- " << fActiveArea.depth.upper
    ;
  
} // geo::PixelPlaneGeo::UpdateActiveArea()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::firstPixelCenter() const {
  return fDecompPixel.ReferencePoint();
} // geo::PixelPlaneGeo::firstPixelCenter()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::fromCenterToFirstPixel
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
  
  assert(getNsensElem(ixMain) > 0);
  assert(getNsensElem(ixSec) > 0);
  assert(getSensElemPitch(ixMain) > 0.0);
  assert(getSensElemPitch(ixSec) > 0.0);
  assert(!isNull(getSensElemDir(ixMain)));
  assert(!isNull(getSensElemDir(ixSec)));
  
  geo::Point_t center = pixelPlaneCenter;
  
  for (DirIndex_t dir = 0; dir < geo::pixel::NCoords; ++dir) {
    center += getSensElemDir(dir) * (
      (1.0 - static_cast<double>(getNsensElem(dir)))
      * getSensElemPitch(dir) / 2.0
      );
  } // for
  
  return center;
} // geo::PixelPlaneGeo::fromCenterToFirstPixel()


// -----------------------------------------------------------------------------
geo::Point_t geo::PixelPlaneGeo::fromFirstPixelToCenter
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
  
  assert(getNsensElem(ixMain) > 0);
  assert(getNsensElem(ixSec) > 0);
  assert(getSensElemPitch(ixMain) > 0.0);
  assert(getSensElemPitch(ixSec) > 0.0);
  assert(!isNull(getSensElemDir(ixMain)));
  assert(!isNull(getSensElemDir(ixSec)));
  
  geo::Point_t center = pixelCenter;
  
  for (DirIndex_t dir = 0; dir < geo::pixel::NCoords; ++dir) {
    center -= getSensElemDir(dir) * (
      (1.0 - static_cast<double>(getNsensElem(dir)))
      * getSensElemPitch(dir) / 2.0
      );
  } // for
  
  return center;
} // geo::PixelPlaneGeo::fromFirstPixelToCenter()


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
