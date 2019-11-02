/**
 * @file   larcorealg/Geometry/PlaneGeo.cxx
 * @brief  Interface for classes describing the geometry of a TPC readout plane.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   April 17, 2109
 * @ingroup Geometry
 * 
 */

// class header
#include "larcorealg/Geometry/PlaneGeo.h"

// LArSoft includes
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::round01()
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug namespace

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// ROOT includes
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TClass.h" // TClass::IsA()

// Boost libraries
#include "boost/iterator/transform_iterator.hpp"

// C/C++ standard library
#include <sstream> // std::ostringstream
#include <array>
#include <utility> // std::move()
#include <cmath> // std::abs()
#include <cassert>


namespace {

  /// Returns the offset to apply to value to move it inside [ -limit, +limit ].
  template <typename T>
  T symmetricCapDelta(T value, T limit) {

    return (value < -limit)
      ? -limit - value
      : (value > +limit)
        ? +limit - value
        : 0.0
      ;

  } // symmetricCapDelta()


  /// Returns a value shifted to fall into [ -limit; +limit ] interval.
  template <typename T>
  T symmetricCap(T value, T limit) {

    return value + symmetricCapDelta(value, limit);

  } // symmetricCap()

} // local namespace



//------------------------------------------------------------------------------
//---  geo::PlaneGeo
//------------------------------------------------------------------------------
geo::PlaneGeo::PlaneGeo(
  TGeoNode const& node,
  geo::TransformationMatrix&& trans
  )
  : fTrans(std::move(trans))
  , fVolume(node.GetVolume())
  , fView(geo::kUnknown)
  , fOrientation(geo::kVertical)
  , fCenter{}
  , fDecompFrame()
{

  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PlaneGeo` class. If changing *what* is being initialized,
   *      please also update that documentation at the top of `geo::PlaneGeo`
   *      class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`, section
   *      "Initialization details").
   */
  
  /*
   * REMINDER: do not rely on virtual methods from derived classes here, as
   *           they might not be available yet (the derived class constructor
   *           hasn't been run yet at this point)
   */
  
  if (!fVolume) {
    throw cet::exception("PlaneGeo")
      << "Plane geometry node " << node.IsA()->GetName()
      << "[" << node.GetName() << ", #" << node.GetNumber()
      << "] has no volume!\n";
  }

  DetectGeometryDirections();

} // geo::PlaneGeo::PlaneGeo()


//------------------------------------------------------------------------------
geo::BoxBoundedGeo geo::PlaneGeo::BoundingBox() const {

  //
  // The algorithm is not very refined...
  //

  TGeoBBox const* pShape = dynamic_cast<TGeoBBox const*>(fVolume->GetShape());
  if (!pShape) {
    throw cet::exception("PlaneGeo")
      << "BoundingBox(): volume " << fVolume->IsA()->GetName()
      << "['" << fVolume->GetName() << "'] has a shape which is a "
      << pShape->IsA()->GetName()
      << ", not a TGeoBBox!";
  }

  geo::BoxBoundedGeo box;
  unsigned int points = 0;
  for (double dx: { -(pShape->GetDX()), +(pShape->GetDX()) }) {
    for (double dy: { -(pShape->GetDY()), +(pShape->GetDY()) }) {
      for (double dz: { -(pShape->GetDZ()), +(pShape->GetDZ()) }) {

        auto const p = toWorldCoords(LocalPoint_t{ dx, dy, dz });

        if (points++ == 0)
          box.SetBoundaries(p, p);
        else
          box.ExtendToInclude(p);

      } // for z
    } // for y
  } // for x
  return box;

} // PlaneGeo::BoundingBox()


//------------------------------------------------------------------------------
std::string geo::PlaneGeo::GenericPlaneInfo(
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1 */
) const {
  
  std::ostringstream out;
  
  //----------------------------------------------------------------------------
  out << "plane " << std::string(ID());

  if (verbosity-- <= 0) return out.str(); // 0

  //----------------------------------------------------------------------------
  out
    << " at " << GetCenter<geo::Vector_t>() << " cm"
    << ", theta: " << ThetaZ() << " rad";

  if (verbosity-- <= 0) return out.str(); // 1

  //----------------------------------------------------------------------------
  unsigned int const nWires = Nwires();

  out << "\n" << indent
    << "normal to wire: " << PhiZ() << " rad"
      << ", with orientation " << OrientationName(Orientation())
      << ", has " << nWires << " wires measuring " << ViewName(View())
      << " with a wire pitch of " << WirePitch() << " cm"
    ;

  if (verbosity-- <= 0) return out.str(); // 2

  //----------------------------------------------------------------------------
  auto const& normal = GetNormalDirection<geo::Vector_t>();
  auto const& incrZdir = GetIncreasingWireDirection<geo::Vector_t>();
  out << "\n" << indent
    << "normal to plane: " << normal
    << ", direction of increasing wire number: " << incrZdir
    << " (" << (WireIDincreasesWithZ()? "increases": "decreases") << " with z)";

  if (verbosity-- <= 0) return out.str(); // 3

  //----------------------------------------------------------------------------

  auto const& wireDir = GetWireDirection<geo::Vector_t>();
  auto const& widthDir = WidthDir<geo::Vector_t>();
  auto const& depthDir = DepthDir<geo::Vector_t>();
  auto const& frameNormalDir = fDecompFrame.NormalDir();

  out << "\n" << indent
    << "wire direction: " << wireDir
    << "; width " << Width() << " cm in direction: " << widthDir
    << ", depth " << Depth() << " cm in direction: " << depthDir
    << " [normal: " << frameNormalDir << "]"
    ;

  if (verbosity-- <= 0) return out.str(); // 4

  //----------------------------------------------------------------------------
  // get the area spanned by the wires
  out << "\n" << indent << "wires cover width "
    << ActiveArea().width.lower << " to " << ActiveArea().width.upper
    << ", depth "
    << ActiveArea().depth.lower << " to " << ActiveArea().depth.upper
    << " cm";
  if (verbosity-- <= 0) return out.str(); // 5

  //----------------------------------------------------------------------------
  // print also the containing box
  auto const box = BoundingBox();
  out << "\n" << indent
    << "bounding box: " << box.Min() << " -- " << box.Max();

//  if (verbosity-- <= 0) return; // 6

  //----------------------------------------------------------------------------
  
  return out.str();
} // geo::PlaneGeo::GenericPlaneInfo()


//------------------------------------------------------------------------------
bool geo::PlaneGeo::isProjectionOnPlane(geo::Point_t const& point) const {

  auto const deltaProj
    = DeltaFromPlane(PointWidthDepthProjection(point));

  return (deltaProj.X() == 0.) && (deltaProj.Y() == 0.);

} // geo::PlaneGeo::isProjectionOnPlane()


//------------------------------------------------------------------------------
auto geo::PlaneGeo::DeltaFromPlane
  (WidthDepthProjection_t const& proj, double wMargin, double dMargin) const
  -> WidthDepthProjection_t 
{

  return {
    symmetricCapDelta(proj.X(), fFrameSize.HalfWidth() - wMargin),
    symmetricCapDelta(proj.Y(), fFrameSize.HalfDepth() - dMargin)
  };

} // geo::PlaneGeo::DeltaFromPlane()


//------------------------------------------------------------------------------
auto geo::PlaneGeo::DeltaFromActivePlane
  (WidthDepthProjection_t const& proj, double wMargin, double dMargin) const
  -> WidthDepthProjection_t
{

  return {
    fActiveArea.width.delta(proj.X(), wMargin),
    fActiveArea.depth.delta(proj.Y(), dMargin)
    };

} // geo::PlaneGeo::DeltaFromActivePlane()


//------------------------------------------------------------------------------
auto geo::PlaneGeo::MoveProjectionToPlane
  (WidthDepthProjection_t const& proj) const
  -> WidthDepthProjection_t 
{
  //
  // We have a more complicate implementation to avoid rounding errors.
  // In this way, the result is really guaranteed to be exactly on the border.
  //
  auto const delta = DeltaFromPlane(proj);
  return {
    (delta.X() == 0.0)
      ? proj.X()
      : ((delta.X() > 0)
        ? -fFrameSize.HalfWidth() // delta positive -> proj on negative side
        : fFrameSize.HalfWidth()
      ),
    (delta.Y() == 0.0)
      ? proj.Y()
      : ((delta.Y() > 0)
        ? -fFrameSize.HalfDepth() // delta positive -> proj on negative side
        : fFrameSize.HalfDepth()
      )
    };

} // geo::PlaneGeo::MoveProjectionToPlane()


//------------------------------------------------------------------------------
geo::Point_t geo::PlaneGeo::MovePointOverPlane(geo::Point_t const& point) const
{

  //
  // This implementation is subject to rounding errors, since the result of
  // the addition might jitter above or below the border.
  //

  auto const deltaProj
    = DeltaFromPlane(PointWidthDepthProjection(point));

  return point + deltaProj.X() * WidthDir<geo::Vector_t>() + deltaProj.Y() * DepthDir<geo::Vector_t>();

} // geo::PlaneGeo::MovePointOverPlane()


//------------------------------------------------------------------------------
void geo::PlaneGeo::UpdateAfterSorting
  (geo::PlaneID const& planeid, geo::BoxBoundedGeo const& TPCbox)
{
  
  UpdateFrameGeometry(planeid, TPCbox);
  
  doUpdateAfterSorting(TPCbox);
  
} // geo::PlaneGeo::UpdateAfterSorting()


//------------------------------------------------------------------------------
std::string geo::PlaneGeo::ViewName(geo::View_t view) {
  switch (view) {
    case geo::kU:       return "U";
    case geo::kV:       return "V";
    case geo::kZ:       return "Z";
    case geo::kY:       return "Y";
    case geo::kX:       return "X";
    case geo::k3D:      return "3D";
    case geo::kUnknown: return "?";
    default:
      return "<UNSUPPORTED (" + std::to_string((int) view) + ")>";
  } // switch
} // PlaneGeo::ViewName()


//------------------------------------------------------------------------------
std::string geo::PlaneGeo::OrientationName(geo::Orient_t orientation) {
  switch (orientation) {
    case geo::kHorizontal: return "horizontal"; break;
    case geo::kVertical:   return "vertical"; break;
    default:               return "unexpected"; break;
  } // switch
} // geo::PlaneGeo::OrientationName()


//------------------------------------------------------------------------------
void geo::PlaneGeo::NotImplemented
  (std::string const& reason /* = "" */, unsigned int stackDepth /* = 5U */)
  const
{
  
  cet::exception e("NotImplemented");
  e << lar::debug::demangle(this)
    << " (from geo::PlaneGeo): call not implemented:\n";
  if (!reason.empty()) e << reason << "\n";
  
  lar::debug::BacktracePrintOptions opts;
  opts.maxLines = stackDepth; // we keep it short...ish
  opts.skipLines = 2U; // skip `printBacktrace()` and `NotImplemented()` calls
  opts.setUniformIndent("  ");
  lar::debug::printBacktrace(e, opts);
  
  throw e;
  
} // geo::PlaneGeo::NotImplemented()


//------------------------------------------------------------------------------
void geo::PlaneGeo::UpdateFrameGeometry
  (geo::PlaneID const& planeid, geo::BoxBoundedGeo const& TPCbox)
{
  // the order here matters

  // reset our ID
  fID = planeid;

  UpdatePlaneNormal(TPCbox);
  UpdateWidthDepthDir();
  UpdateOrientation();
  
} // geo::PlaneGeo::UpdateFrameGeometry()


//------------------------------------------------------------------------------
void geo::PlaneGeo::DetectGeometryDirections() {

  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PlaneGeo` class. If changing *what* is being initialized,
   *      please also update that documentation at the top of `geo::PlaneGeo`
   *      class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`, section
   *      "Initialization details").
   */
  
  /*
   * REMINDER: do not rely on virtual methods from derived classes here, as
   *           they might not be available yet (this method is called by
   *           `geo::PlaneGeo` constructor, when the derived class constructor
   *           hasn't been run yet)
   */
  
  //
  // We need to identify which are the "long" directions of the plane.
  // We assume it is a box, and the shortest side is excluded.
  // The first direction ("width") is given by preference to z.
  // If z is the direction of the normal to the plane... oh well.
  // Let's say privilege to the one which comes from local z, then y.
  // That means: undefined.
  //
  // Requirements:
  //  - ROOT geometry information (shapes and transformations)
  //  - the shape must be a box (an error is PRINTED if not)
  //
  
  assert(fVolume);
  
  //
  // how do they look like in the world?
  //
  TGeoBBox const* pShape = dynamic_cast<TGeoBBox const*>(fVolume->GetShape());
  if (!pShape) {
    mf::LogError("BoxInfo")
      << "Volume " << fVolume->IsA()->GetName() << "['" << fVolume->GetName()
      << "'] has a shape which is a " << pShape->IsA()->GetName()
      << ", not a TGeoBBox! Dimensions won't be available.";
    // set it invalid
    fDecompFrame.SetReferencePoint(geo::origin());
    fDecompFrame.SetMainDir({ 0., 0., 0. });
    fDecompFrame.SetSecondaryDir({ 0., 0., 0. });
    fFrameSize = { 0.0, 0.0 };
    return;
  }

  std::array<geo::Vector_t, 3U> sides;
  size_t iSmallest = 3;
  {

    size_t iSide = 0;

    sides[iSide] = toWorldCoords(LocalVector_t{ pShape->GetDX(), 0.0, 0.0 });
    iSmallest = iSide;
    ++iSide;

    sides[iSide] = toWorldCoords(LocalVector_t{ 0.0, pShape->GetDY(), 0.0 });
    if (sides[iSide].Mag2() < sides[iSmallest].Mag2()) iSmallest = iSide;
    ++iSide;

    sides[iSide] = toWorldCoords(LocalVector_t{ 0.0, 0.0, pShape->GetDZ() });
    if (sides[iSide].Mag2() < sides[iSmallest].Mag2()) iSmallest = iSide;
    ++iSide;

  }

  //
  // which are the largest ones?
  //
  size_t kept[2];
  {
    size_t iKept = 0;
    for (size_t i = 0; i < 3; ++i) if (i != iSmallest) kept[iKept++] = i;
  }

  //
  // which is which?
  //
  // Pick width as the most z-like.
  //
  size_t const iiWidth =
    std::abs(sides[kept[0]].Unit().Z()) > std::abs(sides[kept[1]].Unit().Z())
    ? 0: 1;
  size_t const iWidth = kept[iiWidth];
  size_t const iDepth = kept[1 - iiWidth]; // the other

  // rather than trusting the box, which has been wildly rotated, we pick the
  // width direction to go toward positive _z_ (if possible)
  if (sides[iWidth].Z() < 0.0) sides[iWidth] = -sides[iWidth];
  // this is rhetorical, since will be fixed to have the normal point inward:
  if (sides[iDepth].Y() < 0.0) sides[iDepth] = -sides[iDepth];

  fDecompFrame.SetMainDir(geo::vect::rounded01(sides[iWidth].Unit(), 1e-4));
  fDecompFrame.SetSecondaryDir
    (geo::vect::rounded01(sides[iDepth].Unit(), 1e-4));
  fDecompFrame.SetReferencePoint(GetBoxCenter<geo::Point_t>());
  fFrameSize.halfWidth = sides[iWidth].R();
  fFrameSize.halfDepth = sides[iDepth].R();

} // geo::PlaneGeo::DetectGeometryDirections()


//------------------------------------------------------------------------------
void geo::PlaneGeo::UpdatePlaneNormal(geo::BoxBoundedGeo const& TPCbox) {
  
  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PlaneGeo` class. If changing *what* is being initialized,
   *      please also update that documentation at the top of `geo::PlaneGeo`
   *      class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`, section
   *      "Initialization details").
   */
  
  /*
   * Updates the normal to the plane in `fNormal` (which is independent from
   * the one of the decomposition frames and is the "official" one).
   * 
   * Requirements:
   *  * frame base directions being set and not parallel (just enough to
   *    define the geometric plane they lie on)
   *  * the TPC box needs to be thick enough (see below).
   * 
   * The prescription demands the normal to the plane to face the volume where
   * the drift field is and where the ionization electrons come from.
   * The algorithm requires that the center of the TPC box be in that volume.
   * This turns out to be a problem in "fake" TPCs that host wire planes with
   * "wrapped" wires but no electric drift field. In this case, we require the
   * geometry to have the fake TPC's containing the wrapper wire planes thick
   * enough that the said prescription may be fulfilled.
   * Failure to do so will result in the normal to the plane be defined with the
   * wrong verse, and possibly not consistently between wire planes, with subtle
   * errors later on.
   */
  assert(WidthDir<geo::Vector_t>().Mag2() > 0.0);
  assert(DepthDir<geo::Vector_t>().Mag2() > 0.0);
  
  //
  // direction normal to the plane, points toward the center of TPC
  //

  // start from the direction orthogonal to the plane (verse doesn't matter)
  fNormal = WidthDir<geo::Vector_t>().Cross(DepthDir<geo::Vector_t>());

  // now evaluate where normal should be pointing
  auto const towardCenter = TPCbox.Center() - GetBoxCenter<geo::Point_t>();

  // if they are pointing in opposite directions, flip the normal
  if (fNormal.Dot(towardCenter) < 0) fNormal = -fNormal;
  geo::vect::round01(fNormal, 1e-3);
  
  assert(fNormal.Mag2() != 0.0); // must not be null at this time
  
} // geo::PlaneGeo::UpdatePlaneNormal()


//------------------------------------------------------------------------------
void geo::PlaneGeo::UpdateWidthDepthDir() {

  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PlaneGeo` class. If changing *what* is being initialized,
   *      please also update that documentation at the top of `geo::PlaneGeo`
   *      class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`, section
   *      "Initialization details").
   */
  
  //
  // fix the positiveness of the width/depth/normal frame
  //
  assert(GetNormalDirection<geo::Vector_t>().Mag2() > 0.0);

  // The basis is already set and orthonormal, with only the width
  // and depth verse arbitrary.
  // We choose the direction of the secondary axis ("depth")
  // so that the frame normal is oriented in the general direction of the
  // plane normal (the latter is computed independently).
  if (WidthDir().Cross(DepthDir()).Dot(GetNormalDirection()) < 0.0) {
    fDecompFrame.SetSecondaryDir
      (geo::vect::rounded01(-fDecompFrame.SecondaryDir(), 1e-4));
  }

} // geo::PlaneGeo::UpdateWidthDepthDir()


//------------------------------------------------------------------------------
void geo::PlaneGeo::UpdateOrientation() {

  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PlaneGeo` class. If changing *what* is being initialized,
   *      please also update that documentation at the top of `geo::PlaneGeo`
   *      class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h`, section
   *      "Initialization details").
   */
  
  //
  // this algorithm needs to know about the axis;
  // the normal is expected to be already updated.
  //

  auto const& normal = GetNormalDirection<geo::Vector_t>();
  // sanity check
  if (normal.Mag2() == 0.0) {
    // this likely means construction is not complete yet
    throw cet::exception("Geometry")
      << "geo::PlaneGeo::UpdateOrientation(): no plane normal set yet!\n";
  } // if


  if (std::abs(std::abs(normal.X()) - 1.) < 1e-3)
    fOrientation = kVertical;
  else if (std::abs(std::abs(normal.Y()) - 1.) < 1e-3)
    fOrientation = kHorizontal;
  else {
    // at this point, the only problem is the lack of a label for this
    // orientation; probably introducing a geo::kOtherOrientation would
    // suffice
    throw cet::exception("Geometry")
      << "Plane with unsupported orientation (normal: " << normal << ")\n";
  }

} // geo::PlaneGeo::UpdateOrientation()


//------------------------------------------------------------------------------

