////////////////////////////////////////////////////////////////////////
/// \file larcorealg/Geometry/WirePlaneGeo.cxx
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/WirePlaneGeo.h"

// LArSoft includes
#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::convertTo()
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/SortByPointers.h" // util::makePointerVector()
#include "larcorealg/CoreUtils/counter.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// ROOT includes
#include "TMath.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TGeoMatrix.h"
#include "TClass.h"

// C/C++ standard library
#include <sstream> // std::ostringstream
#include <array>
#include <functional> // std::less<>, std::greater<>, std::transform()
#include <iterator> // std::back_inserter()
#include <utility> // std::move()
#include <type_traits> // std::is_same<>, std::decay_t<>
#include <cassert>



namespace geo{

  namespace details {

    //......................................................................
    /**
     * @brief Class computing the active area of the plane.
     *
     * The active area is defined in the width/depth space which include
     * approximatively all wires.
     *
     * This algorithm requires the frame reference and the wire pitch to be
     * already defined.
     *
     * That area is tuned so that all its points are closer than half a wire
     * pitch from a wire.
     *
     */
    struct ActiveAreaCalculator {

      ActiveAreaCalculator
        (geo::WirePlaneGeo const& plane, double wMargin, double dMargin)
        : plane(plane)
        , wMargin(wMargin)
        , dMargin(dMargin)
        {}

      ActiveAreaCalculator(geo::WirePlaneGeo const& plane, double margin = 0.0)
        : ActiveAreaCalculator(plane, margin, margin)
        {}

      operator geo::WirePlaneGeo::Rect()
        { return recomputeArea(); }

        private:
      using Projection_t = ROOT::Math::PositionVector2D
        <ROOT::Math::Cartesian2D<double>, geo::WirePlaneGeo::WidthDepthReferenceTag>;
      using Vector_t = geo::WirePlaneGeo::WidthDepthDisplacement_t;

      static_assert(
        !std::is_same<Projection_t, geo::WirePlaneGeo::WidthDepthProjection_t>::value,
        "Necessary maintenance: remove the now optional conversions"
        );

      static constexpr std::size_t kFirstWireStart = 0;
      static constexpr std::size_t kFirstWireEnd   = 1;
      static constexpr std::size_t kLastWireStart  = 2;
      static constexpr std::size_t kLastWireEnd    = 3;

      geo::WirePlaneGeo const& plane; ///< Plane to work on.
      double const wMargin = 0.0; ///< Margin subtracted from each side of width.
      double const dMargin = 0.0; ///< Margin subtracted from each side of depth.
      geo::WirePlaneGeo::Rect activeArea; ///< Result.

      /// Cache: wire end projections.
      Projection_t wireEnds[4];

      void initializeWireEnds()
        {
          //
          // Collect the projections of the relevant points.
          //
          // Sorted so that start points have width not larger than end points.
          //
          // PointWidthDepthProjection() erroneously returns a vector rather
          // than a point, so a conversion is required
          auto makeProjection
            = [](auto v){ return Projection_t(v.X(), v.Y()); };

          wireEnds[kFirstWireStart] = makeProjection
            (plane.PointWidthDepthProjection(plane.FirstWire().GetStart()));
          wireEnds[kFirstWireEnd] = makeProjection
            (plane.PointWidthDepthProjection(plane.FirstWire().GetEnd()));
          if (wireEnds[kFirstWireStart].X() > wireEnds[kFirstWireEnd].X())
            std::swap(wireEnds[kFirstWireStart], wireEnds[kFirstWireEnd]);
          wireEnds[kLastWireStart] = makeProjection
            (plane.PointWidthDepthProjection(plane.LastWire().GetStart()));
          wireEnds[kLastWireEnd] = makeProjection
            (plane.PointWidthDepthProjection(plane.LastWire().GetEnd()));
          if (wireEnds[kLastWireStart].X() > wireEnds[kLastWireEnd].X())
            std::swap(wireEnds[kLastWireStart], wireEnds[kLastWireEnd]);
        } // initializeWireEnds()

      void includeAllWireEnds()
        {
          //
          // Find the basic area containing all the coordinates.
          //

          // include all the coordinates of the first and last wire
          for (auto const& aWireEnd: wireEnds) {
            activeArea.width.extendToInclude(aWireEnd.X());
            activeArea.depth.extendToInclude(aWireEnd.Y());
          }

        } // includeAllWireEnds()

      void adjustCorners()
        {
          //
          // Modify the corners so that none is father than half a pitch from all
          // wires.
          //
          // directions in wire/depth plane
          Vector_t const widthDir = { 1.0, 0.0 };
          Vector_t const depthDir = { 0.0, 1.0 };
          Vector_t wireCoordDir = plane.VectorWidthDepthProjection
            (plane.GetIncreasingWireDirection());
          double const hp = plane.WirePitch() / 2.0; // half pitch

          //
          // The plan: identify if wires are across or corner, and then:
          // - across:
          //   - identify which sides
          //   - set the farther end of the wire from the side to be p/2 from
          //     its corner
          // - corner:
          //   - identify which corners
          //   - move the corners to p/2 from the wires
          //

          //
          // are the wires crossing side to side, as opposed to cutting corners?
          //

          // these are the angles of the original wire coordinate direction
          double const cosAngleWidth = geo::vect::dot(wireCoordDir, widthDir);
          double const cosAngleDepth = geo::vect::dot(wireCoordDir, depthDir);
          // if the wire coordinate direction is on first or third quadrant:
          bool const bPositiveAngle
            = none_or_both((wireCoordDir.X() >= 0), (wireCoordDir.Y() >= 0));

          // now we readjust the wire coordinate direction to always point
          // toward positive width; this breaks the relation between
          // wireCoordDir and which is the first/last wire
          if (cosAngleWidth < 0) wireCoordDir = -wireCoordDir;

          // let's study the first wire (ends are sorted by width)
          assert(wireEnds[kFirstWireEnd].X() >= wireEnds[kFirstWireStart].X());
          bool const bAlongWidth // horizontal
            =  equal(wireEnds[kFirstWireEnd].X(), activeArea.width.upper)
            && equal(wireEnds[kFirstWireStart].X(), activeArea.width.lower);
          bool const bAlongDepth = !bAlongWidth && // vertical
               equal(std::max(wireEnds[kFirstWireStart].Y(), wireEnds[kFirstWireEnd].Y()), activeArea.depth.upper)
            && equal(std::min(wireEnds[kFirstWireStart].Y(), wireEnds[kFirstWireEnd].Y()), activeArea.depth.lower);
          assert(!(bAlongWidth && bAlongDepth));

          if (bAlongWidth) { // horizontal

            // +---------+
            // |   ___,--|  upper width bound
            // |--'      |

            // find which is the wire with higher width:
            // the last wire is highest if the wire coordinate direction (which
            // is defined by what is first and what is last) is parallel to the
            // width direction
            std::size_t const iUpperWire
              = (cosAngleDepth > 0)? kLastWireStart: kFirstWireStart;
            // largest distance from upper depth bound of the two ends of wire
            double const maxUpperDistance = activeArea.depth.upper
              - std::min
                (wireEnds[iUpperWire].Y(), wireEnds[iUpperWire ^ 0x1].Y())
              ;
            // set the upper side so that the maximum distance is p/2
            // (it may be actually less if the wire is not perfectly horizontal)
            activeArea.depth.upper += (hp - maxUpperDistance);

            // |--.___   |
            // |      `--|  deal with the lower bound now
            // +---------+

            std::size_t const iLowerWire
              = (cosAngleDepth > 0)? kFirstWireStart: kLastWireStart;
            // largest distance from lower depth bound of the two ends of wire
            double const maxLowerDistance
              = std::max
                (wireEnds[iLowerWire].Y(), wireEnds[iLowerWire ^ 0x1].Y())
              - activeArea.depth.lower
              ;
            // set the upper side so that the minimum distance is p/2
            activeArea.depth.lower -= (hp - maxLowerDistance);

          } // horizontal wires
          else if (bAlongDepth) { // vertical
            // --,---+
            //   |   |
            //    \  |
            //    |  |  upper depth bound
            //     \ |
            //     | |
            // ------+

            // find which is the wire with higher depth:
            // the last wire is highest if the wire coordinate direction (which
            // is defined by what is first and what is last) is parallel to the
            // depth direction
            std::size_t const iUpperWire
              = (cosAngleWidth > 0)? kLastWireStart: kFirstWireStart;
            // largest distance from upper depth bound of the two ends of wire
            double const maxUpperDistance = activeArea.width.upper
              - std::min
                (wireEnds[iUpperWire].X(), wireEnds[iUpperWire ^ 0x1].X())
              ;
            // set the upper side so that the minimum distance is p/2
            activeArea.width.upper += (hp - maxUpperDistance);

            // +-,----
            // | |
            // |  \     .
            // |  |     deal with the lower bound now
            // |   \    .
            // |   |
            // +------
            std::size_t const iLowerWire
              = (cosAngleWidth > 0)? kFirstWireStart: kLastWireStart;
            // largest distance from lower width bound of the two ends of wire
            double const maxLowerDistance
              = std::max
                (wireEnds[iLowerWire].X(), wireEnds[iLowerWire ^ 0x1].X())
              - activeArea.width.lower
              ;
            // set the upper side so that the minimum distance is p/2
            activeArea.width.lower -= (hp - maxLowerDistance);

          } // vertical wires
          else if (bPositiveAngle) { // wires are not going across: corners!
            // corners at (lower width, lower depth), (upper width, upper depth)

            // -._------+
            //    `-._  |  upper width corner (upper depth)
            //        `-|

            // start of the wire on the upper corner
            // (width coordinate is lower for start than for end)
            std::size_t const iUpperWire
              = (cosAngleWidth > 0)? kLastWireStart: kFirstWireStart;

            double const upperDistance = geo::vect::dot(
              Vector_t(activeArea.width.upper - wireEnds[iUpperWire].X(), 0.0),
              wireCoordDir
              );
            // make the upper distance become p/2
            auto const upperDelta = (hp - upperDistance) * wireCoordDir;
            activeArea.width.upper += upperDelta.X();
            activeArea.depth.upper += upperDelta.Y();

            // |-._
            // |   `-._    lower width corner (lower depth)
            // +-------`-

            // end of the wire on the lower corner
            // (width coordinate is lower than the end)
            std::size_t const iLowerWire
              = (cosAngleWidth > 0)? kFirstWireEnd: kLastWireEnd;

            double const lowerDistance = geo::vect::dot(
              Vector_t(wireEnds[iLowerWire].X() - activeArea.width.lower, 0.0),
              wireCoordDir
              );
            // make the lower distance become p/2 (note direction of wire coord)
            auto const lowerDelta = (hp - lowerDistance) * wireCoordDir;
            activeArea.width.lower -= lowerDelta.X();
            activeArea.depth.lower -= lowerDelta.Y();

          }
          else { // !bPositiveAngle
            // corners at (lower width, upper depth), (upper width, lower depth)

            //       _,-|
            //   _,-'   |  upper width corner (lower depth)
            // -'-------+

            // start of the wire on the upper corner
            // (width coordinate is lower than the end)
            std::size_t const iUpperWire
              = (cosAngleWidth > 0)? kLastWireStart: kFirstWireStart;

            double const upperDistance = geo::vect::dot(
              Vector_t(activeArea.width.upper - wireEnds[iUpperWire].X(), 0.0),
              wireCoordDir
              );
            // make the upper distance become p/2
            auto const upperDelta = (hp - upperDistance) * wireCoordDir;
            activeArea.width.upper += upperDelta.X();
            activeArea.depth.lower += upperDelta.Y();

            // +------_,-
            // |  _,-'     lower width corner (upper depth)
            // |-'

            // end of the wire on the lower corner
            // (width coordinate is lower than the end)
            std::size_t const iLowerWire
              = (cosAngleWidth > 0)? kFirstWireEnd: kLastWireEnd;

            double const lowerDistance = geo::vect::dot(
              Vector_t(wireEnds[iLowerWire].X() - activeArea.width.lower, 0.0),
              wireCoordDir
              );
            // make the lower distance become p/2 (note direction of wire coord)
            auto const lowerDelta = (hp - lowerDistance) * wireCoordDir;
            activeArea.width.lower -= lowerDelta.X();
            activeArea.depth.upper -= lowerDelta.Y();

          } // if ...

        } // adjustCorners()


      void applyMargin()
        {
          if (wMargin != 0.0) {
            activeArea.width.lower += wMargin;
            activeArea.width.upper -= wMargin;
          }
          if (dMargin != 0.0) {
            activeArea.depth.lower += dMargin;
            activeArea.depth.upper -= dMargin;
          }
        } // applyMargin()


      geo::WirePlaneGeo::Rect recomputeArea()
        {
          activeArea = {};

          //
          // 0. collect the projections of the relevant point
          //
          initializeWireEnds();

          //
          // 1. find the basic area containing all the coordinates
          //
          includeAllWireEnds();

          //
          // 2. adjust area so that no corner is father than half a wire pitch
          //
          adjustCorners();

          //
          // 3. apply an absolute margin
          //
          applyMargin();

          return activeArea;
        } // computeArea()

      /// Returns true if a and b are both true or both false (exclusive nor).
      static bool none_or_both(bool a, bool b) { return a == b; }

      /// Returns whether the two numbers are the same, lest a tolerance.
      template <typename T>
      static bool equal(T a, T b, T tol = T(1e-5))
        { return std::abs(a - b) <= tol; }

    }; // struct ActiveAreaCalculator

    //......................................................................

  } // namespace details


  //----------------------------------------------------------------------------
  //---  geo::WirePlaneGeo
  //---
  WirePlaneGeo::WirePlaneGeo(
    TGeoNode const& node,
    geo::TransformationMatrix&& trans,
    WireCollection_t&& wires
    )
    : geo::PlaneGeo(node, std::move(trans))
    , fWire(std::move(wires))
    , fWirePitch(0.)
    , fSinPhiZ(0.)
    , fCosPhiZ(0.)
    , fDecompWire()
  {
    
    UpdateWirePitchSlow();
    
  } // WirePlaneGeo::WirePlaneGeo()

  
  //......................................................................

  geo::WireGeo WirePlaneGeo::doWire(unsigned int iWire) const {
    
    if (!doHasWire(iWire)) {
      throw cet::exception("WireOutOfRange")
        << "Request for non-existing wire " << iWire << " on plane " << ID()
        << " (there are " << doNwires() << ")\n";
    }
    return { *this, iWire };
  } // WirePlaneGeo::doWire(int)

  //......................................................................

  // sort the WireGeo objects
  void WirePlaneGeo::doSortElements(geo::GeoObjectSorter const& sorter )
  {
    /*
     * The design of geometry sorting is one of the things in the geometry
     * system that has given me the most grief.
     * Each time something is changed in the underlying representation, it has
     * catastrophic consequences on that side.
     * The following is a very expensive and very generic approach to the task.
     * It is recommended that rather than coping with this an experiment derives
     * their own plane doing the sorting in its best way.
     * Maybe a bit overkill (it needs a custom builder too), but quite
     * effective.
     */
    
    /*
     * (1) prepare a vector of `geo::WireGeo`
     * (2) prepare a vector of _pointers_ to `geo::WireGeo`
     * (3) sort that vector with the geometry sorter
     * (4) move the `geo::SenseWireGeo` objects in a new container, using
     *     the wire number from `geo::WireGeo` ID as new position
     * (5) replace the official vector with the one just created
     * 
     * Three more vectors, for Bjarne sake!!
     */
    
    // (1) prepare a vector of `geo::WireGeo`
    std::vector<geo::WireGeo> storage;
    storage.reserve(fWire.size());
    for (std::size_t index: util::counter(fWire.size()))
      storage.push_back(doWire(index));
    
    // (2) prepare a vector of _pointers_ to `geo::WireGeo`
    //     (`SortWires()` interface requires non-const pointers)
    std::vector<geo::WireGeo*> wirePtrs = util::makePointerVector(storage);
    
    // (3) sort that vector with the geometry sorter
    sorter.SortWires(wirePtrs);
    
    // (4) move the `geo::SenseWireGeo` objects in a new container, using
    //     the wire number from `geo::WireGeo` ID as new position
    WireCollection_t wires;
    wires.reserve(wirePtrs.size());
    for (geo::WireGeo const* wirePtr: wirePtrs)
      wires.push_back(std::move(fWire[wirePtr->ID().Wire]));
    
    // (5) replace the official vector with the one just created
    fWire = std::move(wires);
    
  } // doSortElements()


  //......................................................................
  bool WirePlaneGeo::doWireIDincreasesWithZ() const {
    return GetIncreasingWireDirection().Z() > 0.;
  } // WirePlaneGeo::doWireIDincreasesWithZ()


  //......................................................................
  lar::util::simple_geo::Volume<> WirePlaneGeo::doCoverage() const {

    // add both coordinates of first and last wire
    std::array<double, 3> A, B;

    FirstWire().GetStart(A.data());
    LastWire().GetEnd(B.data());

    return { A.data(), B.data() };
  } // WirePlaneGeo::doCoverage()


  //......................................................................
  geo::WireID WirePlaneGeo::doNearestWireID(geo::Point_t const& pos) const {

    //
    // 1) compute the wire coordinate of the point
    // 2) get the closest wire number
    // 3) check if the wire does exist
    // 4) build and return the wire ID
    //

    // this line merges parts (1) and (2); add 0.5 to have the correct rounding:
    int nearestWireNo = int(0.5 + WireCoordinate(pos));

    // if we are outside of the wireplane range, throw an exception
    if ((nearestWireNo < 0) || ((unsigned int) nearestWireNo >= Nwires())) {

      auto wireNo = nearestWireNo; // save for the output

      if (nearestWireNo < 0 ) wireNo = 0;
      else                    wireNo = Nwires() - 1;

      throw InvalidWireError("Geometry", ID(), nearestWireNo, wireNo)
        << "Can't find nearest wire for position " << pos
        << " in plane " << std::string(ID()) << " approx wire number # "
        << wireNo << " (capped from " << nearestWireNo << ")\n";
    } // if invalid

    return { ID(), (geo::WireID::WireID_t) nearestWireNo };

  } // WirePlaneGeo::doNearestWireID()


  //......................................................................
  geo::WireGeo WirePlaneGeo::doNearestWire(geo::Point_t const& point) const {

    //
    // Note that this code is ready for when NearestWireID() will be changed
    // to return an invalid ID instead of throwing.
    // As things are now, `NearestWireID()` will never return an invalid ID,
    // but it will throw an exception similar to this one.
    //

    geo::WireID const wireID = NearestWireID(point);
    if (wireID) return Wire(wireID); // we have that wire, so we return it

    // wire ID is invalid, meaning it's out of range. Throw an exception!
    geo::WireID const closestID = ClosestWireID(wireID);
    throw InvalidWireError("Geometry", ID(), closestID.Wire, wireID.Wire)
      << "Can't find nearest wire for position " << point
      << " in plane " << std::string(ID()) << " approx wire number # "
      << closestID.Wire << " (capped from " << wireID.Wire << ")\n";

  } // WirePlaneGeo::doNearestWire()


  //......................................................................
  geo::WireID geo::WirePlaneGeo::doClosestWireID
    (geo::WireID::WireID_t wireNo) const
  {
    return { ID(), boundedValue<geo::WireID::WireID_t>(wireNo, 0, Nwires()) };
  }

  
  //......................................................................
  double WirePlaneGeo::doThetaZ() const { return FirstWire().ThetaZ(); }

  //......................................................................
  double WirePlaneGeo::doPhiZ() const { return std::atan2(fSinPhiZ, fCosPhiZ); }

  //......................................................................
  void WirePlaneGeo::doUpdateAfterSorting(geo::BoxBoundedGeo const& TPCbox) {
    
    // the order here matters

    UpdateIncreasingWireDir();

    // update wires
    geo::WireID::WireID_t wireNo = 0;
    for (auto& wire: fWire) {

      wire->UpdateAfterSorting(geo::WireID(fID, wireNo), shouldFlipWire(*wire));

      ++wireNo;
    } // for wires

    UpdateDecompWireOrigin();
    UpdateWireDir();
    UpdateWirePlaneCenter();
    UpdateWirePitch();
    UpdateActiveArea();
    UpdatePhiZ();
    UpdateView();

  } // WirePlaneGeo::doUpdateAfterSorting()

  //......................................................................
  std::string WirePlaneGeo::doPlaneInfo
    (std::string indent /* = "" */, unsigned int verbosity /* = 1U */) const
  {
    std::ostringstream out;
    PrintWirePlaneInfo(out, indent, verbosity);
    return out.str();
  } // WirePlaneGeo::doPlaneInfo()
  
  
  //......................................................................
  geo::Vector_t WirePlaneGeo::GetNormalAxis() const {
    
    //
    // This used to be called by the implementation of UpdatePlaneNormal().
    // Now that implementation is in the base class and this one is pretty
    // much unused.
    //
    
    const unsigned int NWires = Nwires();
    if (NWires < 2) return {}; // why are we even here?

    // 1) get the direction of the middle wire
    auto const WireDir = Wire(NWires / 2).Direction<geo::Vector_t>();

    // 2) get the direction between the middle wire and the next one
    auto const ToNextWire = Wire(NWires / 2 + 1).GetCenter<geo::Point_t>()
      - Wire(NWires / 2).GetCenter<geo::Point_t>();

    // 3) get the direction perpendicular to the plane
    // 4) round it
    // 5) return its norm
    return geo::vect::rounded01(WireDir.Cross(ToNextWire).Unit(), 1e-4);

  } // WirePlaneGeo::GetNormalAxis()


  //......................................................................
  void WirePlaneGeo::UpdateWirePitch() {
    fWirePitch = geo::WireGeo::WirePitch(Wire(0), Wire(1));
  } // WirePlaneGeo::UpdateWirePitch()

  //......................................................................
  void WirePlaneGeo::UpdatePhiZ() {
    TVector3 const& wire_coord_dir = GetIncreasingWireDirection();
  /*
    TVector3 const& normal = GetNormalDirection();
    TVector3 z(0., 0., 1.);

    // being defined in absolute terms as angle respect to z axis,
    // we take the z component as cosine, and all the rest as sine
    fCosPhiZ = wire_coord_dir.Dot(z);
    fSinPhiZ = wire_coord_dir.Cross(z).Dot(normal);
  */
    fCosPhiZ = wire_coord_dir.Z();
    fSinPhiZ = wire_coord_dir.Y();
  } // WirePlaneGeo::UpdatePhiZ()

  void WirePlaneGeo::UpdateView()
  {
    /*
     * This algorithm assigns views according to the angle the wire axis cuts
     * with y axis ("thetaY"), but from the point of view of the center of the
     * TPC.
     * A special case is when the drift axis is on y axis.
     *
     * In the normal case, the discrimination happens on the the arctangent of
     * the point { (y,w), (y x n,w) }, where w is the wire direction, y is the
     * coordinate axis and n the normal to the wire plane. This definition gives
     * the same value regardless of the direction of w on its axis.
     *
     * If thetaY is 0, wires are parallel to the y axis:
     * the view is assigned as kX or kZ depending on whether the plane normal is
     * closer to the z axis or the x axis, respectively (the normal describes
     * a direction _not_ measured by the wires).
     *
     * If thetaY is a right angle, the wires are orthogonal to y axis and view
     * kY view is assigned.
     * If thetaY is smaller than 0, the view is called "U".
     * If thetaY is larger than 0, the view is called "V".
     *
     * The special case where the drift axis is on y axis is treated separately.
     * In that case, the role of y axis is replaced by the z axis and the
     * discriminating figure is equivalent to the usual ThetaZ().
     *
     * If thetaZ is 0, the wires are measuring x and kX view is chosen.
     * If thetaZ is a right angle, the wires are measuring z and kZ view is
     * chosen.
     * If thetaZ is smaller than 0, the view is called "U".
     * If thetaZ is larger than 0, the view is called "V".
     *
     */

    auto const& normalDir = GetNormalDirection<geo::Vector_t>();
    auto const& wireDir = GetWireDirection<geo::Vector_t>();

    // normal direction has been rounded, so exact comparison can work
    if (std::abs(normalDir.Y()) != 1.0) {
      //
      // normal case: drift direction is not along y (vertical)
      //

      // yw is pretty much GetWireDirection().Y()...
      // thetaY is related to atan2(ynw, yw)
      double const yw = geo::vect::dot(wireDir, geo::Yaxis());
      double const ynw = geo::vect::mixedProduct
        (geo::Yaxis(), normalDir, wireDir);

      if (std::abs(yw) < 1.0e-4) { // wires orthogonal to y axis
        double const closeToX
          = std::abs(geo::vect::dot(normalDir, geo::Xaxis()));
        double const closeToZ
          = std::abs(geo::vect::dot(normalDir, geo::Zaxis()));
        SetView((closeToZ > closeToX)? geo::kX: geo::kY);
      }
      else if (std::abs(ynw) < 1.0e-4) { // wires parallel to y axis
        SetView(geo::kZ);
      }
      else if ((ynw * yw) < 0) SetView(geo::kU); // different sign => thetaY > 0
      else if ((ynw * yw) > 0) SetView(geo::kV); // same sign => thetaY < 0
      else assert(false); // logic error?!

    }
    else { // if drift is vertical
      //
      // special case: drift direction is along y (vertical)
      //

      // zw is pretty much GetWireDirection().Z()...
      double const zw = geo::vect::dot(wireDir, geo::Zaxis());
      // while GetNormalDirection() axis is on y, its direction is not fixed:
      double const znw = geo::vect::mixedProduct
        (geo::Zaxis(), normalDir, wireDir);

      // thetaZ is std::atan(znw/zw)

      if (std::abs(zw) < 1.0e-4) { // orthogonal to z, orthogonal to y...
        // this is equivalent to thetaZ = +/- pi/2
        SetView(geo::kZ);
      }
      else if (std::abs(znw) < 1.0e-4) { // parallel to z, orthogonal to y...
        // this is equivalent to thetaZ = 0
        SetView(geo::kX);
      }
      else if ((znw * zw) < 0) SetView(geo::kU); // different sign => thetaZ > 0
      else if ((znw * zw) > 0) SetView(geo::kV); // same sign => thetaZ < 0
      else assert(false); // logic error?!

    } // if drift direction... else

  } // UpdateView()


  //......................................................................
  void WirePlaneGeo::UpdateIncreasingWireDir() {

    //
    // Direction measured by the wires, pointing toward increasing wire number;
    // requires:
    // - the normal to the plane to be correct
    // - wires to be sorted
    //

    // 1) get the direction of the middle wire
    auto refWireNo = Nwires() / 2;
    if (refWireNo == Nwires() - 1) --refWireNo;
    auto const& refWire = Wire(refWireNo);
    auto const& WireDir = geo::vect::toVector(refWire.Direction()); // we only rely on the axis


    // 2) get the axis perpendicular to it on the wire plane
    //    (arbitrary direction)
    auto wireCoordDir = GetNormalDirection<geo::Vector_t>().Cross(WireDir).Unit();

    // 3) where is the next wire?
    auto toNextWire
      = geo::vect::toVector(Wire(refWireNo + 1).GetCenter() - refWire.GetCenter());

    // 4) if wireCoordDir is pointing away from the next wire, flip it
    if (wireCoordDir.Dot(toNextWire) < 0) {
      wireCoordDir = -wireCoordDir;
    }
    fDecompWire.SetSecondaryDir(geo::vect::rounded01(wireCoordDir, 1e-4));

  } // WirePlaneGeo::UpdateIncreasingWireDir()

  
  //......................................................................
  void WirePlaneGeo::UpdateWireDir() {

    fDecompWire.SetMainDir(geo::vect::rounded01(geo::vect::toVector(FirstWire().Direction()), 1e-4));
    
    //
    // check that the resulting normal matches the plane one
    //
    assert(lar::util::makeVector3DComparison(1e-5)
      .equal(fDecompWire.NormalDir(), GetNormalDirection<geo::Vector_t>()));

  } // WirePlaneGeo::UpdateWireDir()


  //......................................................................
  void WirePlaneGeo::UpdateWirePitchSlow() {

    //
    // Compare one wire (the first one, for convenience) with all other wires;
    // the wire pitch is the smallest distance we find.
    //
    // This algorithm assumes wire pitch is constant, but it does not assume
    // wire ordering (which UpdateWirePitch() does).
    //
    auto firstWire = fWire.cbegin(), wire = firstWire, wend = fWire.cend();
    fWirePitch = geo::SenseWireGeo::WirePitch(**firstWire, **(++wire));

    while (++wire != wend) {
      auto wirePitch = geo::SenseWireGeo::WirePitch(**firstWire, **wire);
      if (wirePitch < 1e-4) continue; // it's 0!
      if (wirePitch < fWirePitch) fWirePitch = wirePitch;
    } // while

  } // WirePlaneGeo::UpdateWirePitchSlow()


  //......................................................................
  void WirePlaneGeo::UpdateDecompWireOrigin() {

    //
    // update the origin of the reference frame (the middle of the first wire)
    //
    fDecompWire.SetOrigin(geo::vect::toPoint(FirstWire().GetCenter()));

  } // WirePlaneGeo::UpdateDecompWireOrigin()

  //......................................................................
  void WirePlaneGeo::UpdateActiveArea() {

    //
    // The active area is defined in the width/depth space which include
    // approximatively all wires.
    //
    // See `ActiveAreaCalculator` for details of the algorithm.
    //

    // we scratch 1 um from each side to avoid rounding errors later
    fActiveArea = details::ActiveAreaCalculator(*this, 0.0001);

  } // WirePlaneGeo::UpdateActiveArea()


  //......................................................................
  void WirePlaneGeo::UpdateWirePlaneCenter() {

    //
    // The center of the wire plane is defined as the center of the plane box,
    // translated to the plane the wires lie on.
    // This assumes that the thickness direction of the box is aligned with
    // the drift direction, so that the translated point is still in the middle
    // of width and depth dimensions.
    // It is possible to remove that assumption by translating the center of the
    // box along the thickness direction enough to bring it to the wire plane.
    // The math is just a bit less straightforward, so we don't bother yet.
    //
    // Requirements:
    //  * the wire decomposition frame must be set up (at least its origin and
    //    normal direction)
    //

    fCenter = GetBoxCenter<geo::Point_t>();

    DriftPoint(fCenter, fDecompWire.PointNormalComponent(fCenter));
    
    geo::vect::round0(fCenter, 1e-7); // round dimensions less than 1 nm to 0
    
    fDecompFrame.SetOrigin(fCenter); // equivalent to GetCenter() now

  } // WirePlaneGeo::UpdateWirePlaneCenter()


  //......................................................................
  bool WirePlaneGeo::shouldFlipWire(geo::SenseWireGeo const& wire) const {
    //
    // The correct orientation is so that:
    //
    // (direction) x (wire coordinate direction) . (plane normal)
    //
    // is positive; it it's negative, then we should flip the wire.
    //
    // Note that the increasing wire direction comes from the wire frame, while
    // the normal direction is computed independently by geometry.
    // The resulting normal in the wire frame is expected to be the same as the
    // plane normal from GetNormalDirection(); if this is not the case, flipping
    // the wire direction should restore it.
    //

    return wire.Direction()
      .Cross(GetIncreasingWireDirection())
      .Dot(GetNormalDirection())
      < +0.5; // should be in fact exactly +1

  } // WirePlaneGeo::shouldFlipWire()

  //......................................................................


} // namespace geo
////////////////////////////////////////////////////////////////////////
