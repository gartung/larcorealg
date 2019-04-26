////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/WirePlaneGeo.h
/// \brief Encapsulate the construction of a single detector wire plane.
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARCOREALG_GEOMETRY_WIREPLANEGEO_H
#define LARCOREALG_GEOMETRY_WIREPLANEGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/SenseWireGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/Decomposer.h"

// C/C++ standard libraries
#include <vector>
#include <string>


namespace geo {

  namespace details {
    struct ActiveAreaCalculator;
  } // namespace details

  //......................................................................

  /**
   * @brief Geometry information for a single wire plane.
   * @ingroup Geometry
   * 
   * The plane is represented in the geometry by a solid which contains wires.
   * The box which is representation of the plane has some thickness, and it
   * should not be assumed that the wires are in the median section of it,
   * that is, the center of the box may not lie on the plane defined by the
   * wires.
   *
   * The wire plane defines two local reference frames.
   * The first, depending on wire directions and therefore called "wire base",
   * is defined by the normal to the plane (pointing toward the center of the
   * TPC), the direction of the wires, and the direction that the wires measure.
   * This is a positive orthogonal base.
   * Note that for this base to be correctly defined, the Geometry service has
   * to provide external information (for example, where the center of the
   * TPC is).
   *
   * The second, depending only on the shape of the plane and called "frame
   * base", is defined by the normal (the same as for the previous one), and two
   * orthogonal axes, "width" and "depth", aligned with the sides of the plane.
   * If the plane has not the shape of a box, this reference frame is not
   * available. This coordinate system is also positive defined.
   * These components are all measured in centimeters.
   *
   */
  class WirePlaneGeo: public geo::PlaneGeo {

  public:
    
    using WireCollection_t = std::vector<std::unique_ptr<geo::SenseWireGeo>>;

    /// Construct a representation of a single plane of the detector
    WirePlaneGeo(
      TGeoNode const& node,
      geo::TransformationMatrix&& trans,
      WireCollection_t&& wires
      );
    
    
    /**
     * @brief Prints information about this plane.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
     * @param indent prepend each line with this string
     * @param verbosity amount of information printed
     *
     * Information on single wires is not printed.
     * Note that the first line out the output is _not_ indented.
     *
     * Verbosity levels
     * -----------------
     *
     * * 0: only plane ID
     * * 1 _(default)_: also center and wire angle
     * * 2: also information about wires
     * * 3: also information about normal and increasing coordinate direction
     * * 4: also information about wire direction, width and depth
     * * 5: also coverage
     * * 6: also bounding box
     *
     */
    template <typename Stream>
    void PrintWirePlaneInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1U)
      const;
    
    
  protected:
    
    // --- BEGIN -- Polymorphic implementation: anode plane --------------------
    /**
     * @name Polymorphic implementation: anode plane.
     */
    /// @{
    
    virtual double doThetaZ() const override;
    virtual double doPhiZ() const override;
    virtual double doSinPhiZ() const override
      { return fSinPhiZ; }
    virtual double doCosPhiZ() const override
      { return fCosPhiZ; }
    virtual unsigned int doNwires() const override
      { return fWire.size(); }
    virtual geo::WireGeo doWire(unsigned int iWire) const override;
    virtual double doWirePitch() const override
      { return fWirePitch; }
    virtual bool doWireIDincreasesWithZ() const override;
    virtual geo::Vector_t doGetIncreasingWireDirection() const override
      { return fDecompWire.SecondaryDir(); }
    virtual geo::Vector_t doGetWireDirection() const override
      { return fDecompWire.MainDir(); }
    virtual geo::WireID doNearestWireID
      (geo::Point_t const& pos) const override;
    virtual geo::WireGeo doNearestWire
      (geo::Point_t const& pos) const override;
    virtual geo::WireID doClosestWireID
      (geo::WireID::WireID_t wireNo) const override;
    virtual lar::util::simple_geo::Volume<> doCoverage() const override;
    virtual double doPlaneCoordinateFrom
      (geo::Point_t const& point, geo::WireGeo const& refWire) const override
      {
        return
          fDecompWire.VectorSecondaryComponent(point - refWire.GetCenter<geo::Point_t>());
      }
    virtual double doPlaneCoordinate(geo::Point_t const& point) const override
      { return fDecompWire.PointSecondaryComponent(point); }
    virtual double doWireCoordinate(geo::Point_t const& point) const override
      { return PlaneCoordinate(point) / WirePitch(); }
    virtual WireDecomposedVector_t doDecomposePoint
      (geo::Point_t const& point) const override
      { return fDecompWire.DecomposePoint(point); }
    virtual geo::Point_t doProjectionReferencePoint() const override
      { return fDecompWire.ReferencePoint(); }
    virtual WireCoordProjection_t doProjection
      (geo::Point_t const& point) const override
      { return fDecompWire.ProjectPointOnPlane(point); }
    virtual WireCoordProjection_t doProjection
      (geo::Vector_t const& v) const override
      { return fDecompWire.ProjectVectorOnPlane(v); }
    virtual geo::Point_t doComposePoint
      (WireDecomposedVector_t const& decomp) const override
      { return fDecompWire.ComposePoint(decomp); }
    virtual geo::Point_t doComposePoint
      (double distance, WireCoordProjection_t const& proj) const override
      { return fDecompWire.ComposePoint(distance, proj); }
    
    virtual std::string doPlaneInfo
      (std::string indent = "", unsigned int verbosity = 1U) const override;

    
    /// Applies the wire sorting algorithm from the specified sorter.
    virtual void doSortElements(geo::GeoObjectSorter const& sorter) override;
    
    /**
     * Called by `geo::PlaneGeo::UpdateAfterSorting()` after setting the new
     * plane ID and geometry decomposition frame.
     * 
     * This implementation defines the wide decomposition base and does plenty
     * of sorting and updates.
     */
    virtual void doUpdateAfterSorting(geo::BoxBoundedGeo const& box) override;
    
    /// @}
    // --- END -- Polymorphic implementation: anode plane ----------------------
    
    
    // --- BEGIN -- Polymorphic implementation: wire abstraction ---------------
    /**
     * @name Polymorphic implementation: wire abstraction.
     * 
     * Since in this implementation we have actual wire objects stored in memory
     * (not `geo::WireGeo` though!) we refer to them for the information.
     */
    /// @{
    
    virtual double doWireRMax(WireLocator const& wloc) const override
      { return getWire(wloc).RMax(); }
    virtual double doWireRMin(WireLocator const& wloc) const override
      { return getWire(wloc).RMin(); }
    virtual double doWireHalfL(WireLocator const& wloc) const override
      { return getWire(wloc).HalfL(); }
    virtual void doWireFillCenterXYZ
      (WireLocator const& wloc, double* xyz, double localz = 0.0) const override
      { getWire(wloc).GetCenter(xyz, localz); }
    virtual void doWireFillStartXYZ
      (WireLocator const& wloc, double* xyz) const override
      { getWire(wloc).GetStart(xyz); }
    virtual void doWireFillEndXYZ
      (WireLocator const& wloc, double* xyz) const override
      { getWire(wloc).GetEnd(xyz); }
    virtual geo::Point_t doWireGetPositionFromCenter
      (WireLocator const& wloc, double localz) const override
      { return getWire(wloc).GetPositionFromCenter<geo::Point_t>(localz); }
    virtual geo::Point_t doWireGetPositionFromCenterUnbounded
      (WireLocator const& wloc, double localz) const override
      {
        return getWire(wloc)
          .GetPositionFromCenterUnbounded<geo::Point_t>(localz);
      }
    virtual geo::Point_t doWireGetCenter(WireLocator const& wloc) const override
      { return getWire(wloc).GetCenter<geo::Point_t>(); }
    virtual geo::Point_t doWireGetStart(WireLocator const& wloc) const override
      { return getWire(wloc).GetStart<geo::Point_t>(); }
    virtual geo::Point_t doWireGetEnd(WireLocator const& wloc) const override
      { return getWire(wloc).GetEnd<geo::Point_t>(); }
    virtual double doWireLength(WireLocator const& wloc) const override
      { return getWire(wloc).Length(); }
    virtual double doWireThetaZ(WireLocator const& wloc) const override
      { return getWire(wloc).ThetaZ(); }
    virtual bool doWireIsParallelTo
      (WireLocator const& wloc, geo::WireGeo const& wire) const override
      { return getWire(wloc).isParallelTo(wire); }
    virtual geo::Vector_t doWireDirection
      (WireLocator const& wloc) const override
      { return getWire(wloc).Direction<geo::Vector_t>(); }
    virtual double doWireDistanceFrom
      (WireLocator const& wloc, geo::WireGeo const& wire) const override
      { return getWire(wloc).DistanceFrom(wire); }
    virtual TGeoNode const* doWireNode(WireLocator const& wloc) const override
      { return getWire(wloc).Node(); }
    virtual double doWireComputeZatY0(WireLocator const& wloc) const override
      { return getWire(wloc).ComputeZatY0(); }
    
    /// Returns the transformation from local to global coordinates.
    /// @returns the transformation object, or `nullptr` if doesn't exist
    virtual WireLocalTransformation_t const* doWireTrans
      (WireLocator const& wloc) const override
      { return &(getWire(wloc).GetTransformation()); }
    
    virtual std::string doWireInfo(
      WireLocator const& wloc,
      std::string indent = "", unsigned int verbosity = 1
      ) const override
      { return getWire(wloc).WireInfo(indent, verbosity); }
    
    virtual void doWireUpdateAfterSorting
      (WireLocator const& wloc, geo::WireID const& wireid, bool flip) override
      { getWire(wloc).UpdateAfterSorting(wireid, flip); }
    
    /// @}
    // --- END -- Polymorphic implementation: wire abstraction -----------------
    
    
  private:
    
    // @{
    /// Return the stored wire for the specified location. Unchecked!
    geo::SenseWireGeo const& getWire(WireLocator const& wloc) const
      { return *(fWire[wloc.wireNo]); }
    geo::SenseWireGeo& getWire(WireLocator const& wloc)
      { return *(fWire[wloc.wireNo]); }
    // @}
    
    /// Returns a direction normal to the plane (pointing is not defined).
    geo::Vector_t GetNormalAxis() const;

    /// Updates the cached direction to increasing wires.
    void UpdateIncreasingWireDir();

    /// Updates the cached direction to wire.
    void UpdateWireDir();

    /// Updates the stored wire pitch.
    void UpdateWirePitch();

    /// Updates the stored wire plane center.
    void UpdateWirePlaneCenter();

    /// Updates the stored @f$ \phi_{z} @f$.
    void UpdatePhiZ();

    /// Updates the stored view
    void UpdateView();

    /// Updates the stored wire pitch with a slower, more robust algorithm.
    void UpdateWirePitchSlow();

    /// Updates the position of the wire coordinate decomposition.
    void UpdateDecompWireOrigin();

    /// Updates the internally used active area.
    void UpdateActiveArea();

    /// Whether the specified wire should have start and end swapped.
    bool shouldFlipWire(geo::SenseWireGeo const& wire) const;

  private:

    WireCollection_t fWire;        ///< List of wires in this plane.

    double           fWirePitch;   ///< Pitch of wires in this plane.
    double           fSinPhiZ;     ///< Sine of @f$ \phi_{z} @f$.
    double           fCosPhiZ;     ///< Cosine of @f$ \phi_{z} @f$.

    /// Decomposition on wire coordinates; the main direction is along the wire,
    /// the secondary one is the one measured by the wire, the normal matches
    /// the plane's normal.
    WireDecomposer_t fDecompWire;

    friend struct details::ActiveAreaCalculator;

  }; // class WirePlaneGeo

} // namespace geo



//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::WirePlaneGeo::PrintWirePlaneInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1U */
) const {
  
  //----------------------------------------------------------------------------
  out << "plane " << std::string(ID());
  
  if (verbosity-- <= 0) return; // 0
  
  //----------------------------------------------------------------------------
  out
    << " at " << GetCenter<geo::Vector_t>() << " cm"
    << ", theta: " << ThetaZ() << " rad";
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  unsigned int const nWires = Nwires();
  
  out << "\n" << indent
    << "normal to wire: " << PhiZ() << " rad"
      << ", with orientation " << OrientationName(Orientation())
      << ", has " << nWires << " wires measuring " << ViewName(View())
      << " with a wire pitch of " << WirePitch() << " cm"
    ;
  
  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  auto const& normal = GetNormalDirection<geo::Vector_t>();
  auto const& incrZdir = GetIncreasingWireDirection<geo::Vector_t>();
  auto const& wireNormalDir = fDecompWire.NormalDir();
  out << "\n" << indent
    << "normal to plane: " << normal
    << ", direction of increasing wire number: " << incrZdir
    << " [wire frame normal: " << wireNormalDir << "]"
    << " (" << (WireIDincreasesWithZ()? "increases": "decreases") << " with z)";
    
  if (verbosity-- <= 0) return; // 3
  
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
    
  if (verbosity-- <= 0) return; // 4
  
  //----------------------------------------------------------------------------
  // get the area spanned by the wires
  out << "\n" << indent << "wires cover width "
    << ActiveArea().width.lower << " to " << ActiveArea().width.upper
    << ", depth "
    << ActiveArea().depth.lower << " to " << ActiveArea().depth.upper
    << " cm";
  if (verbosity-- <= 0) return; // 5
  
  //----------------------------------------------------------------------------
  // print also the containing box
  auto const box = BoundingBox();
  out << "\n" << indent
    << "bounding box: " << box.Min() << " -- " << box.Max();
  
//  if (verbosity-- <= 0) return; // 6
  
  //----------------------------------------------------------------------------
} // geo::PlaneGeo::PrintWirePlaneInfo()


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_WIREPLANEGEO_H
////////////////////////////////////////////////////////////////////////
