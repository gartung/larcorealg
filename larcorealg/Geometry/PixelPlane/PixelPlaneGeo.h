/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeo.h
 * @brief   Representation of a readout plane with pixels as sensitive elements.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 23, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeo.cxx`
 * @ingroup Geometry
 */

#ifndef LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H
#define LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/PixelPlane/PixelCoord.h" // geo::pixel::PixelCoordT<>
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/Decomposer.h"

// C/C++ standard libraries
#include <string>
#include <array>
#include <limits> // std::numeric_limits<>
#include <utility> // std::pair<>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace geo {
  class PixelPlaneGeo;
} // namespace geo


/**
 * @brief Geometry information for a single pixel plane.
 * @ingroup Geometry
 * 
 * 
 * Interface as a wire plane
 * ==========================
 * 
 * This object provides a compatibility layer answering to calls that are worded
 * in terms of TPC wires. In translating the wire language into the pixel one,
 * the following choices are made:
 * 
 * * the wire coordinate (i.e. the only coordinate that a wire plane claims to
 *     measure) is translated into one of the two coordinates the pixel plane
 *     measures, namely the "main" coordinate. In this sense, the pixel plane
 *     can be seen as a wire plane where each wire is made of a row of pixels,
 *     all at the same secondary coordinate (see also the definitions section
 *     below);
 * * the angle @f$ \varphi_{z} @f$ (`geo::PlaneGeo::PhiZ()`) matches the
 *     secondary direction of the wire frame; this is a consequence of the
 *     definitions of the two entities.
 * 
 * 
 * Definitions and initialization
 * ===============================
 * 
 * 
 * The plane is represented in the geometry by a box.
 * 
 * If the box contains a volume with the name in
 * `geo::PixelPlaneGeo::SensitivePixelPlaneVolumeName`, this volume is used to
 * determine the sensitive area, as follows:
 * 
 * * the volume is required to be a box, i.e. with three orthogonal sides
 * * it is expected and assumed that the plane is not orthogonal to the beam
 *   direction, defined in LArSoft as _z_;
 * * the shortest side is considered the depth of the pixels, and ignored;
 * * the other of the two sides most aligned to _z_ defines the secondary
 *   direction of the plane; that direction is chosen toward positive _z_;
 *   the remaining side defines the main direction and its direction is taken so
 *   that the tern main, secondary and normal directions is positive defined,
 *   following the standard `geo::PlaneGeo` convention (that is also the
 *   convention imposing the secondary direction to be the one most oriented
 *   toward _z_ and therefore the one the wires measure);
 * * these two sides define the borders of the pixels;
 *   pixels are defined aligned to these sides and to cover the entire area
 *   with the possible exception of a slice at the far end too narrow to host
 *   a further row of pixels;
 * * the pixels are considered to be thickless and laid on the face of the box
 *   closest to the cathode
 * * the number and size of pixels is determined TODO
 * * the first pixel is defined as the one with the lowest main and secondary
 *   coordinates
 * * each pixel is assigned:
 *     * a pixel coordinate identifier (`PixelCoordID_t`) of a pair of
 *       coordinate indices, sequence number respectively in the main and in
 *       the secondary direction, each starting from `0`; this index is
 *       represented with a signed type (`PixelCoordIndex_t`);
 *     * an index, absolute sequence number within the plane, so that most
 *       pixels consecutive in the main direction have consecutive indices;
 *       i.e., the absolute index is the secondary index multiplied by a stride
 *       (the number of assigned main indices) added to the main index 
 *       of the pixel; this results in the index being (weakly) sorted in the
 *       secondary direction; the index is represented with an _unsigned_ type;
 *     * a pixel coordinate set in the center of the pixel;
 * * the plane view is always defined as `geo::k3D`.
 * 
 * If there is no volume named after
 * `geo::PixelPlaneGeo::SensitivePixelPlaneVolumeName` in the plane, the plane
 * volume itself is used instead, with the same rules.
 * 
 * The plane defines two local reference frames, as in the standard
 * `geo::PlaneGeo` protocol, but these two match except for the origin.
 * 
 * The first, depending on pixel orientation, assumes the role of the
 * "wire base". It is defined by the normal to the plane (pointing toward the
 * center of the TPC), the main direction of the plane and the secondary one.
 * This is a positive orthogonal base.
 * The origin of this base is in the center of the pixel with index `{ 0, 0 }`.
 * Note that for this base to be correctly defined, the Geometry service has
 * to provide external information (for example, where the center of the
 * TPC is).
 *
 * The second frame plays the role of the "frame base" of `geo::PlaneGeo`.
 * It shares the directions with the wire base, with the following differences:
 * * the center is translated on the main and secondary directions to match the
 *     center of the sensitive plane center instead of the center of the first
 *     sensitive element;
 * * the main direction of the "wire" frame and the depth direction result
 *     antiparallel, i.e. the two bases are rotated by a right angle on the
 *     normal axis. This arises from the following prescriptions:
 *     * ( width x depth x normal ) and ( main x secondary x normal ) bases
 *         are defined positive
 *     * normal direction (and verse) in the two bases is the same
 *     * the width axis is the most aligned to the _z_ axis, therefore
 *         the width axis matches the secondary direction of the "wire" frame
 * 
 * Note that the center of this frame and the one of the "wire frame" also share
 * the coordinate along the drift direction.
 * In summary, the width direction is matching the secondary direction of the
 * wire base, including its verse; the normal is also shared; the depth
 * direction is defined according to the general prescription of this base,
 * resulting opposite to the main direction of the "wire frame".
 * 
 * As usual, all measures are in centimeters.
 *
 * @todo A protocol to pass the number of pixels in the two directions needs to
 *       be established. Options:
 *       * FHiCL configuration from the builder, passed to the object via
 *         constructor
 *       * custom properties of the GDML material the plane is made of;
 *         it's three quarters hack, but it has great reach; not proven to be
 *         feasible!
 *       * hard coded; wait.. whaaatt???
 * 
 */
class geo::PixelPlaneGeo: public geo::PlaneGeo {
  
  using DirIndex_t = std::size_t; ///< Type used for identifying a direction.
  
    public:
  
  /// The name of the solid representing the sensitive surface.
  static std::string const SensitivePixelPlaneVolumeName;
  
  // --- BEGIN -- Types for pixel identification -------------------------------
  /// @name Types for pixel identification
  /// @{
  
  /// Type of absolute index of the pixel in the plane.
  using PixelIndex_t = unsigned int;
  
  /// Type of each coordinate in the coordinate ID.
  using PixelCoordIndex_t = int;
  
  /// Type of coordinate identifier for a pixel on the plane.
  using PixelCoordID_t = geo::pixel::PixelCoordT<PixelCoordIndex_t>;
    
  /// Type of coordinates (in pixel units) for a point in the pixel plane.
  using PixelCoords_t = geo::pixel::PixelCoordT<double>;
  
  /// @}
  // --- END -- Types for pixel identification ---------------------------------
  
  
  /// A value for an invalid pixel index.
  static constexpr PixelIndex_t InvalidPixelIndex
    = std::numeric_limits<PixelIndex_t>::max();
  
  
  /// Constructor: a representation of a single pixel plane of the detector.
  PixelPlaneGeo(TGeoNode const& node, geo::TransformationMatrix&& trans);
  
  
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
   * * 1 _(default)_: also center and number of pixels
   * * 2: also information about resolution in the two directions
   * * 3: also information about width and depth
   * * 3: also information about normal directions (consistency cross check)
   * * 5: also coverage (relative to center)
   * * 6: also coverage (absolute coordinates)
   * * 7: also bounding box
   *
   */
  template <typename Stream>
  void PrintPixelPlaneInfo
    (Stream&& out, std::string indent = "", unsigned int verbosity = 1U)
    const;
  
  
  /**
   * @brief Prints information about the specified pixel.
   * @tparam Stream type of output stream to use
   * @param out stream to send the information to
   * @param coords coordinate ID of the pixel to be printed
   * @param indent prepend each line with this string
   * @param verbosity amount of information printed
   *
   * Note that the first line out the output is _not_ indented.
   *
   * Verbosity levels
   * -----------------
   *
   * * 0: only coordinate ID
   * * 1 _(default)_: also center
   * * 2: also start and end
   * * 3: also secondary angle (`ThetaZ()`)
   * 
   */
  template <typename Stream>
  void PrintPixelInfo(
    Stream&& out, PixelCoordID_t const coords,
    std::string indent = "", unsigned int verbosity = 1U
    ) const;
  
  
  /**
   * @brief Renders the information about the specified pixel into a string.
   * @param coords coordinate ID of the pixel to be printed
   * @param indent prepend each line with this string
   * @param verbosity amount of information printed
   * @return a string with all the information
   * @see `PrintPixelInfo()`
   *
   * See `PrintPixelInfo()` for the details of the arguments.
   */
  std::string PixelInfo(
    PixelCoordID_t const coords,
    std::string indent = "", unsigned int verbosity = 1U
    ) const;
  
  
    protected:
  
  /// Pixel initialization helper data structure.
  struct SquarePixelGeometry_t {
    
    using Vector2_t
      = ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<double>>;
    
    struct AxisInfo_t {
      /// Direction of this side in local (GDML) plane coordinates;
      /// modulus is the space covered by pixels [cm]
      LocalVector_t dir;
      double side; ///< Size of the side of each pixel (i.e. the pitch).
    };
    
    static constexpr std::size_t NSides = 2U;
    
    /// Information about one side of the pixels.
    std::array<AxisInfo_t, NSides> sides;
    
    LocalPoint_t center; ///< Center of the pixelized area.
    
  }; // struct SquarePixelGeometry_t
  
  
  // --- BEGIN -- Polymorphic implementation: anode plane --------------------
  /**
   * @name Polymorphic implementation: anode plane.
   */
  /// @{
  
  /// Returns the angle of the main direction from _z_ [rad]
  virtual double doThetaZ() const override;
  
  /// Returns the angle of the secondary direction from _z_ [rad]
  virtual double doPhiZ() const override;
  
  /// Returns the sine of the angle of the secondary direction from _z_.
  virtual double doSinPhiZ() const override;
  
  /// Returns the cosine of the angle of the secondary direction from _z_.
  virtual double doCosPhiZ() const override;
  
  /// Returns the total number of pixels.
  virtual unsigned int doNwires() const override;
  
  /**
   * @brief Returns a handle to the specified sensitive element on this plane.
   * @param iWire index of the sensitive element in [ `0`, `Nwires()` [
   * @return a `geo::WireGeo` handle for that wire
   * 
   * If the wire index is invalid, the result is undefined.
   * 
   * Technically, the handle refers back to this `geo::PixelPlaneGeo` object.
   */
  virtual geo::WireGeo doWire(unsigned int iWire) const override;
  
  /// Returns the pitch (pixel size) on the main direction.
  virtual double doWirePitch() const override;
  
  virtual bool doWireIDincreasesWithZ() const override;
  virtual geo::Vector_t doGetIncreasingWireDirection() const override;
  virtual geo::Vector_t doGetWireDirection() const override;
  virtual geo::WireID doNearestWireID
    (geo::Point_t const& pos) const override;
  virtual geo::WireGeo doNearestWire
    (geo::Point_t const& pos) const override;
  virtual geo::WireID doClosestWireID [[noreturn]]
    (geo::WireID::WireID_t wireNo) const override;
  /**
   * @brief Returns a "volume" enclosing the full pixel surface.
   * 
   * The returned volume has zero thickness and includes the surface of all
   * sensitive pixels.
   */
  virtual lar::util::simple_geo::Volume<> doCoverage() const override;
  /**
   * @brief Returns the relative coordinate in wire direction [cm]
   * @param point location in space
   * @param refWire wire to compute the coordinate from
   * @see `doWireCoordinate()`
   * 
   * Compared to `doWireCoordinate()`, this returns the result in distance
   * rather than sensitive element units.
   */
  virtual double doPlaneCoordinateFrom
    (geo::Point_t const& point, geo::WireGeo const& refWire) const override;
  /**
   * @brief Returns the relative coordinate in wire direction [cm]
   * @param point location in space
   * @see `doWireCoordinate()`
   * 
   * Compared to `doWireCoordinate()`, this returns the result in distance
   * rather than sensitive element units.
   */
  virtual double doPlaneCoordinate(geo::Point_t const& point) const override;
  /**
   * @brief Returns the relative coordinate in wire direction [cm]
   * @param point location in space
   * @see `doPlaneCoordinate()`
   * 
   * Compared to `doPlaneCoordinate()`, this returns the result in units of
   * sensitive element size, centered on the element number `0`.
   */
  virtual double doWireCoordinate(geo::Point_t const& point) const override;
  
  // the following methods directly refer to the decomposition object
  virtual WireDecomposedVector_t doDecomposePoint
    (geo::Point_t const& point) const override;
  virtual geo::Point_t doProjectionReferencePoint() const override;
  virtual WireCoordProjection_t doProjection
    (geo::Point_t const& point) const override;
  virtual WireCoordProjection_t doProjection
    (geo::Vector_t const& v) const override;
  virtual geo::Point_t doComposePoint
    (WireDecomposedVector_t const& decomp) const override;
  virtual geo::Point_t doComposePoint
    (double distance, WireCoordProjection_t const& proj) const override;
  
  /// Prints the custom pixel plane information.
  virtual std::string doPlaneInfo
    (std::string indent = "", unsigned int verbosity = 1U) const override;
  
  /// Does nothing since there are no elements to sort.
  virtual void doSortElements(geo::GeoObjectSorter const& sorter) override {}
  
  
  /**
   * @brief Completes the initialization after subelements are sorted.
   * 
   * Called by `geo::PlaneGeo::UpdateAfterSorting()` after setting the
   * new plane ID and geometry decomposition frame.
   * 
   * This implementation establishes the legacy orientation quantities.
   */
  virtual void doUpdateAfterSorting(geo::BoxBoundedGeo const& box) override;
  
  /// @}
  // --- END -- Polymorphic implementation: anode plane ----------------------
  
  
  // --- BEGIN -- Polymorphic implementation: wire abstraction ---------------
  /**
   * @name Polymorphic implementation: wire abstraction.
   * 
   * The plane will answer for all its "wires" (which are in fact rectangular
   * pixels).
   * 
   * The locator object hosts in `wireNo` the actual pixel index; the locator
   * can be formally transformed into a pixel index via `wireToPixelIndex()`.
   */
  /// @{
  
  /**
   * @brief Returns half of the longest of the two sides of the pixel.
   * @param wloc locator of the wire _(unused)_
   * @return half the lenth of the largest of the two pixel sides [cm]
   * 
   * This is a somehow arbitrary interpretation of the maximum radius of a wire.
   */
  virtual double doWireRMax(WireLocator const& wloc) const override;
  /**
   * @brief Returns half of the shortest of the two sides of the pixel.
   * @param wloc locator of the wire _(unused)_
   * @return half the lenth of the shortest of the two pixel sides [cm]
   * 
   * This is a quite arbitrary interpretation of the minimum radius of a wire.
   */
  virtual double doWireRMin(WireLocator const& wloc) const override;
  /**
   * @brief Returns the length of half the pitch on the main plane coordinate.
   * @param wloc locator of the wire _(unused)_
   * @return half the pitch on the main direction [cm]
   * 
   * The main coordinate on the pixel plane is equivalent to the wire direction
   * in a wire plane.
   */
  virtual double doWireHalfL(WireLocator const& wloc) const override;
  /**
   * @brief Fills an array with the coordinates of a pixel, shifted [cm]
   * @param wloc locator of the pixel
   * @param xyz pointer to the memory to be filled with the three coordinates
   * @param localz shift fraction along the "wire direction" (main coordinate)
   * @see `doWireGetPositionFromCenter()`
   * 
   * The position is stored in world coordinate, in centimeters.
   * The `xyz` location must have at least three `double` cells available to
   * store the result.
   * 
   * The shifting parameter `localz` describes how far from the center the
   * requested position is, along the main direction of the plane.
   * The unit of `localz` is the distance from the center to the border. It is
   * clamped to the range within the pixel (i.e. `-1` to `+1`).
   * Therefore, `doWireFillCenterXYZ(wloc, xyz, 1.0)` sets `xyz` on the center
   * of one border of the pixel, `doWireFillCenterXYZ(wloc, xyz, -1.0)` on the
   * opposite one (they are on the center of the border because the secondary
   * coordinate is not changed and still sits in the middle of its direction).
   * 
   * @deprecated Use `doWireGetPositionFromCenter()` instead.
   */
  virtual void doWireFillCenterXYZ
    (WireLocator const& wloc, double* xyz, double localz = 0.0) const override;
  /**
   * @brief Fills an array with the coordinates of a corner of a pixel [cm]
   * @param wloc locator of the pixel
   * @param xyz pointer to the memory to be filled with the three coordinates
   * @see `doWireGetStart()`
   * 
   * The position is stored in world coordinate, in centimeters.
   * The `xyz` location must have at least three `double` cells available to
   * store the result.
   * 
   * The "start" is defined as the corner with the lowest main and secondary
   * coordinates on the pixel.
   * 
   * @deprecated Use `doWireGetStart()` instead of `xyz` arrays.
   */
  virtual void doWireFillStartXYZ
    (WireLocator const& wloc, double* xyz) const override;
  /**
   * @brief Fills an array with the coordinates of a corner of a pixel [cm]
   * @param wloc locator of the pixel
   * @param xyz pointer to the memory to be filled with the three coordinates
   * @see `doWireGetEnd()`
   * 
   * The position is stored in world coordinate, in centimeters.
   * The `xyz` location must have at least three `double` cells available to
   * store the result.
   * 
   * The "end" is defined as the corner with the largest main and secondary
   * coordinates on the pixel.
   * 
   * @deprecated Use `doWireGetEnd()` instead of `xyz` arrays.
   */
  virtual void doWireFillEndXYZ
    (WireLocator const& wloc, double* xyz) const override;
  /**
   * @brief Returns the center of a pixel, shifted along main direction [cm]
   * @param wloc locator of the pixel
   * @param localz shift fraction along the "wire direction" (main coordinate)
   * @return the position of the shifted pixel center
   * 
   * The position is returned in world coordinate, in centimeters.
   * 
   * The shifting parameter `localz` describes how far from the center the
   * requested position is, along the main direction of the plane.
   * The unit of `localz` is the distance from the center to the border.
   * It is clamped to the range within the pixel (i.e. `-1` to `+1`).
   * Therefore, `doWireGetPositionFromCenter(wloc, 1.0)` returns the center
   * of one border of the pixel, `doWireGetPositionFromCenter(wloc, -1.0)` the
   * center of the opposite border (they are on the center of the border because
   * the secondary coordinate is not changed and still sits in the middle of its
   * direction).
   */
  virtual geo::Point_t doWireGetPositionFromCenter
    (WireLocator const& wloc, double localz) const override;
  /**
   * @brief Returns the center of a pixel, shifted along main direction [cm]
   * @param wloc locator of the pixel
   * @param localz shift fraction along the "wire direction" (main coordinate)
   * @return the position of the shifted pixel center
   * 
   * The position is returned in world coordinate, in centimeters.
   * 
   * The shifting parameter `localz` describes how far from the center the
   * requested position is, along the main direction of the plane.
   * The unit of `localz` is the distance from the center to the border.
   * Therefore, `doWireGetPositionFromCenter(wloc, 1.0)` returns the center
   * of one border of the pixel, `doWireGetPositionFromCenter(wloc, -1.0)` the
   * center of the opposite border (they are on the center of the border because
   * the secondary coordinate is not changed and still sits in the middle of its
   * direction).
   */
  virtual geo::Point_t doWireGetPositionFromCenterUnbounded
    (WireLocator const& wloc, double localz) const override;
  /// Returns the center of the specified pixel in world coordinates [cm]
  virtual geo::Point_t doWireGetCenter(WireLocator const& wloc) const override;
  /**
   * @brief Returns the position of a border of a pixel [cm]
   * @param wloc locator of the pixel
   * @return the position of a border of a pixel
   * 
   * The position is stored in world coordinate, in centimeters.
   * 
   * The "start" is defined as the corner with the lowest
   * main and secondary coordinates on the pixel.
   */
  virtual geo::Point_t doWireGetStart(WireLocator const& wloc) const override;
  /**
   * @brief Returns the position of a border of a pixel [cm]
   * @param wloc locator of the pixel
   * @return the position of a border of a pixel
   * 
   * The position is stored in world coordinate, in centimeters.
   * 
   * The "end" is defined as the corner with the highest
   * main and secondary coordinates on the pixel.
   */
  virtual geo::Point_t doWireGetEnd(WireLocator const& wloc) const override;
  /// Returns the size of the pixel on the main ("wire") direction.
  virtual double doWireLength(WireLocator const&) const override;
  /// Returns plane's `ThetaZ()`.
  virtual double doWireThetaZ(WireLocator const&) const override;
  /// Returns if the specified pixel main direction is aligned with another one.
  virtual bool doWireIsParallelTo
    (WireLocator const& wloc, geo::WireGeo const& wire) const override;
  /// Returns the main direction of the pixel.
  virtual geo::Vector_t doWireDirection(WireLocator const& wloc) const override;
  /**
   * @brief Returns the distance between pixels on the secondary direction.
   * @param wloc identifier of the local wire
   * @param wire reference pixel (or any wire object)
   * @return projection of distance from `wire` to `wloc` on secondary direction
   * 
   * The result may be negative if the pixel center is on the negative side of
   * the secondary direction with respect to the reference `wire` center.
   * 
   * Note that this is not a 2D distance, but instead the projection of that
   * distance on the secondary ("wire coordinate") direction.
   * This choice implements the `geo::WireGeo` protocol.
   */
  virtual double doWireDistanceFrom
    (WireLocator const& wloc, geo::WireGeo const& wire) const override;
  /// No geometry node is available for pixels.
  virtual TGeoNode const* doWireNode(WireLocator const&) const override
    { return nullptr; }
  
  /**
   * @brief Returns the width coordinate of the intersection of the pixel main
   *        axis with the width axis [cm]
   * @param wloc identifier of the pixel
   * @return the width coordinate
   * 
   * Liberal interpretation of the protocol: find the intersection of the wire
   * with a standard line with a constant coordinate in the frame reference,
   * and return the other coordinate (also from a standard reference) of the
   * point found.
   * 
   * That means that if we take a standard vertical wire plane with a inclined
   * wire, the plane has the _y_/_z_ plane frame reference, and _y_ = 0 can be
   * taken as the "standard" line. So we can intersect a wire with that
   * reference ("Y0") and get the other coordinate ("Z"), measuring it against
   * the "standard" _z_ = 0. If the wire is vertical, that is the world frame
   * _z_ as any point on the wire.
   * 
   * The generalization in sensitive plane terms is to get the frame coordinate
   * "width" of intersection of the prolongation of the pixel in its main
   * ("wire") direction with the width axis at depth 0.
   * 
   * In the end it does not really matter because this is used to relate
   * different wire planes in the same TPC, and there is only one pixel plane
   * in each TPC.
   */
  virtual double doWireComputeZatY0(WireLocator const& wloc) const override;
  
  /// Returns `nullptr`, as no transformation is available at pixel level.
  virtual WireLocalTransformation_t const* doWireTrans
    (WireLocator const&) const override
    { return nullptr; }
  
  virtual std::string doWireInfo(
    WireLocator const& wloc,
    std::string indent = "", unsigned int verbosity = 1
    ) const override;
  
  /// No update required. It does not even make sense here.
  virtual void doWireUpdateAfterSorting
    (WireLocator const& wloc, geo::WireID const&, bool) override
    {}
  
  /// @}
  // --- END -- Polymorphic implementation: wire abstraction -------------------
  
  
  // --- BEGIN -- Plane coordinate and ID conversions --------------------------
  /// @name Plane coordinate and ID conversions
  /// @{
  
  /**
   * @brief Converts a position on wire frame into pixel index.
   * @param point point projected into the "wire" reference frame [cm]
   * @return index of the pixel at the specified point
   * 
   * If the point projection is not on the plane, the result is undefined.
   */
  PixelIndex_t indexAt(WireCoordProjection_t const& point) const;
  
  /**
   * @brief Converts a position on wire frame into index coordinates.
   * @param point point projected into the "wire" reference frame [cm]
   * @return index coordinates of the pixel at the specified point
   * 
   * The index coordinates are not guaranteed to be on the plane.
   */
  PixelCoordID_t coordsOf(WireCoordProjection_t const& point) const;
  
  /**
   * @brief Converts a wire locator into pixel index coordinates.
   * @param wloc wire locator
   * @return pixel coordinate ID for `wloc`, undefined if index is invalid
   * 
   * The result is undefined if `wloc` is not valid (see `isPixelOnPlane()`).
   */
  PixelCoordID_t coordsOf(WireLocator const& wloc) const;
  
  /**
   * @brief Converts a pixel index into index coordinates.
   * @param index the pixel index
   * @return pixel coordinate ID for `index`, undefined if index is invalid
   * 
   * The result is undefined if `index` is not valid (see `isPixelOnPlane()`).
   */
  PixelCoordID_t coordsOf(PixelIndex_t const index) const;
  
  /**
   * @brief Returns the nearest coordinate index for a coordinate value.
   * @param coord the coordinate value on the "wire" frame [cm]
   * @param dir the direction of the coordinate
   * @return the nearest coordinate index (may be off plane)
   * 
   * The returned coordinate index is not guaranteed to be covered by the plane
   * (see `isPixelOnPlane()`).
   * 
   * The rounding is done so that the pixel _x_ covers the coordinate range
   * @f$ \left[ x - p/2, x + p/2 \right[ @f$ where _p_ is the pitch in the
   * coordinate direction.
   */
  PixelCoordIndex_t roundCoord(double coord, DirIndex_t const dir) const;
  
  /**
   * @brief Translates a coordinate into a the pixel coordinate.
   * @param coord the coordinate on "wire" frame to be translated [cm]
   * @param dir the direction of the coordinate
   * @return the coordinate index of `coord`
   * 
   * This method takes the coordinate on `dir` of the projection of a point
   * on the "wire" frame, and it converts to a pixel coordinate, that is a
   * value in units of pixel size.
   * 
   * The returned coordinate index is not guaranteed to be on the plane.
   */
  double pixelCoord(double coord, DirIndex_t const dir) const;
  
  /**
   * @brief Returns the nearest coordinate index for a pixel coordinate value.
   * @param coord the pixel coordinate value
   * @return the nearest coordinate index (may be off plane)
   * 
   * The returned coordinate index is not guaranteed to be covered by the plane
   * (see `isPixelOnPlane()`).
   * 
   * The rounding is done do that the pixel _x_ covers the pixel coordinate
   * range @f$ \left[ x - 0.5, x + 0.5 \right[ @f$.
   */
  PixelCoordIndex_t roundPixelCoord(double const coord) const;
  
  /**
   * @brief Converts index coordinates into an index.
   * @param coords index coordinates of the pixel
   * @return the index corresponding to `coords`, undefined if not on plane
   */
  PixelIndex_t indexOf(PixelCoordID_t const& coords) const;
  
  /// Converts a wire number into a pixel index.
  inline PixelIndex_t wireToPixelIndex(geo::WireID::WireID_t const wire) const;
  
  /// Converts a wire locator into a pixel index.
  inline PixelIndex_t wireToPixelIndex(WireLocator const& wloc) const;
  
  /// Converts a wire ID into a pixel index.
  inline PixelIndex_t wireToPixelIndex(geo::WireID const& wireid) const;
  
  /// Returns whether the projection `point` is on the sensitive plane.
  bool isOnPlane(WireCoordProjection_t const& point) const;
  
  /// Returns whether the coordinate `coord` is on the sensitive plane.
  bool isOnPlane(double const coord, DirIndex_t const dir) const;
  
  /// Returns whether the pixel at `coords` is on the sensitive plane.
  bool isPixelOnPlane(PixelCoords_t const& coords) const;
  
  /// Returns whether the index coordinates `coords` are on the sensitive plane.
  bool isPixelOnPlane(PixelCoordID_t const& coords) const;
  
  /// Returns whether the pixel coordinate `coord` is on the sensitive plane.
  bool isPixelOnPlane(double const coord, DirIndex_t const dir) const;
  
  /// Returns whether the coordinate index `coord` is on the sensitive plane.
  bool isPixelOnPlane(PixelCoordIndex_t const coord, DirIndex_t const dir) const;
  
  /// Returns whether the pixel index is on the sensitive plane.
  bool isPixelOnPlane(PixelIndex_t const index) const;
  
  /// Returns whether the wire ID is on the plane (plane ID is ignored).
  bool isWireIDvalid(geo::WireID const& wireid) const;
  
  /// Returns whether the specified wire number is on the plane.
  bool isWireIDvalid(geo::WireID::WireID_t const wire) const;
  
  
  // Returns whether the specified direction is supported.
  static constexpr bool isDirIndex(DirIndex_t dir)
    { return dir < geo::pixel::NCoords; }
  
  /// @}
  // --- END -- Plane coordinate and ID conversions ----------------------------
  
  
    private:
  
  /// Number of pixels along the main and secondary directions.
  std::array<unsigned int, geo::pixel::NCoords> fNPixels;
  
  /// Pixel pitch for the two directions [cm]
  /// @todo Must be read/deduced from somewhere
  std::array<double, geo::pixel::NCoords> fPitches;
  
  /// Decomposition on pixel coordinates (the "wire frame" of `geo::PlaneGeo`).
  WireDecomposer_t fDecompPixel;
  
  /// Displacement from the center to the center of the two positive sides.
  std::array<geo::Vector_t, geo::pixel::NCoords> fPixelDirs;
  
  double fThetaZ; ///< Cached value of @f$ \theta_{z} @f$ direction.
  double fPhiZ; ///< Cached value of @f$ \varphi_{z} @f$ direction.
  
  
  /// Return a temporary wire object handle.
  inline geo::WireGeo getWire(WireLocator const& wloc) const;
  
  /// Return a temporary wire object handle.
  geo::WireGeo getWire(geo::WireID::WireID_t const iWire) const;
  
  
  // --- BEGIN --- Candidate extensions to `geo::PlaneGeo` interface -----------
  // the following methods do not override anything --- yet.
  
  /**
   * @brief  Returns the direction of the sensitive elements along the specified
   *         plane direction
   * @param dir index of the plane direction to query
   * @return a unit vector pointing to the the requested direction
   * 
   * The return value is undefined if the direction is not supported.
   */
  virtual geo::Vector_t doSensElemDir(DirIndex_t const dir) const;
  
  /// Returns the main direction of sensitive elements on the plane.
  virtual geo::Vector_t doSensElemMainDir() const;
  
  /// Returns the main direction of sensitive elements on the plane.
  /// direction [cm]
  virtual geo::Vector_t doSensElemSecondaryDir() const;
  
  
  /**
   * @brief  Returns the length of the sensitive plane along the specified
   *         direction [cm]
   * @param dir index of the plane direction to query
   * @return a unit vector pointing to the the requested direction
   * 
   * The return value is undefined if the direction is not supported.
   */
  virtual double doSensElemDirSize(DirIndex_t const dir) const;
  
  /// Returns the length of the sensitive plane along the main direction [cm]
  virtual double doSensElemMainDirSize() const;
  
  /// Returns the length of the sensitive plane along the secondary direction
  /// [cm]
  virtual double doSensElemSecondaryDirSize() const;
  
  
  /**
   * @brief  Returns the number of sensitive elements along the specified plane
   *         direction.
   * @param dir index of the plane direction to query
   * @return number of sensitive elements along the direction `dir`
   * 
   * The return value is undefined if the direction is not supported.
   */
  virtual unsigned int doNsensElem(DirIndex_t const dir) const;
  
  /// Returns the number of sensitive elements along the main plane direction.
  virtual unsigned int doNsensElemMain() const;
  
  /// Returns the number of sensitive elements along the secondary plane
  /// direction.
  // This should return `1` for wire planes.
  virtual unsigned int doNsensElemSecondary() const;
  
  
  /**
   * @brief  Returns the pitch between sensitive elements along the specified
   *         plane direction [cm]
   * @param dir index of the plane direction to query
   * @return pitch between sensitive elements along the direction `dir`
   * 
   * The return value is undefined if the direction is not supported.
   */
  virtual double doSensElemPitch(DirIndex_t const dir) const;
  
  /// Returns the pitch of sensitive elements along the main plane direction [cm]
  virtual double doSensElemPitchMain() const;
  
  /// Returns the pitch of sensitive elements along the secondary plane
  /// direction [cm]
  // This should return `Depth()` for wire planes.
  virtual double doSensElemPitchSecondary() const;
  
  
  // --- END --- Candidate extensions to `geo::PlaneGeo` interface -------------
  
  
  // --- BEGIN -- Implementation of the candidate interface extensions ---------
  /**
   * @brief Returns the specified direction of sensitive elements.
   * @param dir index of the plane direction to query
   * @return a unit vector pointing to the requested direction
   * 
   * The return value is a null vector if the direction is not supported.
   */
  // this is to prove that we don't have to support `TVector3` forever...
  geo::Vector_t getSensElemDir(DirIndex_t const dir) const;
  
  
  /**
   * @brief Returns size of the sensitive plane on the specified direction [cm]
   * @param dir index of the plane direction to query
   * @return the size of the sensitive plane along the specified direction [cm]
   * 
   * The return value is `0` if the direction is not supported.
   */
  double getSensElemDirSize(DirIndex_t const dir) const;
  
  
  /**
   * @brief Returns half a pixel step in the specified direction.
   * @param dir index of the plane direction to query
   * @return vector from the center of a pixel to ts border following `dir`
   * 
   * The return value is undefined if the direction is not supported.
   */
  geo::Vector_t getSensElemHalfStepDir(DirIndex_t const dir) const;
  
  
  /**
   * @brief Returns the number of sensitive elements along the specified plane
   *        direction
   * @param dir index of the plane direction to query
   * @return number of sensitive elements along the direction `dir`
   * 
   * The return value is undefined if the direction is not supported.
   */
  unsigned int getNsensElem(DirIndex_t const dir) const;
  
  
  /// Returns the total number of sensitive elements on the plane.
  unsigned int getNsensElem() const;
  
  
  /**
   * @brief Returns the pitch of sensitive elements along the specified plane
   *        direction.
   * @param dir index of the plane direction to query
   * @return pitch between sensitive elements along the direction `dir`
   * 
   * The return value is undefined if the direction is not supported.
   */
  double getSensElemPitch(DirIndex_t const dir) const;
  
  
  /**
   * @brief Returns the coordinate `dir` of `point` with respect to a reference.
   * @param point the location to learn the coordinate of, in world frame [cm]
   * @param ref the sensitive element whose center is used as reference
   * @param dir index of the plane direction to query
   * @return pitch between sensitive elements along the direction `dir`
   * 
   * The returned coordinate may be out of the plane.
   */
  double getPlaneCoordinateFrom
    (geo::Point_t const& point, geo::WireGeo const& ref, DirIndex_t const dir)
    const;
  
  /**
   * @brief Returns the coordinates of a point on the two directions [cm]
   * @param point the location to learn the coordinates of, in world frame [cm]
   * @return coordinates of a point on the two directions
   * @see `getPlaneCoordinate()`
   * 
   * The returned vector contains a `main()` and a `secondary()` component, both
   * measured in centimeters, which represent the displacement of the point
   * projected on the plane from the center of the first pixel.
   * The returned object is unstable, and it will easily decay into a normal
   * 2D vector (`WireCoordProjection_t`) when transformed or operated upon.
   * 
   * The return value is undefined if the direction is not supported.
   */
  WireCoordProjection_t getPlaneCoordinates(geo::Point_t const& point) const;
  
  /**
   * @brief Returns the coordinate `dir` of `point` [cm]
   * @param point the location to learn the coordinate of, in world frame [cm]
   * @param dir index of the plane coordinate to query
   * @return the coordinate of `point` in the specified coordinate
   * 
   * The returned coordinate may be out of the plane.
   * If `dir` is invalid, the behavior is undefined.
   */
  double getPlaneCoordinate
    (geo::Point_t const& point, DirIndex_t const dir) const;
  
  /// Returns the main coordinate of the specified point (in "wire" frame) [cm]
  double getMainPlaneCoordinate(geo::Point_t const& point) const;
  
  /// Returns secondary coordinate of the specified point (in "wire" frame) [cm]
  double getSecPlaneCoordinate(geo::Point_t const& point) const;
  
  /**
   * @brief Returns the value of the specified point on `dir` coordinate.
   * @param point the point to get the coordinate of
   * @param dir the direction of the coordinate to be extracted
   * @return the coordinate of `point` on `dir` direction, in element units
   * 
   * The returned value is `0.0` for a point exactly at the center of the first
   * sensitive element, `0.5` for a point on the border between the first two
   * sensitive elements, `1.0` for one at the center of the second sensitive
   * element, etc.
   * The coordinate may be outside of the sensitive area of the plane.
   */
  double getWireCoordinate
    (geo::Point_t const& point, DirIndex_t const dir) const;
  
  /// Returns the center of the specified pixel in world coordinates [cm]
  geo::Point_t getSensElemCenter(PixelCoordID_t const& coords) const;
  
  
  /// Returns the "half length" (i.e. in main direction) of the pixel [cm]
  double getPixelHalfL(WireLocator const&) const;
  
  // --- END -- Implementation of the candidate interface extensions -----------
  
  
  // --- BEGIN --- Initialization procedures -----------------------------------
  
  /// Extracts in some way the pixel information from GDML.
  virtual SquarePixelGeometry_t completePixelGeometry
    (SquarePixelGeometry_t const& info) const;
  
  /// Utilizes the `pixelGeometry` information for "wire" frame initialization.
  void initializePixelGeometry(SquarePixelGeometry_t const& pixelGeometry);
  
  /// Updates the number of pixels and the pitches.
  void UpdateNpixelsAndPitches();
  
  /// Updates the pixel ("wire") decomposition frame.
  void UpdateDecompPixel();
  
  /// Updates the cached pixel directions (half steps).
  void UpdatePixelDirs();
  
  /// Aligns the origin of the pixel decomposition frame with the inner face.
  void UpdateDecompPixelOrigin();
  
  /// Shifts the formal pixel plane center to the plane of sensitive elements.
  void UpdatePlaneCenter();
  
  /// Updates the cached angles related to the pixel grid orientation.
  void UpdateAngles();
  
  /// Updates the internally used active area.
  void UpdateActiveArea();

  // --- END --- Initialization procedures -------------------------------------
  
  
  /// Returns world position of the center of the first sensitive element [cm]
  geo::Point_t firstPixelCenter() const;
  
  /// Applies a transformation from the specified plane center to pixel #0.
  geo::Point_t fromCenterToFirstPixel
    (geo::Point_t const& pixelPlaneCenter) const;
  
  /// Extracts the information about the size of pixels.
  void discoverPitches();

  /// Returns the plane box thickness (not width, not depth... the other one!).
  double extractPlaneThickness() const;
  
  
  /// Returns the name of the specified direction.
  static std::string getDirectionName(DirIndex_t const dir);
  
}; // class PixelPlaneGeo


//------------------------------------------------------------------------------
//--- inline implementation
// -----------------------------------------------------------------------------
inline auto geo::PixelPlaneGeo::wireToPixelIndex
  (geo::WireID const& wireid) const -> PixelIndex_t
{
  return wireToPixelIndex(wireid.Wire);
} // geo::PixelPlaneGeo::wireToPixelIndex(WireID)


// -----------------------------------------------------------------------------
inline auto geo::PixelPlaneGeo::wireToPixelIndex(WireLocator const& wloc) const
  -> PixelIndex_t
{
  return wireToPixelIndex(wloc.wireNo);
} // geo::PixelPlaneGeo::wireToPixelIndex(WireLocator)


// -----------------------------------------------------------------------------
inline auto geo::PixelPlaneGeo::wireToPixelIndex
  (geo::WireID::WireID_t const wire) const -> PixelIndex_t
{
  PixelIndex_t const index = wire; // just the same value
  assert(isPixelOnPlane(index));
  return index;
} // geo::PixelPlaneGeo::wireToPixelIndex()


//------------------------------------------------------------------------------
inline geo::WireGeo geo::PixelPlaneGeo::getWire(WireLocator const& wloc) const {
  return getWire(wireToPixelIndex(wloc));
} // geo::PixelPlaneGeo::getWire()


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void geo::PixelPlaneGeo::PrintPixelPlaneInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1U */
) const {
  
  using namespace geo::pixel;
  
  //----------------------------------------------------------------------------
  out << "plane " << std::string(ID());
  
  if (verbosity-- <= 0) return; // 0
  
  //----------------------------------------------------------------------------
  out
    << " at " << GetCenter<geo::Vector_t>() << " cm"
    << " with " << getNsensElem(ixMain) << "x" << getNsensElem(ixSec)
    << " pixels";
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  
  out << "\n" << indent
    << "plane orientation " << OrientationName(Orientation())
    << ", measuring " << ViewName(View())
    ;
  
  for (DirIndex_t dir: { ixMain, ixSec }) {
    out << "\n" << indent
      << "- " << getDirectionName(dir) << " direction: "
      << getSensElemDir(dir) << " with " << getNsensElem(dir) << " pixel "
      << getSensElemPitch(dir) << " mm each"
      ;
  } // for dir
  
  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  
  auto const& widthDir = WidthDir<geo::Vector_t>();
  auto const& depthDir = DepthDir<geo::Vector_t>();
  auto const& frameNormalDir = fDecompFrame.NormalDir();
  
  out << "\n" << indent
    << "width " << Width() << " cm in direction: " << widthDir
    << ", depth " << Depth() << " cm in direction: " << depthDir
    << " [normal: " << frameNormalDir << "]"
    ;
    
  if (verbosity-- <= 0) return; // 3
  
  //----------------------------------------------------------------------------
  auto const& normal = GetNormalDirection<geo::Vector_t>();
  auto const& wireNormalDir = fDecompPixel.NormalDir();
  out << "\n" << indent
    << "normal to plane: " << normal
    << " [pixel frame normal: " << wireNormalDir << "]"
    ;
  
  if (verbosity-- <= 0) return; // 4
  
  //----------------------------------------------------------------------------
  // get the area spanned by the pixels in local coordinates
  out << "\n" << indent << "pixel cover width "
    << ActiveArea().width.lower << " to " << ActiveArea().width.upper
    << ", depth "
    << ActiveArea().depth.lower << " to " << ActiveArea().depth.upper
    << " cm";
  if (verbosity-- <= 0) return; // 5
  
  //----------------------------------------------------------------------------
  // get the area spanned by the pixels in local coordinates
  auto const coverage = Coverage();
  out << "\n" << indent << "  absolute coordinates: ("
    << coverage.Min().x << ", " << coverage.Min().y << ", " << coverage.Min().z
    << ") to ("
    << coverage.Max().x << ", " << coverage.Max().y << ", " << coverage.Max().z
    << ") cm";
  if (verbosity-- <= 0) return; // 6
  
  //----------------------------------------------------------------------------
  // print also the containing box
  auto const box = BoundingBox();
  out << "\n" << indent
    << "bounding box: " << box.Min() << " -- " << box.Max();
  
//  if (verbosity-- <= 0) return; // 7
  
  //----------------------------------------------------------------------------
  
} // geo::PixelPlaneGeo::PrintPixelPlaneInfo()


//------------------------------------------------------------------------------
template <typename Stream>
void geo::PixelPlaneGeo::PrintPixelInfo(
  Stream&& out, PixelCoordID_t const coords,
  std::string indent /* = "" */, unsigned int verbosity /* = 1U */
  ) const
{

  //----------------------------------------------------------------------------
  out << "pixel ( " << coords.main() << " ; " << coords.secondary() << " )";

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  geo::WireGeo const wire = getWire({ indexOf(coords) });
  
  out << " at " << wire.GetCenter<geo::Point_t>();

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------
  out << "; start: " << wire.GetStart<geo::Point_t>() << " cm, end: "
    << wire.GetEnd<geo::Point_t>() << " cm";

  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  out << "; theta(z)=" << wire.ThetaZ() << " rad";

//  if (verbosity-- <= 0) return; // 3
  
} // geo::PixelPlaneGeo::PrintPixelInfo()


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEO_H

