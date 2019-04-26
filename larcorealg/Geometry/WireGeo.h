/**
 * @file    larcorealg/Geometry/WireGeo.h
 * @brief   Interface for a active readout element on a TPC plane.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    April 24, 2019
 * @see     larcorealg/Geometry/PlaneGeo.h, larcorealg/Geometry/WireGeo.cxx
 * @ingroup Geometry
 * 
 * Whaaaat?!? where is everything??
 * 
 * `geo::WireGeo` was suppressed and absorbed by `geo::PlaneGeo`, Inc.
 * A fork of `geo::WireGeo` is `geo::SenseWireGeo`, which sheltered and hid
 * under the protective wing of `geo::WirePlaneGeo`.
 * 
 */

#ifndef LARCOREALG_GEOMETRY_WIREGEO_H
#define LARCOREALG_GEOMETRY_WIREGEO_H

// LArSoft
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/GeoElementTraits.h" // geo::element_traits
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util ns

// C/C++ standard libraries
#include <string>
#include <utility> // std::forward()
#include <limits>
#include <type_traits> // std::is_aggregate_v
#include <cmath> // std::abs(), ...

// ROOT
#include "TVector3.h"


// forward declarations
class TGeoNode;


namespace geo {
  
  namespace details {
    
    //--------------------------------------------------------------------------
    /// Information needed to locate a wire on the plane. Simple so far.
    struct WireLocator {
    
      /// Number of the wire on the plane.
      geo::WireID::WireID_t wireNo
        = std::numeric_limits<geo::WireID::WireID_t>::max();
      
    }; // struct WireLocator
    
    // not necessary for it to be POD, but if not POD any more we want to know:
    static_assert(std::is_aggregate_v<WireLocator>);
    
    
    //--------------------------------------------------------------------------
    
  } // namespace details
  
  
  //----------------------------------------------------------------------------
  // forward declaration;
  // this `WireGeo` is nothing without its plane; it needs to know it
  // and be liked and befriended by it
  class PlaneGeo;
  
  class WirePtr; // defined later
  
  //----------------------------------------------------------------------------
  /// A wrapper collecting all information of a wire and making it look like
  /// an object.
  class WireGeo {
    
    using WireLocator = details::WireLocator;
    
    template <typename Ret, typename... Args>
    using ConstPlaneFunc_t
      = Ret (geo::PlaneGeo::*)(WireLocator const&, Args...) const;
    
    template<typename Ret, typename... FuncArgs, typename... Args>
    decltype(auto) askThePlane
      (ConstPlaneFunc_t<Ret, FuncArgs...> func, Args&&... args) const
      { return (plane().*func)(location(), std::forward<Args>(args)...); }
    
    
      public:
    
    using DefaultVector_t = TVector3; ///< :-(
    using DefaultPoint_t = TVector3; ///< :-(
    
    // --- BEGIN -- Types for geometry-local reference vectors ---------------
    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     *
     * These types represents points and displacement vectors in the reference
     * frame defined in the wire geometry "box" from the GDML geometry
     * description.
     *
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     *
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::WireGeo` have the same type but are not compatible.
     */

    /// Tag for vectors in the "local" GDML coordinate frame of the plane.
    struct WireGeoCoordinatesTag {};

    /// Type of points in the local GDML wire plane frame.
    using LocalPoint_t = geo::Point3DBase_t<WireGeoCoordinatesTag>;

    /// Type of displacement vectors in the local GDML wire plane frame.
    using LocalVector_t = geo::Vector3DBase_t<WireGeoCoordinatesTag>;

    using LocalTransformation_t = geo::LocalTransformationGeo
        <ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;
    
    ///@}
    // --- END -- Types for geometry-local reference vectors -----------------
    
    /**
     * @brief Constructor from the reference plane and the wire index.
     * @param plane the plane this object belongs to
     * @param wireNo the number of this wire within that plane
     * 
     * The object `plane` is responsible of providing all the necessary
     * information to implement the wire functionality.
     */
    WireGeo
      (geo::PlaneGeo const& plane, geo::WireID::WireID_t wireNo) noexcept;

    /**
     * @brief Constructor from the reference plane and the wire index.
     * @param plane the plane this object belongs to
     * @param wireid the ID of this wire within that plane
     * @throw cet::exception (category: `"Geometry"`) if `wireid` is on a
     *                       different `plane`
     */
    WireGeo(geo::PlaneGeo const& plane, geo::WireID const& wireid);
    
    
    // --- BEGIN -- Introspection ----------------------------------------------
    /// @name Introspection
    /// @{
    
    using ID_t = geo::WireID; ///< Type of identifier of this element.
    
    /// Returns the complete ID of this wire.
    geo::WireID ID() const;
    
    /// Returns whether this wire is valid (includes a check on the plane).
    bool IsValid() const;
    
    /// @}
    // --- END -- Introspection ------------------------------------------------
    
    // --- BEGIN -- Size and coordinates -------------------------------------
    /// @{
    /// @name Size and coordinates

    /// Returns the outer half-size of the wire [cm]
    double RMax() const;

    /// Returns the inner radius of the wire (usually 0) [cm]
    double RMin() const;

    /// Returns half the length of the wire [cm]
    double HalfL() const;

    /**
     * @brief Fills the world coordinate of a point on the wire
     * @param xyz _(output)_ the position to be filled, as [ x, y, z ] (in cm)
     * @param localz distance of the requested point from middle of the wire
     * @see `GetCenter()`, `GetStart()`, `GetEnd()`, `GetPositionFromCenter()`
     *
     * The center of the wires corresponds to `localz` equal to `0`; negative
     * positions head toward the start of the wire, positive toward the end.
     *
     * If the `localz` position would put the point outside the wire, the
     * returned position is the wire end closest to the requested position.
     *
     * @deprecated Use the version returning a vector instead.
     */
    void GetCenter(double* xyz, double localz = 0.0) const;

    /// Fills the world coordinate of one end of the wire
    /// @deprecated Use the version returning a vector instead.
    void GetStart(double* xyz) const;

    /// Fills the world coordinate of one end of the wire
    /// @deprecated Use the version returning a vector instead.
    void GetEnd(double* xyz) const;
    
    
    //@{
    /**
     * @brief Returns the position (world coordinate) of a point on the wire
     * @tparam Point type of vector to be returned
     * @param localz distance of the requested point from the middle of wire
     * @return the position of the requested point (in cm)
     * @see `GetCenter()`, `GetStart()`, `GetEnd()`,
     *      `GetPositionFromCenterUnbounded()`
     *
     * The center of the wires corresponds to `localz` equal to `0`; negative
     * positions head toward the start of the wire, positive toward the end.
     *
     * If the `localz` position would put the point outside the wire, the
     * returned position is the wire end closest to the requested position.
     */
    template <typename Point>
    Point GetPositionFromCenter(double localz) const
      { return geo::vect::convertTo<Point>(GetPositionFromCenterImpl(localz)); }
    DefaultPoint_t GetPositionFromCenter(double localz) const
      { return GetPositionFromCenter<DefaultPoint_t>(localz); }
    //@}

    //@{
    /**
     * @brief Returns the position (world coordinate) of a point on the wire
     * @tparam Point type of vector to be returned
     * @param localz distance of the requested point from the middle of wire
     * @return the position of the requested point (in cm)
     * @see `GetCenter()`, `GetStart()`, `GetEnd()`, `GetPositionFromCenter()`
     *
     * The center of the wires corresponds to `localz` equal to `0`; negative
     * positions head toward the start of the wire, positive toward the end.
     *
     * If the `localz` position would put the point outside the wire, the
     * returned position will lie beyond the end of the wire.
     */
    template <typename Point>
    Point GetPositionFromCenterUnbounded(double localz) const
      {
        return geo::vect::convertTo<Point>
          (GetPositionFromCenterUnboundedImpl(localz)); 
      }
    DefaultPoint_t GetPositionFromCenterUnbounded(double localz) const
      { return GetPositionFromCenterUnbounded<DefaultPoint_t>(localz); }
    //@}

    //@{
    /// Returns the world coordinate of the center of the wire [cm]
    /// @tparam Point type of the point being returned
    template <typename Point>
    Point GetCenter() const
      { return geo::vect::convertTo<Point>(GetCenterImpl()); }
    DefaultPoint_t GetCenter() const { return GetCenter<DefaultPoint_t>(); }
    //@}

    //@{
    /// Returns the world coordinate of one end of the wire [cm]
    /// @tparam Point type of the point being returned
    template <typename Point>
    Point GetStart() const
      { return geo::vect::convertTo<Point>(GetStartImpl()); }
    DefaultPoint_t GetStart() const { return GetStart<DefaultPoint_t>(); }
    //@}

    //@{
    /// Returns the world coordinate of one end of the wire [cm]
    /// @tparam Point type of the point being returned
    template <typename Point>
    Point GetEnd() const
      { return geo::vect::convertTo<Point>(GetEndImpl()); }
    DefaultPoint_t GetEnd() const { return GetEnd<DefaultPoint_t>(); }
    //@}

    /// Returns the wire length in centimeters
    double Length() const;

    /// Returns the z coordinate, in centimetres, at the point where y = 0.
    /// Assumes the wire orthogonal to x axis and the wire not parallel to z.
    double ComputeZatY0() const;

    /**
     * @brief Returns 3D distance from the specified wire
     * @return the signed distance in centimetres (0 if wires not parallel)
     *
     * If the specified wire is "ahead" in z respect to this, the distance is
     * returned negative.
     */
    double DistanceFrom(WireGeo const& wire) const;
    
    /// @}
    // --- END -- Size and coordinates ---------------------------------------
    
    
    // --- BEGIN -- Orientation and angles -----------------------------------
    /// @{
    /// @name Orientation and angles

    /// Returns angle of wire w.r.t. _z_ axis in the _y_ - _z_ plane [radians]
    double ThetaZ() const;

    /**
     * Returns angle of wire with respect to _z_ axis.
     * @param degrees return the angle in degrees rather than radians
     * @return wire angle
     */
    double ThetaZ(bool degrees) const
      { return degrees? util::RadiansToDegrees(ThetaZ()): ThetaZ(); }
    
    
    //@{
    /// Returns trigonometric operations on `ThetaZ()`.
    double CosThetaZ() const;
    double SinThetaZ() const;
    double TanThetaZ() const;
    //@}

    /// Returns if this wire is horizontal (@f$ \theta_{z} \approx 0 @f$).
    /// @deprecated The method name is confusing, and it was unused as of `v08_16_00` anyway.
    [[deprecated]] bool isHorizontal() const { return std::abs(SinThetaZ()) < 1e-5; }

    /// Returns if this wire is vertical (@f$ \theta_{z} \approx \pi/2 @f$).
    /// @deprecated The method name is confusing, and it was unused as of `v08_16_00` anyway.
    [[deprecated]] bool isVertical() const { return std::abs(CosThetaZ()) < 1e-5; }

    /// Returns if this wire is parallel to another.
    bool isParallelTo(WireGeo const& wire) const;

    //@{
    /// Returns the wire direction as a norm-one vector.
    /// @tparam Vector type of the vector being returned
    template <typename Vector>
    Vector Direction() const
      { return geo::vect::convertTo<Vector>(DirectionImpl()); }
    DefaultVector_t Direction() const { return Direction<DefaultVector_t>(); }
    //@}
    /// @}
    // --- END -- Orientation and angles -------------------------------------
    
    // --- BEGIN -- Coordinate conversion ------------------------------------
    /// @{
    /**
     * @name Coordinate conversion
     *
     * Local points and displacement vectors are described by the types
     * `WireGeo::LocalPoint_t` and `WireGeo::LocalVector_t`,
     * respectively.
     */

    /// Transform point from local wire frame to world frame.
    void LocalToWorld(const double* wire, double* world) const
      { trans()->LocalToWorld(wire, world); }

    /// Transform point from local wire frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
      { return trans()->toWorldCoords(local); }

    /// Transform direction vector from local to world.
    void LocalToWorldVect(const double* wire, double* world) const
      { trans()->LocalToWorldVect(wire, world); }

    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
      { return trans()->toWorldCoords(local); }

    /// Transform point from world frame to local wire frame.
    void WorldToLocal(const double* world, double* wire) const
      { trans()->WorldToLocal(world, wire); }

    /// Transform point from world frame to local wire frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
      { return trans()->toLocalCoords(world); }

    /// Transform direction vector from world to local.
    void WorldToLocalVect(const double* world, double* wire) const
      { trans()->WorldToLocalVect(world, wire); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
      { return trans()->toLocalCoords(world); }

    /// @}
    // --- END -- Coordinate conversion --------------------------------------


    /**
     * @brief Returns the geometry node representing this wire.
     * @return a pointer to the `TGeoNode`, or `nullptr` if not available.
     * 
     * If the geometry description does not explicitly include this element
     * (for example because wires are created by procedure), the returned
     * value is `nullptr`.
     */
    TGeoNode const* Node() const;

    
    // --- BEGIN -- Printout of wire information -----------------------------
    /// @name Printout of wire information
    /// @{
    /**
     * @brief Prints information about this wire.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
     * @param indent prepend each line with this string
     * @param verbosity amount of information printed
     * @see `WireInfo()`
     *
     * All arguments are equivalent to the ones of `WireInfo()`.
     */
    template <typename Stream>
    void PrintWireInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1)
      const
      { out << WireInfo(indent, verbosity); }

    /**
     * @brief Returns a string with all the information of the wire.
     * @see `PrintWireInfo()`
     *
     * Note that the first line out the output is _not_ indented.
     */
    std::string WireInfo
      (std::string indent = "", unsigned int verbosity = 1) const;
    
    
    /**
     * @brief Prints information about this wire.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
     * @param indent prepend each line with this string
     * @param verbosity amount of information printed
     * @see `PrintWireInfo()`, `GenericWireInfo()`
     *
     * The information printed by this method follows a standard format for
     * the sensitive elements of all `geo::PlaneGeo` classes.
     * This is currently driven by the assumption that the element is actually
     * a wire.
     * For a specific printout, use `PrintWireInfo()` instead.
     * 
     * Note that the first line out the output is _not_ indented.
     *
     * Verbosity levels
     * -----------------
     *
     * * 0: only start and end
     * * 1 _(default)_: also length
     * * 2: also angle with z axis
     * * 3: also center
     * * 4: also direction
     *
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintGenericWireInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;
    
    
    /**
     * @brief Returns a string with all the information of the wire.
     * @see `PrintWireInfo()`
     *
     * Note that the first line out the output is _not_ indented.
     * All arguments are equivalent to the ones of `PrintGenericWireInfo()`.
     */
    std::string GenericWireInfo
      (std::string indent = "", unsigned int verbosity = 1) const;
    
    
    /// Maximum verbosity supported by `PrintWireInfo()`.
    static constexpr unsigned int MaxVerbosity = 4;
    
    /// @}
    // --- END -- Printout of wire information -------------------------------
    
    
    /// Returns a ("smart") pointer to this wire.
    // I expect this will backfire sooner or later.
    constexpr geo::WirePtr operator& () const;
    
    
    /// Returns the pitch (distance on y/z plane) between two wires, in cm
    static double WirePitch(WireGeo const& w1, WireGeo const& w2)
      { return std::abs(w2.DistanceFrom(w1)); }
    
    
    /// Returns an invalid wire.
    static constexpr WireGeo InvalidWire() { return {}; }
    
    
      private:
    
    struct WireRefData {
      /// Pointer to the containing plane. Never null.
      geo::PlaneGeo const* plane = nullptr;
      
      /// Information for the plane to locate this wire.
      WireLocator wireLoc;
    }; // struct WireRefData
    
    static_assert(std::is_aggregate_v<WireRefData>);
    
    
    WireRefData fData; ///< All data membership club of this object.
    
    
    /// Private creation of an invalid wire. Exposed via `InvalidWire()`.
    constexpr WireGeo() noexcept = default;
    
    
    geo::PlaneGeo const& plane() const { return *(fData.plane); }
    WireLocator const& location() const { return fData.wireLoc; }
    
    
    LocalTransformation_t const* trans() const;
    
    
    /// `GetPositionFromCenter()` returning a `geo::Point_t`.
    geo::Point_t GetPositionFromCenterImpl(double localz) const;
    
    /// `GetPositionFromCenterUnbounded()` returning a `geo::Point_t`.
    geo::Point_t GetPositionFromCenterUnboundedImpl(double localz) const;
    
    /// `GetCenter()` returning a `geo::Point_t`.
    geo::Point_t GetCenterImpl() const;
    
    /// `GetStart()` returning a `geo::Point_t`.
    geo::Point_t GetStartImpl() const;
    
    /// `GetEnd()` returning a `geo::Point_t`.
    geo::Point_t GetEndImpl() const;
    
    /// `Direction()` returning a `geo::Vector_t`.
    geo::Vector_t DirectionImpl() const;
    
    
  }; // class WireGeo
  
  static_assert(std::is_copy_constructible_v<WireGeo>);
  static_assert(std::is_move_constructible_v<WireGeo>);
  static_assert(std::is_copy_assignable_v   <WireGeo>);
  static_assert(std::is_move_assignable_v   <WireGeo>);
  
  //----------------------------------------------------------------------------
  
  /// A replacement for `geo::WireGeo const*`.
  /// Smart pointer, if you have a low bar for "smart".
  class WirePtr {
    
    /// Copy of the wire proxy.
    geo::WireGeo fWire = geo::WireGeo::InvalidWire();
    
      public:
    
    using element_type = geo::WireGeo const;
    
    /// Creates an invalid wire pointer.
    constexpr WirePtr() = default;
    
    constexpr WirePtr(geo::WireGeo const* from) noexcept: fWire(*from) {}
    
    
    geo::WireGeo const& operator*() const { return fWire; }
    
    geo::WireGeo const* operator->() const { return std::addressof(fWire); }
    
    
    /// Returns whether the pointer wire is valid.
    operator bool() const { return fWire.IsValid(); }
    
    /// Returns whether the pointer wire is invalid.
    bool operator! () const { return !(fWire.IsValid()); }
    
    static constexpr WirePtr pointer_to(geo::WireGeo const& from) noexcept
      { return &from; }
    
  }; // class WirePtr
  
  static_assert(std::is_copy_constructible_v<WirePtr>);
  static_assert(std::is_move_constructible_v<WirePtr>);
  static_assert(std::is_copy_assignable_v   <WirePtr>);
  static_assert(std::is_move_assignable_v   <WirePtr>);
  
  /// `WirePtr` *is* a constant pointer already.
  using WireConstPtr = WirePtr;
  
  static_assert(std::is_copy_constructible_v<WireConstPtr>);
  static_assert(std::is_move_constructible_v<WireConstPtr>);
  static_assert(std::is_copy_assignable_v   <WireConstPtr>);
  static_assert(std::is_move_assignable_v   <WireConstPtr>);
  
  /// "Reference" to a `geo::WireGeo` is actually... a copy of it???
  /// Meh... temporaries.
  using WireRef = WireGeo const;
  
  static_assert( std::is_copy_constructible_v<WireRef>);
  static_assert( std::is_move_constructible_v<WireRef>);
  static_assert(!std::is_copy_assignable_v   <WireRef>); // because constant
  static_assert(!std::is_move_assignable_v   <WireRef>); // because constant
  
  constexpr WirePtr InvalidWirePtr {}; // default constructed
  
  
  //----------------------------------------------------------------------------
  // traits specialization
  template <>
  struct element_traits<geo::WireGeo::ID_t> {
    
    /// The type representing the geometry of the element.
    using geometry_type = geo::WireGeo;
    
    /// The type representing the ID of the element.
    using id_type = geo::WireGeo::ID_t;
    
    /// Type used as reference to the geometry element object.
    using geometry_reference = WireRef;
    
    /// Type used as pointer to the geometry element object.
    using geometry_pointer = WirePtr;
    
    
  }; // geo::details::default_element_traits<>
  
  //----------------------------------------------------------------------------
  
  
} // namespace geo


//------------------------------------------------------------------------------
//--- inline implementation
//------------------------------------------------------------------------------
inline geo::WireGeo::WireGeo
  (geo::PlaneGeo const& plane, geo::WireID::WireID_t wireNo) noexcept
  : fData{ &plane, { wireNo } }
{}


//------------------------------------------------------------------------------
inline constexpr geo::WirePtr geo::WireGeo::operator& () const
  { return { this }; }


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void geo::WireGeo::PrintGenericWireInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1 */
) const {

  //----------------------------------------------------------------------------
  out << "wire from " << GetStart<geo::Point_t>()
    << " to " << GetEnd<geo::Point_t>();

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  out << " (" << Length() << " cm long)";

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------
  out << ", theta(z)=" << ThetaZ() << " rad";

  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------
  out << "\n" << indent
    << "  center at " << GetCenter<geo::Point_t>() << " cm";

  if (verbosity-- <= 0) return; // 3

  //----------------------------------------------------------------------------
  out << ", direction: " << Direction<geo::Vector_t>();

//  if (verbosity-- <= 0) return; // 4

  //----------------------------------------------------------------------------
} // geo::WireGeo::PrintGenericWireInfo()


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_WIREGEO_H
