/**
 * @file    larcorealg/Geometry/WireGeo.h
 * @brief   Interface for a active readout element on a TPC plane.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    April 24, 2019
 * @see     larcorealg/Geometry/WireGeo.cxx
 * @ingroup Geometry
 */

#ifndef LARCOREALG_GEOMETRY_WIREGEO_H
#define LARCOREALG_GEOMETRY_WIREGEO_H

// LArSoft
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect
#include "larcorealg/CoreUtils/UncopiableAndUnmovableClass.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::WireID
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::*

// ROOT
#include "TVector3.h"
#include "Math/GenVector/Transform3D.h"

// C/C++ libraries
#include <vector>
#include <string>
#include <cmath> // std::sin(), ...


// forward declarations
class TGeoNode;
class TVector3;


namespace geo {


  /** **************************************************************************
   * @brief Interface for a geometry description of a TPC readout element.
   * @ingroup Geometry
   * 
   * A `geo:WireGeo` represents an element of a TPC plane that is connected to
   * some readout, and is therefore instrumental to the creation of a waveform.
   * While the legacy name of this interface is suggestive of a physical wire,
   * this object can as well represent a deposited strip (as in some dual-phase
   * LAr TPCs) or a pixel. Nevertheless, for historical reasons, the main
   * interface is worded and designed to well describe an actual wire.
   * 
   * Different "wires" may be connected to the same readout channel. That is of
   * no relevance for the geometry description.
   *
   * The wire object has a start and an end point. Their definition of them is
   * related to the other wires objects in the plane and to the TPC itself.
   */
  class WireGeo: private lar::PolymorphicUncopiableClass {
    
    using DefaultVector_t = TVector3; // ... not for long
    using DefaultPoint_t = TVector3; // ... not for long

  public:
    
    // --- BEGIN -- Types for geometry-local reference vectors -----------------
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

    ///@}
    // --- END -- Types for geometry-local reference vectors -------------------
    
    
    // --- BEGIN -- Size and coordinates ---------------------------------------
    /// @{
    /// @name Size and coordinates

    /// Returns the outer half-size of the wire [cm]
    double RMax() const { return doRMax(); }

    /// Returns half the length of the wire [cm]
    double HalfL() const { return doHalfL(); }

    /// Returns the inner radius of the wire (usually 0) [cm]
    double RMin() const { return doRMin(); }

    /**
     * @brief Fills the world coordinate of a point on the wire
     * @param xyz _(output)_ the position to be filled, as [ x, y, z ] (in cm)
     * @param localz distance of the requested point from the middle of the wire
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
    void GetCenter(double* xyz, double localz=0.0) const
      { doGetCenter(xyz, localz); }

    /// Fills the world coordinate of one end of the wire
    /// @deprecated Use the version returning a vector instead.
    void GetStart(double* xyz) const { doGetStart(xyz); }

    /// Fills the world coordinate of one end of the wire
    /// @deprecated Use the version returning a vector instead.
    void GetEnd(double* xyz) const { doGetEnd(xyz); }

    //@{
    /**
     * @brief Returns the position (world coordinate) of a point on the wire
     * @tparam Point type of vector to be returned (current default: `TVector3`)
     * @param localz distance of the requested point from the middle of the wire
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
      { return geo::vect::convertTo<Point>(doGetPositionFromCenter(localz)); }
    DefaultPoint_t GetPositionFromCenter(double localz) const
      { return GetPositionFromCenter<DefaultPoint_t>(localz); }
    //@}

    //@{
    /**
     * @brief Returns the position (world coordinate) of a point on the wire
     * @tparam Point type of vector to be returned (current default: `TVector3`)
     * @param localz distance of the requested point from the middle of the wire
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
        return 
          geo::vect::convertTo<Point>(doGetPositionFromCenterUnbounded(localz)); 
      }
    DefaultPoint_t GetPositionFromCenterUnbounded(double localz) const
      { return GetPositionFromCenterUnbounded<DefaultPoint_t>(localz); }
    //@}

    //@{
    /// Returns the world coordinate of the center of the wire [cm]
    /// @tparam Point type of the point being returned
    template <typename Point>
    Point GetCenter() const
      { return geo::vect::convertTo<Point>(doGetCenter()); }
    DefaultPoint_t GetCenter() const { return GetCenter<DefaultPoint_t>(); }
    //@}

    //@{
    /// Returns the world coordinate of one end of the wire [cm]
    /// @tparam Point type of the point being returned
    template <typename Point>
    Point GetStart() const { return geo::vect::convertTo<Point>(doGetStart()); }
    DefaultPoint_t GetStart() const { return GetStart<DefaultPoint_t>(); }
    //@}

    //@{
    /// Returns the world coordinate of one end of the wire [cm]
    /// @tparam Point type of the point being returned
    template <typename Point>
    Point GetEnd() const { return geo::vect::convertTo<Point>(doGetEnd()); }
    DefaultPoint_t GetEnd() const { return GetEnd<DefaultPoint_t>(); }
    //@}

    /// Returns the wire length in centimeters
    double Length() const { return doLength(); }

    /// Returns the z coordinate, in centimetres, at the point where y = 0.
    /// Assumes the wire orthogonal to x axis and the wire not parallel to z.
    double ComputeZatY0() const { return doComputeZatY0(); }

    /**
     * @brief Returns 3D distance from the specified wire
     * @return the signed distance in centimetres (0 if wires are not parallel)
     *
     * If the specified wire is "ahead" in z respect to this, the distance is
     * returned negative.
     */
    double DistanceFrom(geo::WireGeo const& wire) const
      { return doDistanceFrom(wire); }
    
    /// @}
    // --- END -- Size and coordinates -----------------------------------------
    
    
    // --- BEGIN -- Orientation and angles -------------------------------------
    /// @{
    /// @name Orientation and angles

    /// Returns angle of wire with respect to z axis in the Y-Z plane [radians]
    double ThetaZ() const { return doThetaZ(); }

    /**
     * Returns angle of wire with respect to _z_ axis.
     * @param degrees return the angle in degrees rather than radians
     * @return wire angle
     */
    double ThetaZ(bool degrees) const
      { return degrees? util::RadiansToDegrees(ThetaZ()): ThetaZ(); }
    
    
    //@{
    /// Returns trigonometric operations on `ThetaZ()`.
    double CosThetaZ() const { return doCosThetaZ(); }
    double SinThetaZ() const { return doSinThetaZ(); }
    double TanThetaZ() const { return doTanThetaZ(); }
    //@}

    /// Returns if this wire is horizontal (@f$ \theta_{z} \approx 0 @f$).
    /// @deprecated The method name is confusing, and it was unused as of `v08_16_00` anyway.
    [[deprecated]] bool isHorizontal() const { return std::abs(SinThetaZ()) < 1e-5; }

    /// Returns if this wire is vertical (@f$ \theta_{z} \approx \pi/2 @f$).
    /// @deprecated The method name is confusing, and it was unused as of `v08_16_00` anyway.
    [[deprecated]] bool isVertical() const { return std::abs(CosThetaZ()) < 1e-5; }

    /// Returns if this wire is parallel to another.
    bool isParallelTo(geo::WireGeo const& wire) const
      { return doIsParallelTo(wire); }

    //@{
    /// Returns the wire direction as a norm-one vector.
    /// @tparam Vector type of the vector being returned
    template <typename Vector>
    Vector Direction() const
      { return geo::vect::convertTo<Vector>(doDirection()); }
    DefaultVector_t Direction() const { return Direction<DefaultVector_t>(); }
    //@}
    /// @}
    // --- END -- Orientation and angles ---------------------------------------
    
    // --- BEGIN -- Coordinate conversion --------------------------------------
    /// @{
    /**
     * @name Coordinate conversion
     *
     * Local points and displacement vectors are described by the types
     * `geo::WireGeo::LocalPoint_t` and `geo::WireGeo::LocalVector_t`,
     * respectively.
     */

    /// Transform point from local wire frame to world frame.
    void LocalToWorld(const double* wire, double* world) const
      { doTrans()->LocalToWorld(wire, world); }

    /// Transform point from local wire frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
      { return doTrans()->toWorldCoords(local); }

    /// Transform direction vector from local to world.
    void LocalToWorldVect(const double* wire, double* world) const
      { doTrans()->LocalToWorldVect(wire, world); }

    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
      { return doTrans()->toWorldCoords(local); }

    /// Transform point from world frame to local wire frame.
    void WorldToLocal(const double* world, double* wire) const
      { doTrans()->WorldToLocal(world, wire); }

    /// Transform point from world frame to local wire frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
      { return doTrans()->toLocalCoords(world); }

    /// Transform direction vector from world to local.
    void WorldToLocalVect(const double* world, double* wire) const
      { doTrans()->WorldToLocalVect(world, wire); }

    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
      { return doTrans()->toLocalCoords(world); }

    /// @}
    // --- END -- Coordinate conversion ----------------------------------------


    /**
     * @brief Returns the geometry node representing this wire.
     * @return a pointer to the `TGeoNode`, or `nullptr` if not available.
     * 
     * If the geometry description does not explicitly include this element
     * (for example because wires are created by procedure), the returned
     * value is `nullptr`.
     */
    TGeoNode const* Node() const { return doNode(); }

    
    // --- BEGIN -- Printout of wire information -------------------------------
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
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1) const
      { out << WireInfo(indent, verbosity); }

    /**
     * @brief Returns a string with all the information of the wire.
     * @see `PrintWireInfo()`
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
    std::string WireInfo
      (std::string indent = "", unsigned int verbosity = 1) const
      { return doWireInfo(indent, verbosity); }

    /**
     * @brief Returns a string with all the information of the wire.
     * @see `PrintWireInfo()`
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
    std::string GenericWireInfo
      (std::string indent = "", unsigned int verbosity = 1) const;
    
    
    /// Maximum verbosity supported by `PrintWireInfo()`.
    static constexpr unsigned int MaxVerbosity = 4;
    
    /// @}
    // --- END -- Printout of wire information ---------------------------------
    
    
    /// Internal updates after the relative position of the wire is known
    /// (currently no-op)
    void UpdateAfterSorting(geo::WireID const& wireid, bool flip)
      { doUpdateAfterSorting(wireid, flip); }
    
    
    /// Returns the pitch (distance on y/z plane) between two wires, in cm
    static double WirePitch(geo::WireGeo const& w1, geo::WireGeo const& w2)
      { return std::abs(w2.DistanceFrom(w1)); }
    
  protected:
  
    using LocalTransformation_t = geo::LocalTransformationGeo
      <ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;
    
    
    // --- BEGIN -- Polymorphic implementation ---------------------------------
    /// @name Polymorphic implementation
    /// @{
    
    virtual double doRMax() const = 0;
    virtual double doHalfL() const = 0;
    virtual double doRMin() const = 0;
    virtual void doGetCenter(double* xyz, double localz = 0.0) const = 0;
    virtual void doGetStart(double* xyz) const = 0;
    virtual void doGetEnd(double* xyz) const = 0;
    virtual geo::Point_t doGetPositionFromCenter(double localz) const = 0;
    virtual geo::Point_t doGetPositionFromCenterUnbounded(double localz) const = 0;
    virtual geo::Point_t doGetCenter() const = 0;
    virtual geo::Point_t doGetStart() const = 0;
    virtual geo::Point_t doGetEnd() const = 0;
    virtual double doLength() const = 0;
    virtual double doThetaZ() const = 0;
    virtual double doCosThetaZ() const { return std::cos(doThetaZ()); }
    virtual double doSinThetaZ() const { return std::sin(doThetaZ()); }
    virtual double doTanThetaZ() const { return std::tan(doThetaZ()); }
    virtual bool doIsParallelTo(geo::WireGeo const& wire) const = 0;
    virtual geo::Vector_t doDirection() const = 0;
    virtual double doDistanceFrom(geo::WireGeo const& wire) const = 0;
    virtual TGeoNode const* doNode() const = 0;
    virtual double doComputeZatY0() const = 0;
    
    /// Returns the transformation from local to global coordinates.
    /// @returns the transformation object, or `nullptr` if doesn't exist
    virtual LocalTransformation_t const* doTrans() const = 0;
    
    virtual std::string doWireInfo
      (std::string indent = "", unsigned int verbosity = 1) const
      { return GenericWireInfo(indent, verbosity); }
    
    virtual void doUpdateAfterSorting(geo::WireID const& wireid, bool flip) {}
    
    /// @}
    // --- END -- Polymorphic implementation -----------------------------------
    
    
    /// Throws a `cet::exception` with a short call stack dump.
    /// @throws cet::exception (`"NotImplemented"`) with a short call stack dump
    [[noreturn]] void NotImplemented() const;
    
    
  }; // class WireGeo


} // namespace geo


//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_WIREGEO_H
