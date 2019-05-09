////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/TenseWireGeo.h
/// \brief Encapsulate the geometry of a wire
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_TENSEWIREGEO_H
#define LARCOREALG_GEOMETRY_TENSEWIREGEO_H

// LArSoft
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect

// ROOT
#include "TVector3.h"
#include "Math/GenVector/Transform3D.h"

// C/C++ libraries
#include <vector>
#include <string>
#include <type_traits> // std::is_move_constructible<>, ...
#include <cmath> // std::sin(), ...


// forward declarations
class TGeoNode;


namespace geo {


  /** **************************************************************************
   * @brief Geometry description of a TPC wire.
   * @ingroup Geometry
   * 
   * @note This used to be the one and only `geo::WireGeo` class. Please don't
   *       be too harsh on the new name: it's not easy to pick one!!
   *       (also, since this name is not exposed, if you have a better one,
   *       propose it and it will have a chance to go in; it's a challenge!)
   * 
   * The wire is a single straight segment on a wire plane.
   * Different wires may be connected to the same readout channel. That is of
   * no relevance for the geometry description.
   *
   * The wire has a start and an end point. Their definition of them is related
   * to the other wires in the plane and to the TPC itself.
   *
   * The direction of increasing wire coordinate, defined in the wire plane,
   * is orthogonal to the wire direction and of course points to the direction
   * where the wire number within the plane increases. This direction is
   * indirectly defined when sorting the wires in the plane, which is done by
   * the plane (geo::PlaneGeo). This direction lies by definition on the wire
   * plane. The direction normal to the wire plane is defined by the TPC so that
   * it points inward the TPC rather than outward.
   * Finally, the wire direction is defined so that the triplet of unit vectors
   * direction of the wire @f$ \hat{l} @f$, direction of increasing wire number
   * @f$ \hat{w} @f$, and normal to the plane @f$ \hat{n} @f$ is positively
   * defined (@f$ \hat{l} \times \hat{w} \cdot \hat{n} = +1 @f$).
   * The start @f$ \vec{a}_{w} @f$ and the end of the wire @f$ \vec{b}_{w} @f$
   * are defined so that their difference @f$ \vec{b}_{w} - \vec{a}_{w} @f$
   * points in the same direction as @f$ \hat{l} @f$.
   *
   */
  class TenseWireGeo: public geo::WireGeo {

  public:

    using GeoNodePath_t = std::vector<TGeoNode const*>;

    /**
     * @brief Constructor from a ROOT geometry node and a transformation.
     * @param node ROOT geometry node
     * @param trans transformation matrix (local to world)
     *
     * The node describes the shape of the wire (the only relevant information
     * is in fact the length), while the transformation described its
     * positioning in the world (both position and orientation).
     *
     * A pointer to the node and a copy of the transformation matrix are kept
     * in the `WireGeo` object.
     */
    TenseWireGeo(TGeoNode const& node, geo::TransformationMatrix&& trans);
    
    
    /**
     * @brief Prints information of this wire.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
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
    template <typename Stream>
    void PrintTenseWireInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;
    
    
  protected:
    
    // --- BEGIN -- Polymorphic implementation ---------------------------------
    /// @name Polymorphic implementation
    /// @{
    
    virtual double doRMax() const override;
    virtual double doHalfL() const override { return fHalfL; }
    virtual double doRMin() const override;
    virtual void doGetCenter(double* xyz, double localz = 0.0) const override;
    virtual void doGetStart(double* xyz) const override
      { doGetCenter(xyz, -fHalfL); }
    virtual void doGetEnd(double* xyz) const override
      { doGetCenter(xyz, +fHalfL); }
    virtual geo::Point_t doGetPositionFromCenter(double localz) const override
      { return doGetPositionFromCenterUnbounded(capLength(localz)); }
    virtual geo::Point_t doGetPositionFromCenterUnbounded
      (double localz) const override
      { return toWorldCoords(LocalPoint_t{ 0.0, 0.0, relLength(localz) }); }
    virtual geo::Point_t doGetCenter() const override
      { return fCenter; }
    virtual geo::Point_t doGetStart() const override
      { return doGetPositionFromCenterUnbounded(-fHalfL); }
    virtual geo::Point_t doGetEnd() const override
      { return doGetPositionFromCenterUnbounded(+fHalfL); }
    virtual double doLength() const override
      { return 2.0 * fHalfL; }
    virtual double doThetaZ() const override
      { return fThetaZ; }
    virtual bool doIsParallelTo(geo::WireGeo const& wire) const override
      {
        return // parallel if the dot product of the directions is about +/- 1
          std::abs(std::abs(doDirection().Dot(wire.Direction<geo::Vector_t>())) - 1.) < 1e-5;
      }
    virtual geo::Vector_t doDirection() const override
      { return (doGetCenter() - doGetStart()) / doHalfL(); }
    virtual double doDistanceFrom(geo::WireGeo const& wire) const override;
    virtual TGeoNode const* doNode() const override
      { return fWireNode; }
    virtual double doComputeZatY0() const override
      { return fCenter.Z() - fCenter.Y() / doTanThetaZ(); }
    
    /// Returns the transformation from local to global coordinates.
    /// @returns the transformation object, or `nullptr` if doesn't exist
    virtual LocalTransformation_t const* doTrans() const override
      { return &fTrans; }
    
    virtual std::string doWireInfo
      (std::string indent = "", unsigned int verbosity = 1) const override;
    
    virtual void doUpdateAfterSorting
      (geo::WireID const& wireid, bool flip) override;
    
    /// @}
    // --- END -- Polymorphic implementation -----------------------------------
    

  private:

    const TGeoNode*    fWireNode;  ///< Pointer to the wire node
    double             fThetaZ;    ///< angle of the wire with respect to the z direction
    double             fHalfL;     ///< half length of the wire
    geo::Point_t       fCenter;    ///< Center of the wire in world coordinates.
    LocalTransformation_t fTrans;  ///< Wire to world transform.
    bool               flipped;    ///< whether start and end are reversed

    /// Returns whether ( 0, 0, fHalfL ) identifies end (false) or start (true)
    /// of the wire.
    bool isFlipped() const { return flipped; }

    /// Returns the relative length from center to be used when transforming.
    double relLength(double local) const { return isFlipped()? -local: local; }

    /// Caps the specified local length coordinate to lay on the wire.
    double capLength(double local) const
      { return std::min(+fHalfL, std::max(-fHalfL, local)); }

    /// Stacked `capLength()` and `relLength()`.
    double capRelLength(double local) const
      { return capLength(relLength(local)); }

    /// Set to swap the start and end wire
    void Flip();


  }; // class TenseWireGeo
  
  
  static_assert(!std::is_copy_assignable_v<geo::TenseWireGeo>);
  static_assert(!std::is_copy_constructible_v<geo::TenseWireGeo>);
  static_assert( std::is_move_assignable_v<geo::TenseWireGeo>);
  static_assert( std::is_move_constructible_v<geo::TenseWireGeo>);
  
  
} // namespace geo


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void geo::TenseWireGeo::PrintTenseWireInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1 */
) const {

  //----------------------------------------------------------------------------
  out << "wire from " << doGetStart() << " to " << doGetEnd();

  if (verbosity-- <= 0) return; // 0

  //----------------------------------------------------------------------------
  out << " (" << Length() << " cm long)";

  if (verbosity-- <= 0) return; // 1

  //----------------------------------------------------------------------------
  out << ", theta(z)=" << ThetaZ() << " rad";

  if (verbosity-- <= 0) return; // 2

  //----------------------------------------------------------------------------
  out << "\n" << indent
    << "  center at " << doGetCenter() << " cm";

  if (verbosity-- <= 0) return; // 3

  //----------------------------------------------------------------------------
  out << ", direction: " << doDirection();

//  if (verbosity-- <= 0) return; // 4

  //----------------------------------------------------------------------------
} // geo::TenseWireGeo::PrintTenseWireInfo()


//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_TENSEWIREGEO_H
