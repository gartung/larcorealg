/**
 * @file    larcorealg/Geometry/PixelPlane/PixelCoord.h
 * @brief   Representation of a pixel coordinate.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 30, 2019
 * @ingroup Geometry
 * 
 * This library is header only.
 */

#ifndef LARCOREALG_GEOMETRY_PIXELPLANE_PIXELCOORD_H
#define LARCOREALG_GEOMETRY_PIXELPLANE_PIXELCOORD_H


// C/C++ standard libraries
#include <array>
#include <cstddef> // std::size_t


/// Data structures for description of pixels.
namespace geo::pixel {
  
  using DirIndex_t = std::size_t; ///< Type used for identifying a direction.

  /**
   * @name Mnemonic codes for the indices used to identify coordinates.
   * 
   * The secondary axis (`ixSec`) is the one labeled "wire" or "wire coordinate"
   * in the wire plane legacy interface, while the main one is labeled "wire
   * direction".
   */
  /// @{
  static constexpr DirIndex_t ixMain  = 0U; ///< Index for main coordinate.
  static constexpr DirIndex_t ixSec   = 1U; ///< Index for secondary coordinate.
  static constexpr DirIndex_t NCoords = 2U; ///< Number of coordinates.

  static constexpr DirIndex_t ixWireC = ixSec; ///< Index for "wire" coordinate.
  static constexpr DirIndex_t ixWireD = ixMain; ///< Index for "wire" direction.
  /// @}


  /// Type of coordinate for a pixel on the plane.
  template <typename T>
  class PixelCoordT;
  
} // namespace geo::pixel

// -----------------------------------------------------------------------------
/// Type of coordinate for a pixel on the plane.
template <typename T>
class geo::pixel::PixelCoordT {
  
    public:
  using Coord_t = T; ///< Type of coordinate sported.
  
  
  // inherit pair constructors
  
  /// Default constructor: coordinates set to 0.
  constexpr PixelCoordT() { fCoords.fill(Coord_t{ 0 }); }
  /// Constructor: sets the specified coordinates.
  constexpr PixelCoordT(Coord_t main, Coord_t second)
    : fCoords({ main, second }) {}
  
  
  //@{
  /// Returns the main coordinate index.
  constexpr Coord_t main() const { return get<ixMain>(); }
  Coord_t& main() { return get<ixMain>(); }
  //@}
  
  //@{
  /// Returns the secondary coordinate index.
  constexpr Coord_t secondary() const { return get<ixSec>(); }
  Coord_t& secondary() { return get<ixSec>(); }
  //@}
  
  //@{
  /// Coordinate access by index.
  constexpr Coord_t operator[] (std::size_t i) const { return fCoords[i]; }
  Coord_t& operator[] (std::size_t i) { return fCoords[i]; }
  //@}
  
  //@{
  /// Coordinate access by static index.
  template <std::size_t I>
  constexpr Coord_t get() const { return std::get<I>(fCoords); }
  template <std::size_t I>
  Coord_t& get() { return std::get<I>(fCoords); }
  //@}
  
  
    private:
  std::array<Coord_t, NCoords> fCoords; ///< Value of the two coordinates.
  
}; // PixelCoordT<>


// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_PIXELPLANE_PIXELCOORD_H
