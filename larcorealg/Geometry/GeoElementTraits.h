/**
 * @file   larcorealg/Geometry/GeoElementTraits.h
 * @brief  Classes generalizing traits of geometry elements.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 26, 2019
 *
 * Header only library.
 */

#ifndef LARCOREALG_GEOMETRY_GEOELEMENTTRAITS_H
#define LARCOREALG_GEOMETRY_GEOELEMENTTRAITS_H


// C/C++ standard library
#include <type_traits>


namespace geo {
  
  
  //----------------------------------------------------------------------------
  namespace details {
    
    //--------------------------------------------------------------------------
    template <typename ElementID, typename = void>
    struct default_element_traits;
    
    
    //--------------------------------------------------------------------------
    
  } // namespace details
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Traits of a LArSoft detector geometry element.
   * @tparam ElementID type of identifier of the element the traits pertain to
   * 
   * The relevant geometry elements are currently:
   *  * `geo::CryostatID`
   *  * `geo::TPCID`
   *  * `geo::PlaneID`
   *  * `geo::WireID`
   * 
   */
  template <typename ElementID>
  struct element_traits: details::default_element_traits<ElementID> {};
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Default traits for elements with a standard geometry element class.
   * @tparam ElementGeo geometry element class
   * @tparam ElementID geometry identifier type (used for consistency check)
   * 
   * This can be used as model to implement `geo::element_traits`:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct geo::element_traits<geo::CryostatID>
   *   : geo::default_geo_element_traits<geo::CryostatGeo, geo::CryostatID> {};
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename ElementGeo, typename ElementID = void>
  struct default_geo_element_traits;
  
  
  //----------------------------------------------------------------------------
  
  
} // namespace geo



//------------------------------------------------------------------------------
//--- implementation details
//------------------------------------------------------------------------------
template <typename ElementGeo, typename ElementID /* = void */>
struct geo::default_geo_element_traits {
  
  /// The type representing the geometry of the element.
  using geometry_type = ElementGeo;
  
  /// The type representing the ID of the element.
  using id_type = std::decay_t<decltype(std::declval<geometry_type>().ID())>;
  
  /// Type used as reference to the geometry element object.
  using geometry_reference = std::add_lvalue_reference_t<geometry_type const>;
  
  /// Type used as pointer to the geometry element object.
  using geometry_pointer = std::add_pointer_t<geometry_type const>;
  
  
  static_assert(
    std::is_void_v<ElementID>
    || std::is_same_v
      <ElementID, typename default_geo_element_traits<ElementGeo>::id_type>,
    "ID type specification does not match geometry element ID()"
    );
  
}; // geo::details::default_element_traits<>
  
  
//------------------------------------------------------------------------------
// specialization for IDs not represented by a geometry object
// (e.g. `readout::ROPID`).
template <typename ElementID, typename /* = void */>
struct geo::details::default_element_traits {
  
  /// The type representing the geometry of the element.
  using geometry_type = void;
  
  /// The type representing the ID of the element.
  using id_type = ElementID;
  
  
  struct geometry_reference; // geometry_reference not defined
  
  struct geometry_pointer; // geometry_pointer not defined
  
}; // geo::details::default_element_traits<ElementID>
  
  
//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEOELEMENTTRAITS_H
