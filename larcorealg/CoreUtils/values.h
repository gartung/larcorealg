/**
 * @file   larcorealg/CoreUtils/values.h
 * @brief  Definition of `util::values()` and `util::const_values()`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   May 8, 2019
 * 
 * This is a header-only library.
 */

#ifndef LARCOREALG_COREUTILS_VALUES_H
#define LARCOREALG_COREUTILS_VALUES_H


// LArSoft libraries
#include "larcorealg/CoreUtils/span.h" // util::make_transformed_span()

// C/C++ libraries
#include <map>
#include <unordered_map>
#include <utility> // std::forward(), std::as_const()
#include <tuple> // std::get()
#include <type_traits>
#include <cstddef> // std::size_t


namespace util {
  
  
  // -- BEGIN -- Transformed iterations ----------------------------------------
  /// @name Transformed iterations
  /// @{
  
  /**
   * @brief Range-for loop helper iterating across the values of the specified
   *        collection.
   * @tparam Coll type of the collection to iterate through
   * @param coll the collection to iterate through
   * @return an object suitable for range-for loop
   * @see `util::const_values()`
   * 
   * This function is in most of cases a no-operation, returning the collection
   * just as it was specified, to be iterated on directly.
   * In case of mapping types, though, a different object is returned and the
   * iteration will happen to the value type of the mapping instead than on
   * the key-value pair.
   * 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::map<int, float> data { { 1, 4.0F }, { 3, 12.0F }, { 2, 8.0F } };
   * std::vector<float> values;
   * 
   * for (float value: util::values(data)) 
   *   values.push_back(value);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will result in `values` vector being of size `3` and with values
   * `{ 4.0F, 8.0F, 12.0F }` (the order is the one of iterating through a
   * `std::map`).
   * 
   */
  template <typename Coll>
  decltype(auto) values(Coll&& coll);
  
  /**
   * @brief Range-for loop helper iterating across the constant values of the
   *        specified collection.
   * @see `util::values()`
   * 
   * This function is equivalent to `util::values()` but the values are
   * extracted as if the specified collection were constant.
   */
  template <typename Coll>
  decltype(auto) const_values(Coll&& coll);
  
  
  /// @}
  // -- END -- Transformed iterations ------------------------------------------
  
  
} // namespace util


//==============================================================================
//=== template implementation
//==============================================================================
//------------------------------------------------------------------------------
//--- util::values()
//------------------------------------------------------------------------------
namespace util::details {
  
  //----------------------------------------------------------------------------
  template <typename Coll, typename = void>
  struct values_impl {
    
    template <typename T>
    static constexpr decltype(auto) iterate(T&& coll) noexcept
      { return coll; }
    
  }; // struct values_impl
  
  
  //----------------------------------------------------------------------------
  template <typename Map>
  struct map_values_impl {
    
    static constexpr std::size_t NElement = 1U; // this is the value
    
    template <typename T>
    static constexpr decltype(auto) iterate(T&& coll) noexcept
      {
        auto extractor = [](auto&& value) -> decltype(auto)
          { return std::get<NElement>(value); };
        return util::make_transformed_span(coll, extractor);
      }
    
  }; // map_values_impl
  
  
  //----------------------------------------------------------------------------
  template <typename Key, typename Value, typename... Args>
  struct values_impl<std::map<Key, Value, Args...>>
    : map_values_impl<std::map<Key, Value, Args...>>
  {};
  template <typename Key, typename Value, typename... Args>
  struct values_impl<std::unordered_map<Key, Value, Args...>>
    : map_values_impl<std::unordered_map<Key, Value, Args...>>
  {};
  
  
  //----------------------------------------------------------------------------
  
} // namespace util::details


//------------------------------------------------------------------------------
template <typename Coll>
decltype(auto) util::values(Coll&& coll) {
  using coll_t = std::remove_reference_t<Coll>;
  return details::values_impl<coll_t>::iterate(std::forward<Coll>(coll)); 
} // util::values()


//------------------------------------------------------------------------------
template <typename Coll>
decltype(auto) util::const_values(Coll&& coll)
  { return values(std::as_const(coll)); }


//------------------------------------------------------------------------------


#endif // LARCOREALG_COREUTILS_VALUES_H
