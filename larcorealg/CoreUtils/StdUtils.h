/**
 * @file   larcorealg/CoreUtils/StdUtils.h
 * @brief  Functions pulling in STL customization if available.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   July , 2019
 */

#ifndef LARCOREALG_COREUTILS_STDUTILS_H
#define LARCOREALG_COREUTILS_STDUTILS_H

// C/C++ standard libraries
#include <string> // std::to_string()
#include <algorithm> // std::equal()
#include <iterator> // std::begin(), std::end(), std::next()...

namespace util {
  
  /**
   * @defgroup LArSoft_CoreUtils_StdUtils C++ STL customizations
   * @brief C++ standard library customization for user-defined classes.
   * 
   * There are a number of functions that are provided by C++ standard library
   * for the data types and classes defined in the standard.
   * It is often desirable to have your class react to these standard functions
   * in a standard way, for example for a container to react to `std::begin()`
   * to return its `begin()` iterator. While sometimes this is easy (for example
   * `std::begin()` calls `begin()` member function if available), some other
   * times that is not possible. In that case, since overloading of functions in
   * the `std` namespace is not allowed by C++, the usual pattern is to rely on
   * the argument-dependent lookup (known as "ADL") to have the comnpiler find
   * the overloaded function that is defined in the same namespace as any of
   * the arguments. For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * userns::MyObj obj;
   * // ...
   * using std::to_string;
   * std::string objstr = to_string(obj);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will look in the namespace where the type of `obj` is defined (that is
   * `userns`) for a `userns::to_string`, then will consider `std::to_string`
   * it self.
   * 
   * The utilities provided here provide a transparent way to do that, at the
   * cost of a new header and some non-standard call. The equivalent call of the
   * above would be:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * userns::MyObj obj;
   * // ...
   * std::string objstr = util::to_string(obj);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * @note For customization of templates, like `std::hash` or
   *       `std::numeric_limits`, specialization of classes in `std` is allowed
   *       by the standard, so no particular trick is required.
   */
  /// @{
  
  
  /// ADL-aware version of `std::to_string`.
  template <typename T>
  constexpr decltype(auto) to_string(T&& obj)
    { using std::to_string; return to_string(std::forward<T>(obj)); }
  
  
  // --- BEGIN --- Containers and iterators ------------------------------------
  /// ADL-aware version of `std::begin`.
  template <typename T>
  constexpr decltype(auto) begin(T&& obj)
    { using std::begin; return begin(std::forward<T>(obj)); }
  
  /// ADL-aware version of `std::end`.
  template <typename T>
  constexpr decltype(auto) end(T&& obj)
    { using std::end; return end(std::forward<T>(obj)); }
  
  /// ADL-aware version of `std::cbegin`.
  template <typename T>
  constexpr decltype(auto) cbegin(T&& obj)
    { using std::cbegin; return cbegin(std::forward<T>(obj)); }
  
  /// ADL-aware version of `std::cend`.
  template <typename T>
  constexpr decltype(auto) cend(T&& obj)
    { using std::cend; return cend(std::forward<T>(obj)); }
  
  /// ADL-aware version of `std::size`.
  template <typename T>
  constexpr decltype(auto) size(T&& obj)
    { using std::size; return size(std::forward<T>(obj)); }
  
  /// ADL-aware version of `std::empty`.
  template <typename T>
  constexpr decltype(auto) empty(T&& obj)
    { using std::empty; return empty(std::forward<T>(obj)); }
  
  // --- END --- Containers and iterators --------------------------------------
  
  
  
  // --- BEGIN --- String operations -------------------------------------------
  /**
   * @brief Returns whether the string `s` starts with the string `prefix`.
   * @tparam String string type for the string to be checked
   * @tparam Prefix string type for the prefix
   * @param s string to be checked
   * @param prefix prefix to be sought
   * @return whether `s` starts exactly with `prefix`
   * 
   * Requirements:
   *  * `String` and `Prefix` types must support forward iterators
   *  * the character types in `String` and in `Prefix` must be comparable
   */
  template <typename String, typename Prefix>
  bool starts_with(String const& s, Prefix const& prefix);
  
  /**
   * @brief Returns whether the string `s` ends with the string `suffix`.
   * @tparam String string type for the string to be checked
   * @tparam Suffix string type for the suffix
   * @param s string to be checked
   * @param suffix suffix to be sought
   * @return whether `s` ends exactly with `suffix`
   * 
   * Requirements:
   *  * `String` and `Suffix` types must support forward iterators
   *  * the character types in `String` and in `Suffix` must be comparable
   */
  template <typename String, typename Suffix>
  bool ends_with(String const& s, Suffix const& suffix);
  
  /// @}
  // --- END --- String operations ---------------------------------------------
  
  
} // namespace util


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename String, typename Prefix>
bool util::starts_with(String const& s, Prefix const& prefix) {
  return std::equal(
    util::begin(s), std::next(util::begin(s), prefix.length()),
    util::begin(prefix), util::end(prefix)
    );
} // util::starts_with()


// -----------------------------------------------------------------------------
template <typename String, typename Suffix>
bool util::ends_with(String const& s, Suffix const& suffix)
{
  return (s.length() >= suffix.length()) && std::equal(
    std::next(util::begin(s), s.length() - suffix.length()), util::end(s),
    util::begin(suffix), util::end(suffix)
    );
} // util::ends_with()


// -----------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_STDUTILS_H

