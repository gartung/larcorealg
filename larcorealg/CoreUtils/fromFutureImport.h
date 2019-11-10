/**
 * @file   larcorealg/CoreUtils/fromFutureImport.h
 * @brief  Code that might appear as standard C++ in the future.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   April 25, 2019
 * 
 * This is currently a header-only library.
 */

#ifndef LARCOREALG_COREUTILS_FROMFUTUREIMPORT_H
#define LARCOREALG_COREUTILS_FROMFUTUREIMPORT_H

// C/C++ standard libraries
#include <string_view>
#include <string>
#include <utility> // std::forward()
#include <system_error> // std::errc


/**
 * @defgroup FutureStandards Features from future C++ standards
 * @brief Features expected to be provided by future C++ standards.
 */


/**
 * Namespace anticipating some simple features of C++ standards not yet adopted.
 *
 * It is recommended that whenever all supported compilers support each single
 * feature, that be removed from here, and the standard one immediately adopted.
 * 
 * The aim is that the interface and behaviour here are as similar as possible,
 * so that the update should boil down to a different header inclusion and a
 * different namespace.
 *
 * @addtogroup FutureStandards
 */ 
namespace util::pre_std {
  
#if (__cplusplus < 202000L) // still to be defined, should be C++20
  
  /// Transparent functor that returns its argument just as passed.
  struct identity {
    
    struct is_transparent {}; // STL algorithms will be happy to find this
    
    template <typename T>
    constexpr T&& operator() (T&& t) const noexcept
      { return std::forward<T>(t); }
    
  }; // identity<>
  
#endif // (__cplusplus < 202000L)
  
  
  // --- BEGIN -- charconv -----------------------------------------------------
  /**
   * @name A few functions reminiscent of the ones in `charconv` C++17 library
   * 
   * These functions are complex enough that I am not trying to faithfully
   * implement them here, but rather cheat my way through them.
   */
  /// @{
  
  struct from_chars_result {
    const char* ptr = nullptr;
    std::errc ec = std::errc{};
  }; // struct from_chars_result
  
  
  from_chars_result from_chars
    (const char* first, const char* last, double& value);
  
  from_chars_result from_chars
    (const char* first, const char* last, long double& value);
  
  from_chars_result from_chars
    (const char* first, const char* last, float& value);
  
  from_chars_result from_chars
    (const char* first, const char* last, int& value, int base = 10);
  
  from_chars_result from_chars
    (const char* first, const char* last, unsigned int& value, int base = 10);
  
  /// @}
  // --- END -- charconv -------------------------------------------------------
  
} // namespace pre_std


// -----------------------------------------------------------------------------
/**
 * @defgroup NonStandards Features excluded from C++ standards
 * @brief Features expected not to be provided by C++ standards, which should
 *        have been.
 */


/**
 * @brief Namespace implementing some simple features that haven't made it into
 *        C++ standards and possibly never will.
 *
 * @addtogroup NonStandards
 */ 
namespace util::not_std {}

/**
 * @brief Operations with `std::string_view`.
 * 
 * String and string view concatenation
 * -------------------------------------
 * 
 * For some reason (an hypothesis at
 * https://stackoverflow.com/questions/44636549) concatenation between
 * `std::string` and `std::string_view` is not supported in C++ standard (C++17
 * at least). However, `std::string::operator += ()`, `std::string::append()`,
 * `std::string::assign()`, `std::string::insert()`, ... are able to use
 * anything with iterator range interface, including a `std::string_view`.
 * 
 * An example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * void checkLength(std::string_view const& s, std::size_t const maxLength) {
 *   using namespace std::string_literals;
 *   using namespace util::not_std::string_view_ops;
 *   
 *   if (s.length() <= maxLength) return;
 *   
 *   throw std::logic_error(
 *     "String '"s + s + "' is "
 *     + std::to_string(s.length()) + " characters long!"
 *     );
 *   
 * } // void checkLength()
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * The concatenation is enabled by using the `util::not_std::string_view_ops`
 * namespace. This is a case where we don't have control on what to pass to
 * an interface (`std::logic_error` requires a `std::string`).
 * 
 */
namespace util::not_std::string_view_ops {
  
  /// Returns a new string containing the concatenation of `a` and `b`.
  template <typename CharT, typename Traits, typename Allocator>
  auto operator+ (
    std::basic_string<CharT, Traits, Allocator> const& a,
    std::basic_string_view<CharT, Traits> const& b
    )
    -> std::basic_string<CharT, Traits, Allocator>
  {
    std::basic_string<CharT, Traits, Allocator> ab;
    ab.reserve(a.size() + b.size()); // allocation of memory, once
    return ab.assign(a).append(b);
  }
  
  /// Returns a new string containing the concatenation of `a` and `b`.
  template <typename CharT, typename Traits, typename Allocator>
  auto operator+ (
    std::basic_string_view<CharT, Traits> const& a,
    std::basic_string<CharT, Traits, Allocator> const& b
    )
    -> std::basic_string<CharT, Traits, Allocator>
  {
    std::basic_string<CharT, Traits, Allocator> ab;
    ab.reserve(a.size() + b.size()); // allocation of memory, once
    return ab.assign(a).append(b);
  }
  
  
} // namespace util::not_std::string_view_ops


// -----------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_FROMFUTUREIMPORT_H
