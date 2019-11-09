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

#endif // LARCOREALG_COREUTILS_FROMFUTUREIMPORT_H
