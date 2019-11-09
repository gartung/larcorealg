  /**
 * @file   larcorealg/Geometry/geo_vectors_parse.h
 * @brief  Utitlties to parse strings into vectors.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 8, 2019
 * @see    `larcorealg/Geometry/geo_vectors_parse.cxx`
 * 
 */

#ifndef LARCOREALG_GEOMETRY_GEO_VECTORS_PARSE_H
#define LARCOREALG_GEOMETRY_GEO_VECTORS_PARSE_H

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcorealg/CoreUtils/fromFutureImport.h" // util::pre_std::from_chars
#include "larcorealg/CoreUtils/StdUtils.h" // util::begin(), util::end()

// C/C++ standard libraries
#include <regex>
#include <string>
#include <string_view>
// #include <charconv> // std::from_chars()
#include <stdexcept> // std::runtime_error
#include <system_error> // std::errc, std::make_error_condition()
#include <limits> // std::numeric_limits<>
#include <cstddef> // std::size_t
#include <cassert>


/// Utilities to parse vectors.
namespace geo::vect::parse {
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Class to parse a string into a vector.
   * @tparam Dim dimension of the supported vectors
   * 
   * This class operates with the `parse()` method.
   * It can parse vectors of any type which are compatible with the coordinate
   * managers (see `geo::vect::coord()`), but they are required to be of
   * dimension `Dim`.
   * 
   * Currently, the only supported format is:
   * `( "{" | "(" ) value (( ";" | "," ) value ...) ( "}" | ")" )`
   * with a `value` something that `std::from_chars()` can convert into a number
   * of the proper type for the target vector.
   * 
   * @note This object may be made moderarely more customizable on demand.
   */
  template <std::size_t Dim>
  class VectorParser {
    
    /// Returns the index of the match to the opening bracket of the vector.
    static constexpr std::size_t indexOfOpen() { return 1U; }
    
    /// Returns the index of the match to the `iCoord` coordinate of the vector.
    static constexpr std::size_t indexOfMatch(std::size_t const iCoord)
      { return indexOfOpen() + 1U + iCoord * 2U; };
    
    /// Returns the index of the match to the closing bracket of the vector.
    static constexpr std::size_t indexOfClose()
      { return indexOfMatch(dimension() - 1U) + 1U; }
    
    static std::string OpenOpts() { return "({"; }
    static std::string CloseOpts() { return ")}"; }
    static std::string Blanks() { return "[[:blank:]]*"; }
    static std::string Open() { return "([" + OpenOpts() + "])"; }
    static std::string RealNum() { return "([-+.\\deE]+)"; }
    static std::string Sep() { return "([,;[:blank:]])"; }
    static std::string Close() { return "([" + CloseOpts() + "])"; }
    
    std::string const fPattern; ///< The pattern to apply to the vector.
    std::regex const fRegEx; ///< The regular expression for parsing.
    
    /// Creates the matching pattern.
    static std::string buildPattern();
    
      public:
    
    /// Dimension of the supported vectors.
    static constexpr std::size_t dimension() { return Dim; }
    
    /// Constructor. Just that.
    VectorParser(): fPattern(buildPattern()), fRegEx(fPattern) {}
    
    /**
     * @brief Parses the specified string.
     * @tparam Vector (mandatory) the type of vector to create
     * @tparam String the type of input string
     * @param s the string to be parsed
     * @return an object of type `Vector` with content reflecting `s`
     * @throw std::runtime_error if conversion failed
     * 
     * The `Vector` type must be compatible with the coordinate managers
     * (see `geo::vect::coord()`) and have the presctibed `dimension()`.
     * 
     * The input `String` type must respond to the `begin()`/`end()` STL
     * interface and be a sequence of characters.
     */
    template <typename Vector, typename String>
    Vector parse(String const& s) const;
    
  }; // class VectorParser
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Converts a string to a vector.
   * @tparam Vector (mandatory) the type of vector to create
   * @tparam String the type of input string
   * @param s the string to be parsed
   * @return an object of type `Vector` with content reflecting `s`
   * @throw std::runtime_error if conversion failed
   * 
   * This function creates an ad hoc `VectorParser` object to convert `s`.
   */
  template <typename Vector, typename String>
  Vector parse(String const& s);
  
  
  // ---------------------------------------------------------------------------
  /**
   * @brief Converts a string to a vector.
   * @tparam Vector the type of vector to create
   * @tparam String the type of input string
   * @param v (output) the vector to be filled
   * @param s the string to be parsed
   * @throw std::runtime_error if conversion failed
   * 
   * This function creates an ad hoc `VectorParser` object to convert `s`
   * and stores the result of the conversion in `v`.
   */
  template <typename Vector, typename String>
  void parse(Vector& v, String const& s);
  
  
  // ---------------------------------------------------------------------------
  
} // namespace geo::vect::parse


// -----------------------------------------------------------------------------
// --- template implementation
// -----------------------------------------------------------------------------
namespace geo::vect::parse::details {
  // we don't want to "pollute" the global namespace with some useful things
  // overlooked by C++
  
  template <typename CharT, typename Traits, typename Allocator>
  auto operator+ (
    std::basic_string<CharT, Traits, Allocator> const& a,
    std::basic_string_view<CharT, Traits> const& b
    )
    -> std::basic_string<CharT, Traits, Allocator>
  {
    std::basic_string<CharT, Traits, Allocator> ab;
    ab.reserve(a.size() + b.size()); // allocation of memory, once
    return ab.assign(a).append(b.begin(), b.end());
  }
  
  template <typename CharT, typename Traits, typename Allocator>
  auto operator+ (
    std::basic_string_view<CharT, Traits> const& a,
    std::basic_string<CharT, Traits, Allocator> const& b
    )
    -> std::basic_string<CharT, Traits, Allocator>
  {
    std::basic_string<CharT, Traits, Allocator> ab;
    ab.reserve(a.size() + b.size()); // allocation of memory, once
    return ab.assign(a.begin(), a.end()).append(b);
  }
  
  template <typename CharT, typename Traits, typename Allocator>
  auto operator+= (
    std::basic_string<CharT, Traits, Allocator> const& a,
    std::basic_string_view<CharT, Traits> const& b
    )
    -> std::basic_string<CharT, Traits, Allocator>&
    { return a.append(b.begin(), b.end()); }
  
} // namespace geo::vect::parse::details


// -----------------------------------------------------------------------------
// --- VectorParser
// -----------------------------------------------------------------------------
template <std::size_t Dim>
std::string geo::vect::parse::VectorParser<Dim>::buildPattern() {
  
  using namespace std::string_literals;
  
  assert(OpenOpts().size() == CloseOpts().size());
  
  std::string pattern = "^"s + Blanks() + Open() + Blanks() + RealNum();
  for (std::size_t i = 1U; i < dimension(); ++i)
    pattern += Blanks() + Sep() +  Blanks() + RealNum();
  pattern += Blanks() + Close() + Blanks() + "$"s;
  
  return pattern;
} // geo::vect::parse::VectorParser<>::buildPattern()


// -----------------------------------------------------------------------------
template <std::size_t Dim>
template <typename Vector, typename String>
Vector geo::vect::parse::VectorParser<Dim>::parse(String const& s) const {

  using namespace std::string_literals;
  using namespace details;

  using Vector_t = Vector;
  using Coord_t = geo::vect::coordinate_t<Vector_t>;
  
  using Iter_t = std::decay_t<decltype(util::begin(s))>;
  using MatchResults_t = std::match_results<Iter_t>;
  
  static_assert(geo::vect::dimension<Vector_t>() == dimension(),
    "VectorParser requires vectors matching its own dimension.");
  
  MatchResults_t match;
  if (!std::regex_match(util::begin(s), util::end(s), match, fRegEx)) {
    throw std::runtime_error("VectorParser::parse(): "s
      + " input '" + s + "' does not match dim=" + util::to_string(dimension())
      + " vector pattern '" + fPattern + "'\n"
      );
  }
  
  //
  // check consistency of opening and closing brackets
  //
  std::size_t const iOpen = OpenOpts().find(match[indexOfOpen()]);
  assert(iOpen < OpenOpts().size());
  if (match[indexOfClose()] != CloseOpts()[iOpen]) {
    throw std::runtime_error("VectorParser::parse(): "s
      + "input '" + s + "' opens with '" + match[indexOfOpen()].str()
      + "' but closes with '" + match[indexOfClose()].str()
      + "' (expected: '" + CloseOpts()[iOpen] + "')\n"
      );
  }
  
  //
  // check consistency of separators
  //
  if constexpr (dimension() > 2U) {
    std::string const sep = match[indexOfMatch(0U) + 1U];
    for (std::size_t i = 1; i < dimension(); ++i) {
      std::size_t const sepIndex = indexOfMatch(i) - 1U;
      if (match[sepIndex] == sep) continue;
        throw std::runtime_error("VectorParser::parse(): "s
          + " input '" + s + "' sports inconsistent separators ('"
          + match[sepIndex].str() + "' and '" + sep + "')\n"
          );
    } // for
  } // if dimension > 2
  
  //
  // fill the result
  //
  Vector_t v; // this is the result
  for (std::size_t iCoord = 0U; iCoord < dimension(); ++iCoord) {
    
    std::string const valueStr= match[indexOfMatch(iCoord)];
    assert(!valueStr.empty());
    
    auto begin = valueStr.c_str();
    auto const end = begin + valueStr.length();
    if (*begin == '+') ++begin; // std::from_chars() doesn't support "+"
    
    auto coordValue = std::numeric_limits<Coord_t>::signaling_NaN();
    auto const res = util::pre_std::from_chars(begin, end, coordValue);
    
    if (res.ec != std::errc{}) { // conversion error
      auto err = std::make_error_condition(res.ec);
      throw std::runtime_error("VectorParser::parse(): "s
        + " can't parse coordinate #" + util::to_string(iCoord) + " of input '"
        + s + "' ('" + valueStr + "'): " + err.message() + "\n"
        );
    }
    
    geo::vect::coord(v, iCoord) = coordValue;
    
  } // for
  
  return v;
} // geo::vect::parse::VectorParser<>::parse()


// -----------------------------------------------------------------------------
template <typename Vector, typename String>
Vector geo::vect::parse::parse(String const& s) {
  constexpr std::size_t Dim = geo::vect::dimension<Vector>();
  return VectorParser<Dim>().template parse<Vector>(s);
} // geo::vect::parse::parse()


// -----------------------------------------------------------------------------
template <typename Vector, typename String>
void geo::vect::parse::parse(Vector& v, String const& s) {
  v = parse<Vector>(s);
} // geo::vect::parse::parse(Vector)


// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEO_VECTORS_PARSE_H
