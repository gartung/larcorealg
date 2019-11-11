/**
 * @file    larcorealg/Geometry/GeoMetadataParser.h
 * @brief   Helper to parse metadata from ROOT geometry hierarchy.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    November 10, 2019
 * @see     `larcorealg/Geometry/GeoMetadataParser.cxx`
*/

#ifndef LARCOREALG_GEOMETRY_GEOMETADATAPARSER_H
#define LARCOREALG_GEOMETRY_GEOMETADATAPARSER_H

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_parse.h" // geo::vect::parse namespace
#include "larcorealg/CoreUtils/StdUtils.h" // util::ends_with()
#include "larcorealg/CoreUtils/fromFutureImport.h" // util::pre_std::from_chars

// support libraries
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TMap.h"
#include "TObjString.h"
#include "TString.h"

// C/C++ standard library
#include <string_view>
#include <string>
#include <optional>
#include <utility> // std::pair
#include <cctype> // std::isblank()


// -----------------------------------------------------------------------------
namespace geo { class GeoMetadataParser; }

/**
 * @brief Helper class to parse metadata from ROOT geometry hierarchy.
 * 
 * This helper class facilitates the parsing of metadata put into a `TGeoNode`
 * or `TGeoVolume` by the ROOT GDML parser.
 * 
 * Each class instance is bound to a metadata object (`TMap`), and it can be
 * asked to look for specific keys:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * geo::GeoMetadataParser parser(metadata);
 * 
 * std::optional<unsigned int> nails = parser.getValue<unsigned int>("nails");
 * std::optional<double> nailLength
 *   = parser.getValue<double>("nailLength", "m");
 * std::optional<geo::Vector_t> nailDirection
 *   = parser.getVector<geo::Vector_t>("nailDir");
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * In this example, three optional values are retrieved.
 * The optional objects have values only if there was a corresponding metadata
 * entry (look up the documentation of `std::optional` for details).
 * If a base unit is specified (e.g. `"m"`), the metadata item is required to
 * have a unit specification derived from that base unit (e.g. `"cm"`).
 * The value will be returned in base units, e.g. if the metadata is specifying
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * <auxiliary auxtype="nailLength" auxvalue="13" auxunit="cm" />
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * (13 cm), the returned value will be `0.13`. Avoid using scaled units for
 * base units (e.g. `"kg"`), as it will not work as it should unless the unit
 * specified in the metadata is exactly the base unit.
 * If the base unit is not specified, the metadata must *not* have a unit, or
 * else the parsing will fail.
 * 
 * Currently, the parsing of two types of values is supported:
 * * scalar values, via `getValue()`
 * * vector values, via `getVector()`
 * 
 * In both cases, units are supported.
 * 
 * 
 * About metadata
 * ---------------
 * 
 * The metadata is placed in the ROOT geometry hierarchy by the ROOT GDML
 * parser. This parser extracts from each volume all the `auxiliary` elements,
 * in a form like:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * <auxiliary auxtype="colour" auxvalue="( 0, 0, 1 )" auxunit="RGB" />
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * which results in a metadata entry in the map with key `"colour"` and value
 * `"( 0, 0, 1 )*RGB"` (stored via `TObjString` objects).
 * 
 * This class actually expects the metadata map to be given to it.
 * To extract the metadata from a `TGeoNode` or `TGeoVolume`, one can do:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto pExt = dynamic_cast<TGeoRCExtension const*>(node.GetUserExtension());
 * TMap const* metadata
 *   = pExt? dynamic_cast<TMap const*>(pExt->GetUserObject()): nullptr;
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * which will have `metadata` host the map of metadata, or `nullptr` if either
 * the node has no metadata or the metadata is not in the form of a `TMap`.
 * 
 */
class geo::GeoMetadataParser {
  
  TMap const* fMetadata; ///< A reference to the metadata to be parsed.
  
  /// Returns the `TObjString` for the specified key, if any.
  TObjString const* getMetadataString(std::string const& key) const
    { return getMetadataString(*fMetadata, key); }

  /// Extracts from `metadata` the `TObjString` for the specified key, if any.
  static TObjString const* getMetadataString
    (TMap const& metadata, std::string const& key);

    public:
  
  /**
   * @brief Constructor: parses the specified metadata map.
   * @param metadata the map with the metadata to be parsed
   * 
   * The specified `metadata` must stay valid and accessible while this object
   * operates on it (i.e., it's referenced to rather than copied).
   */
  GeoMetadataParser(TMap const& metadata): fMetadata(&metadata) {}
  
  
  /**
   * 
   * 
   * 
   */
  template <typename T = double>
  std::optional<T> getValue(std::string const& key) const
    { return getValue<T>(*fMetadata, key); }
  
  template <typename T = double, typename UnitT = double>
  std::optional<T> getValue
    (std::string const& key, std::string_view const& baseUnit) const
    { return getValue<T, UnitT>(*fMetadata, key, baseUnit); }
  
  template <typename Vector>
  std::optional<Vector> getVector(std::string const& key) const
    { return getVector<Vector>(*fMetadata, key); }

  template <typename Vector, typename UnitT = double>
  std::optional<Vector> getVector
    (std::string const& key, std::string_view const& baseUnit) const
    { return getVector<Vector, UnitT>(*fMetadata, key, baseUnit); }
  
  auto operator()
    (std::string const& key, std::string_view const& baseUnit) const
    { return getValue(key, baseUnit); }
  
  auto operator() (std::string const& key) const { return getValue(key); }
  
  
  // --- BEGIN --- String manipulation -----------------------------------------
  /// @name String manipulation
  /// @{
  
  /// Returns `s` after depriving it of all the blanks on its left.
  static std::string_view& trimLeft(std::string_view& s);

  /// Returns `s` after depriving it of all the blanks on its right.
  static std::string_view& trimRight(std::string_view& s);
  
  /// Returns `s` after depriving it of all the blanks on left and right.
  static std::string_view& trim(std::string_view& s);

  /// Splits `s` in everything before the first `'*'`, and stuff after (if any).
  std::pair<std::string_view, std::string_view>
  static splitValueAndUnit(std::string_view const& s);
  
  /**
   * @brief Returns the factor to change a `unit` into its unprefixed version.
   * @tparam T (default: `double`) type of number to return the factor in
   * @param unit the unit to be converted
   * @param baseUnit the base unit
   * @return the factor to multiply to a `unit` value to convert to `baseUnit`
   * @throw cet::exception if `unit` does not derive from `baseUnit`
   * @throw cet::exception if the prefix in `unit` is not known or supported
   * 
   * For example, `unitFactor("cm", "m")` returns `0.01`.
   * Note that `baseUnit` can't have prefixes (e.g. `"mm"` won't do).
   */
  template <typename T = double>
  static T unitFactor
    (std::string_view const& unit, std::string_view const& baseUnit);
  
  /**
   * @brief Converts a string into a value (using `std::from_char()`).
   * @tparam T (default: `double`) type of the value to be returned
   * @param s string to be converted
   * @return the value of type `T` extracted from `s`
   * @throw cet::exception (category: `"GeoMetadataParser"`)
   *                       conversion failed or content is left unconverted
   * 
   * This function uses `std::from_char()` to perform the conversion, and
   * throws exceptions on errors.
   */
  template <typename T = double>
  static double parseValue(std::string_view const& s);
  
  /// @}
  // --- END --- String manipulation -------------------------------------------
  
  
  // --- BEGIN --- Parsing static functions ------------------------------------
  /// @name Parsing static functions
  /// @{
  
  /**
   * @brief Retrieves a value of type `T` from metadata element named `key`.
   * @tparam T (default: `double`) type of value to read
   * @param metadata the metadata map to extract metadata from
   * @param key the name of the metadata element to parse
   * @return the value of the metadata element, or no value if no element `key`
   * @throw cet::exception propagated from `parseValue()`
   * 
   * The conversion of the element is done via `parseValue()`.
   * The metadata element must *not* have a unit (no unit is stripped from the
   * value).
   */
  template <typename T = double>
  static std::optional<T> getValue
    (TMap const& metadata, std::string const& key);
  
  /**
   * @brief Retrieves a value of type `T` from metadata element named `key`.
   * @tparam T (default: `double`) type of value to read
   * @tparam UnitT (default: `double`) type of value for the unit scale
   * @param metadata the metadata map to extract metadata from
   * @param key the name of the metadata element to parse
   * @param baseUnit base unit of measurement for the value
   * @return the value of the metadata element in `baseUnit` units
   * @throw cet::exception propagated from `parseValue()`
   * @throw cet::exception (category `"GeoMetadataParser"`) if no unit is
   *                       present in the metadata
   * 
   * In addition to the behaviour of `getValue(TMap const&, std::string const&)`
   * this version of `getValue()` requires that the metadata element has a unit
   * and that unit must be in the form `SIprefix + baseUnit` (e.g. `"mm"`).
   * The base unit needs to be with no prefix (so, use `"g"` rather than `kg`),
   * and _the value is returned in `baseUnit` units_, not in the units specified
   * in the metadata. To be explicit: with
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * <auxiliary auxtype="length" auxvalue="30" auxunit="cm" />
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * the call `getValue(metadata, "length", "m")` will return `0.3`, meaning
   * meters as indicated by the `baseUnit` parameter.
   * 
   * The `UnitT` type is internally used to store the conversion factor.
   * It may need to be set to an integral type in particular cases (e.g. for
   * data size like "123 kB" when forcing integer numbers.
   */
  template <typename T = double, typename UnitT = double>
  static std::optional<T> getValue(
    TMap const& metadata, std::string const& key,
    std::string_view const& baseUnit
    );
  
  /**
   * @brief Retrieves a vector of type `Vector` from metadata element `key`.
   * @tparam Vector type of vector to read
   * @param metadata the metadata map to extract metadata from
   * @param key the name of the metadata element to parse
   * @return the value of the metadata element, or no value if no element `key`
   * @throw std::runtime_error propagated from the vector parser
   * 
   * The conversion of the element is done via `geo::vect::parse::VectorParser`,
   * and the supported vector formats are the same.
   * In doubt: something like `"( 1e-3, 6.0, 8 )"` will work (with semicolons
   * instead of commas, too).
   * 
   * The `Vector` type is any type which supports the coordinate managers
   * (see `geo::vect::coords()`). Types known to work are `geo::Vector_t`,
   * `geo::Point_t`, all the ROOT GenVectors vectors, `TVector3`, etc.
   * The dimensions of the vectors are limited to 2 to 4 because of the
   * coordinate manager limitations.
   * 
   * The metadata element must *not* have a unit (no unit is stripped from the
   * value).
   */
  template <typename Vector>
  static std::optional<Vector> getVector
    (TMap const& metadata, std::string const& key);

  /**
   * @brief Retrieves a vector of type `Vector` from metadata element `key`.
   * @tparam Vector type of vector to read
   * @tparam UnitT (default: `double`) type of value for the unit scale
   * @param metadata the metadata map to extract metadata from
   * @param key the name of the metadata element to parse
   * @param baseUnit base unit of measurement for the value
   * @return the value of the metadata element in `baseUnit` units
   * @throw std::runtime_error propagated from the vector parser
   * @throw cet::exception (category `"GeoMetadataParser"`) if no unit is
   *                       present in the metadata
   * 
   * This methods adds to `getVector(TMap const&, std::string const&)` the unit
   * management features also found in
   * `getValue(TMap const&, std::string const&, std::string_view const&)`.
   * 
   * As that method, this also _requires_ the presence of a unit specification
   * in the metadata, and it applies the unit conversion factor on each of the
   * components of the vector.
   * 
   * The `Vector` type is required to support the product to a scalar, in the
   * form of `Vector operator* (Vector const&, UnitT const&)`.
   */
  template <typename Vector, typename UnitT = double>
  static std::optional<Vector> getVector(
    TMap const& metadata, std::string const& key,
    std::string_view const& baseUnit
    );
  
  /// @}
  // --- END --- Parsing static functions --------------------------------------
  
  
}; // class geo::GeoMetadataParser


// -----------------------------------------------------------------------------
// ---  template implementation
// -----------------------------------------------------------------------------
template <typename T /* = double */>
T geo::GeoMetadataParser::unitFactor
  (std::string_view const& unit, std::string_view const& baseUnit)
{
  if (!util::ends_with(unit, baseUnit)) {
    throw cet::exception("GeoMetadataParser")
      << "unitFactor(): unit '" << unit << "' does not derive from '"
      << baseUnit << "'\n"
      ;
  } // unit check
  
  auto prefix { unit };
  prefix.remove_suffix(baseUnit.length());
  if (prefix.empty())      return T{ 1.0   };
  else if (prefix == "a")  return T{ 1e-18 };
  else if (prefix == "f")  return T{ 1e-15 };
  else if (prefix == "p")  return T{ 1e-12 };
  else if (prefix == "n")  return T{ 1e-09 };
  else if (prefix == "u")  return T{ 1e-06 };
  else if (prefix == "m")  return T{ 1e-03 };
  else if (prefix == "c")  return T{ 1e-02 };
  else if (prefix == "d")  return T{ 1e-01 };
  else if (prefix == "da") return T{ 1e+01 };
  else if (prefix == "h")  return T{ 1e+02 };
  else if (prefix == "k")  return T{ 1e+03 };
  else if (prefix == "M")  return T{ 1e+06 };
  else if (prefix == "G")  return T{ 1e+09 };
  else if (prefix == "T")  return T{ 1e+12 };
  else if (prefix == "P")  return T{ 1e+15 };
  else if (prefix == "E")  return T{ 1e+18 };
  else {
    throw cet::exception("GeoMetadataParser")
      << "unitFactor('" << unit << "', '" << baseUnit << "'): prefix '"
      << prefix << "' unknown.\n";
  }
  
} // geo::GeoMetadataParser::unitFactor()


// -----------------------------------------------------------------------------
template <typename T /* = double */>
double geo::GeoMetadataParser::parseValue(std::string_view const& s) {
  
  T value = std::numeric_limits<T>::has_signaling_NaN
    ? std::numeric_limits<T>::signaling_NaN()
    : std::numeric_limits<T>::max()
    ;
  auto const begin = s.data();
  auto const end = begin + s.length();
  
  if (
    auto [ ptr, ec ] = util::pre_std::from_chars(begin, end, value);
    ec != std::errc{}
    )
  {
    throw cet::exception("GeoMetadataCollector")
      << "Error parsing the value '" << s << "': "
      << std::make_error_condition(ec).message() << "\n"
      ;
  }
  else {
    auto cursor = ptr;
    while (cursor != end) if (!std::isblank(*cursor++)) break;
    if (cursor != end) {
      throw cet::exception("GeoMetadataCollector")
        << "Error parsing the value '" << s
        << "': extra characters after having parsed value '" << value << "' ('"
        << std::string_view(ptr, end - ptr) << "')\n"
        ;
    } // if extra characters
  } // if ... else
  
  return value;
} // geo::GeoMetadataParser::parseValue()


// -----------------------------------------------------------------------------
template <typename T /* = double */>
std::optional<T> geo::GeoMetadataParser::getValue
  (TMap const& metadata, std::string const& key)
{
  using Data_t = T;
  auto pValue = getMetadataString(metadata, key);
  return pValue
    ? std::make_optional(parseValue<Data_t>(pValue->GetString().View()))
    : std::nullopt
    ;
} // geo::GeoMetadataParser::getValue()


// -----------------------------------------------------------------------------
template <typename T /* = double */, typename UnitT /* = double */>
std::optional<T> geo::GeoMetadataParser::getValue(
  TMap const& metadata, std::string const& key,
  std::string_view const& baseUnit
  )
{
  using Data_t = T;
  using UnitFactor_t = UnitT;
  
  auto pValue = getMetadataString(metadata, key);
  if (!pValue) return {};
  
  auto [ valueStr, unitStr ] = splitValueAndUnit(pValue->GetString().View());
  if (unitStr.empty()) {
    throw cet::exception("GeoMetadataParser")
      << "Error parsing metadata '" << key << "' (='" << valueStr
      << "'): unit specification not found, expected one from '" << baseUnit
      << "'.\n"
      ;
  } // if no unit found
  
  auto const value = parseValue<Data_t>(valueStr);
  auto const unitConv = unitFactor<UnitFactor_t>(trim(unitStr), baseUnit);
  return { value * unitConv };
} // geo::GeoMetadataParser::getValue()


// -----------------------------------------------------------------------------
template <typename Vector>
std::optional<Vector> geo::GeoMetadataParser::getVector
  (TMap const& metadata, std::string const& key)
{
  using Vector_t = Vector;
  auto pValue = getMetadataString(metadata, key);
  return pValue
    ? std::make_optional
      (geo::vect::parse::parse<Vector_t>(pValue->GetString().View()))
    : std::nullopt
    ;
} // geo::GeoMetadataParser::getVector()


// -----------------------------------------------------------------------------
template <typename Vector, typename UnitT /* = double */>
std::optional<Vector> geo::GeoMetadataParser::getVector(
  TMap const& metadata, std::string const& key,
  std::string_view const& baseUnit /* = "" */
  )
{
  using Vector_t = Vector;
  using UnitFactor_t = UnitT;
  
  auto pValue = getMetadataString(metadata, key);
  if (!pValue) return {};
  
  auto [ valueStr, unitStr ] = splitValueAndUnit(pValue->GetString().View());
  if (unitStr.empty()) {
    throw cet::exception("GeoMetadataParser")
      << "Error parsing metadata '" << key << "' (='" << valueStr
      << "'): unit specification not found, expected one from '" << baseUnit
      << "'.\n"
      ;
  } // if no unit found
  
  auto const value = geo::vect::parse::parse<Vector_t>(valueStr);
  auto const unitConv = unitFactor<UnitFactor_t>(trim(unitStr), baseUnit);
  return { value * unitConv };
} // geo::GeoMetadataParser::getVector()


// -----------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_GEOMETADATAPARSER_H
