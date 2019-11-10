/**
 * @file    larcorealg/Geometry/GeoMetadataParser.cxx
 * @brief   Helper to parse metadata from ROOT geometry hierarchy.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    November 10, 2019
 * @see     `larcorealg/Geometry/GeoMetadataParser.h`
*/


// library header
#include "larcorealg/Geometry/GeoMetadataParser.h"

// C/C++ standard library
#include <cctype> // std::isblank()


// -----------------------------------------------------------------------------
TObjString const* geo::GeoMetadataParser::getMetadataString
  (TMap const& metadata, std::string const& key)
  { return dynamic_cast<TObjString const*>(metadata.GetValue(key.c_str())); }


// -----------------------------------------------------------------------------
std::string_view& geo::GeoMetadataParser::trimLeft(std::string_view& s) {
  auto const begin = s.cbegin();
  auto const end = s.cend();
  if (begin != end) {
    auto c = begin;
    do {
      if (!std::isblank(*c)) break;
    } while (++c != end);
    s.remove_prefix(c - begin);
  } // if
  return s;
} // geo::GeoMetadataParser::trimLeft()


// -----------------------------------------------------------------------------
std::string_view& geo::GeoMetadataParser::trimRight(std::string_view& s) {
  auto const begin = s.crbegin();
  auto const end = s.crend();
  if (begin != end) {
    auto c = begin;
    do {
      if (!std::isblank(*c)) break;
    } while (++c != end);
    s.remove_suffix(c - begin);
  } // if
  return s;
} // geo::GeoMetadataParser::trimRight()


// -----------------------------------------------------------------------------
std::string_view& geo::GeoMetadataParser::trim(std::string_view& s) {
  return trimLeft(trimRight(s));
} // geo::GeoMetadataParser::trim()


// -----------------------------------------------------------------------------
std::pair<std::string_view, std::string_view>
geo::GeoMetadataParser::splitValueAndUnit(std::string_view const& s) {
  using namespace std::string_view_literals;
  auto const iSep = s.find('*');
  return (iSep == s.npos)
    ? std::pair{ s, ""sv }
    : std::pair{ s.substr(0U, iSep), s.substr(iSep + 1U) }
    ;
} // geo::GeoMetadataParser::splitValueAndUnit()


// -----------------------------------------------------------------------------

