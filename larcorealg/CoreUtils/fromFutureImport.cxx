/**
 * @file   larcorealg/CoreUtils/fromFutureImport.cxx
 * @brief  Code that might appear as standard C++ in the future.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   November 8, 2019
 * @see    `larcorealg/CoreUtils/fromFutureImport.h`
 * 
 * This is currently a header-only library.
 */

// library header
#include "larcorealg/CoreUtils/fromFutureImport.h"

// C/C++ standard library
#include <limits>
#include <cmath> // HUGE_VAL
#include <cerrno>
#include <cstdlib> // std::strtod(), ...


// -----------------------------------------------------------------------------
util::pre_std::from_chars_result util::pre_std::from_chars
  (const char* first, const char* last, float& value)
{
  util::pre_std::from_chars_result res;
  char* str_err = nullptr; // can't send `res.ptr` directly because it's const
  auto const temp = std::strtof(first, &str_err);
  res.ptr = str_err;
  if (res.ptr == first)       res.ec = std::errc::invalid_argument;
  else if (temp == HUGE_VALF) res.ec = std::errc::result_out_of_range;
  else                        value = temp;
  
  return res;
} // pre_std::from_chars(double)


// -----------------------------------------------------------------------------
util::pre_std::from_chars_result util::pre_std::from_chars
  (const char* first, const char* last, double& value)
{
  util::pre_std::from_chars_result res;
  char* str_err = nullptr; // can't send `res.ptr` directly because it's const
  auto const temp = std::strtod(first, &str_err);
  res.ptr = str_err;
  if (res.ptr == first)      res.ec = std::errc::invalid_argument;
  else if (temp == HUGE_VAL) res.ec = std::errc::result_out_of_range;
  else                       value = temp;
  
  return res;
} // pre_std::from_chars(double)


// -----------------------------------------------------------------------------
util::pre_std::from_chars_result util::pre_std::from_chars
  (const char* first, const char* last, long double& value)
{
  util::pre_std::from_chars_result res;
  char* str_err = nullptr; // can't send `res.ptr` directly because it's const
  auto const temp = std::strtold(first, &str_err);
  res.ptr = str_err;
  if (res.ptr == first)       res.ec = std::errc::invalid_argument;
  else if (temp == HUGE_VALL) res.ec = std::errc::result_out_of_range;
  else                        value = temp;
  
  return res;
} // pre_std::from_chars(long double)


// -----------------------------------------------------------------------------
util::pre_std::from_chars_result util::pre_std::from_chars
  (const char* first, const char* last, int& value, int base /* = 10 */)
{
  util::pre_std::from_chars_result res;
  char* str_err = nullptr; // can't send `res.ptr` directly because it's const
  auto const temp = std::strtol(first, &str_err, base);
  res.ptr = str_err;
  if (res.ptr == first)     res.ec = std::errc::invalid_argument;
  else if (errno == ERANGE)
    { res.ec = std::errc::result_out_of_range; errno = 0; }
  else                      value = temp;
  
  return res;
} // pre_std::from_chars(int)


// -----------------------------------------------------------------------------
util::pre_std::from_chars_result util::pre_std::from_chars
  (const char* first, const char* last, unsigned int& value, int base /* = 10 */)
{
  util::pre_std::from_chars_result res;
  char* str_err = nullptr; // can't send `res.ptr` directly because it's const
  auto const temp = std::strtoul(first, &str_err, base);
  res.ptr = str_err;
  if (res.ptr == first)     res.ec = std::errc::invalid_argument;
  else if (errno == ERANGE)
    { res.ec = std::errc::result_out_of_range; errno = 0; }
  else                      value = temp;
  
  return res;
} // pre_std::from_chars(unsigned int)


// -----------------------------------------------------------------------------

