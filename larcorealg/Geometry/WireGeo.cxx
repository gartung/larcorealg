/**
 * @file    larcorealg/Geometry/WireGeo.cxx
 * @brief   Interface for a active readout element on a TPC plane.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    April 24, 2019
 * @see     larcorealg/Geometry/WireGeo.h
 * @ingroup Geometry
 */

// class header
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/CoreUtils/DebugUtils.h" // lar::debug namespace

// Framework includes
#include "cetlib_except/exception.h"

// C/C++ libraries
#include <sstream>


//------------------------------------------------------------------------------
std::string geo::WireGeo::GenericWireInfo
  (std::string indent /* = "" */, unsigned int verbosity /* = 1 */) const
{
  std::ostringstream out;
  
  //----------------------------------------------------------------------------
  out << "active plane element from " << GetStart<geo::Point_t>()
    << " to " << GetEnd<geo::Point_t>();

  if (verbosity-- <= 0) return out.str(); // 0

  //----------------------------------------------------------------------------
  out << " (" << Length() << " cm long)";

  if (verbosity-- <= 0) return out.str(); // 1

  //----------------------------------------------------------------------------
  out << ", theta(z)=" << ThetaZ() << " rad";

  if (verbosity-- <= 0) return out.str(); // 2

  //----------------------------------------------------------------------------
  out << "\n" << indent
    << "  center at " << GetCenter<geo::Point_t>() << " cm";

  if (verbosity-- <= 0) return out.str(); // 3

  //----------------------------------------------------------------------------
  out << ", direction: " << Direction<geo::Vector_t>();

//  if (verbosity-- <= 0) return out.str(); // 4

  //----------------------------------------------------------------------------
  
  return out.str();
  
} // geo::WireGeo::GenericWireInfo()


//------------------------------------------------------------------------------
[[noreturn]] void geo::WireGeo::NotImplemented() const {
  
  cet::exception e("NotImplemented");
  e << lar::debug::demangle(this)
    << " (from geo::WireGeo): call not implemented:\n";
  
  lar::debug::BacktracePrintOptions opts;
  opts.maxLines = 5U; // we keep it short...ish
  opts.skipLines = 2U; // skip `printBacktrace()` and `NotImplemented()` calls
  opts.setUniformIndent("  ");
  lar::debug::printBacktrace(e, opts);
  
  throw e;
  
} // geo::WireGeo::NotImplemented()


//------------------------------------------------------------------------------
