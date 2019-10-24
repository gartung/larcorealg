/**
 * @file   larcorealg/Geometry/GeometryBuilderWireless.cxx
 * @brief  Implementation of wireless geometry extractor (implementation file).
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   February 26, 2019
 * @see    `larcorealg/Geometry/GeometryBuilderWireless.h`
 */

// LArSoft libraries
#include "larcorealg/Geometry/GeometryBuilderWireless.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"


//------------------------------------------------------------------------------
geo::GeometryBuilderWireless::GeometryBuilderWireless
  (Config const& config)
  : geo::GeometryBuilderStandard(config)
{
  MF_LOG_TRACE("GeometryBuilder")
    << "Loading geometry builder: GeometryBuilderWireless";
} // geo::GeometryBuilderWireless::GeometryBuilderWireless()


//------------------------------------------------------------------------------

