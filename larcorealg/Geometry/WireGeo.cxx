/**
 * @file    larcorealg/Geometry/WireGeo.cxx
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    April 26, 2019
 * @see     larcorealg/Geometry/WireGeo.h
 * @ingroup Geometry
 * 
 * Nothing here, except a help to make sure the header correctly compiles.
 */

// class header
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"

// framework libraries
#include "cetlib_except/exception.h"

// C/C++ standard libraries
#include <sstream> // std::ostringstream


//------------------------------------------------------------------------------
geo::WireGeo::WireGeo(geo::PlaneGeo const& plane, geo::WireID const& wireid)
  : WireGeo(plane, wireid.Wire)
{
  if (wireid.asPlaneID() != plane.ID()) {
    throw cet::exception("Geometry")
      << "Attempt to create a geo::WireGeo with ID " << wireid
      << " on the plane " << plane.ID() << "\n";
  }
} // geo::WireGeo::WireGeo()


//------------------------------------------------------------------------------
geo::WireID geo::WireGeo::ID() const
  { return { plane().ID(), location().wireNo }; }


//------------------------------------------------------------------------------
bool geo::WireGeo::IsValid() const
  { return fData.plane && plane().HasWire(location().wireNo); }


//------------------------------------------------------------------------------
double geo::WireGeo::RMax() const 
  { return askThePlane(&geo::PlaneGeo::doWireRMax); }


//------------------------------------------------------------------------------
double geo::WireGeo::RMin() const
  { return askThePlane(&geo::PlaneGeo::doWireRMin); }


//------------------------------------------------------------------------------
double geo::WireGeo::HalfL() const
  { return askThePlane(&geo::PlaneGeo::doWireHalfL); }


//------------------------------------------------------------------------------
void geo::WireGeo::GetCenter(double* xyz, double localz /* = 0.0 */) const
  { askThePlane(&geo::PlaneGeo::doWireFillCenterXYZ, xyz, localz); }


//------------------------------------------------------------------------------
void geo::WireGeo::GetStart(double* xyz) const
  { askThePlane(&geo::PlaneGeo::doWireFillStartXYZ, xyz); }


//------------------------------------------------------------------------------
void geo::WireGeo::GetEnd(double* xyz) const
  { askThePlane(&geo::PlaneGeo::doWireFillEndXYZ, xyz); }


//------------------------------------------------------------------------------
geo::Point_t geo::WireGeo::GetPositionFromCenterImpl(double localz) const
  { return askThePlane(&geo::PlaneGeo::doWireGetPositionFromCenter, localz); }

  
//------------------------------------------------------------------------------
geo::Point_t geo::WireGeo::GetPositionFromCenterUnboundedImpl
  (double localz) const
{
  return 
    askThePlane(&geo::PlaneGeo::doWireGetPositionFromCenterUnbounded, localz);
} // geo::WireGeo::GetPositionFromCenterUnboundedImpl()


//------------------------------------------------------------------------------
geo::Point_t geo::WireGeo::GetCenterImpl() const
  { return askThePlane(&geo::PlaneGeo::doWireGetCenter); }

  
//------------------------------------------------------------------------------
geo::Point_t geo::WireGeo::GetStartImpl() const
  { return askThePlane(&geo::PlaneGeo::doWireGetStart); }

  
//------------------------------------------------------------------------------
geo::Point_t geo::WireGeo::GetEndImpl() const
  { return askThePlane(&geo::PlaneGeo::doWireGetEnd); }

  
//------------------------------------------------------------------------------
double geo::WireGeo::Length() const
  { return askThePlane(&geo::PlaneGeo::doWireLength); }

  
//------------------------------------------------------------------------------
double geo::WireGeo::ComputeZatY0() const
  { return askThePlane(&geo::PlaneGeo::doWireComputeZatY0); }


//------------------------------------------------------------------------------
double geo::WireGeo::DistanceFrom(WireGeo const& wire) const
  { return askThePlane(&geo::PlaneGeo::doWireDistanceFrom, wire); }


//------------------------------------------------------------------------------
double geo::WireGeo::ThetaZ() const
  { return askThePlane(&geo::PlaneGeo::doWireThetaZ); }


//------------------------------------------------------------------------------
double geo::WireGeo::CosThetaZ() const
  { return askThePlane(&geo::PlaneGeo::doWireCosThetaZ); }


//------------------------------------------------------------------------------
double geo::WireGeo::SinThetaZ() const
  { return askThePlane(&geo::PlaneGeo::doWireSinThetaZ); }


//------------------------------------------------------------------------------
double geo::WireGeo::TanThetaZ() const
  { return askThePlane(&geo::PlaneGeo::doWireTanThetaZ); }


//------------------------------------------------------------------------------
bool geo::WireGeo::isParallelTo(WireGeo const& wire) const
  { return askThePlane(&geo::PlaneGeo::doWireIsParallelTo, wire); }


//------------------------------------------------------------------------------
geo::Vector_t geo::WireGeo::DirectionImpl() const
{ return askThePlane(&geo::PlaneGeo::doWireDirection); }


//------------------------------------------------------------------------------
TGeoNode const* geo::WireGeo::Node() const
  { return askThePlane(&geo::PlaneGeo::doWireNode); }


//------------------------------------------------------------------------------
std::string geo::WireGeo::WireInfo
  (std::string indent /* = "" */, unsigned int verbosity /* = 1 */) const
  { return askThePlane(&geo::PlaneGeo::doWireInfo, indent, verbosity); }


//------------------------------------------------------------------------------
std::string geo::WireGeo::GenericWireInfo
  (std::string indent /* = "" */, unsigned int verbosity /* = 1 */) const
{
  std::ostringstream out;
  PrintGenericWireInfo(out, indent, verbosity);
  return out.str();
} // geo::WireGeo::GenericWireInfo()


//------------------------------------------------------------------------------
geo::WireGeo::LocalTransformation_t const* geo::WireGeo::trans() const
  { return askThePlane(&geo::PlaneGeo::doWireTrans); }


//------------------------------------------------------------------------------

