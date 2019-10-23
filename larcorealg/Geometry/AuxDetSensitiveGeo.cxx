////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/AuxDetSensitiveGeo.cxx
/// \brief Encapsulate the geometry of the sensitive portion of an auxilary detector
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect namespace

/// ROOT libraries
#include "TGeoTrd2.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <sstream> // std::ostringstream


namespace geo{

  //-----------------------------------------
  AuxDetSensitiveGeo::AuxDetSensitiveGeo
    (TGeoNode const& node, geo::TransformationMatrix&& trans)
    : fTrans(std::move(trans))
    , fTotalVolume(node.GetVolume())
  {

    MF_LOG_DEBUG("Geometry") << "detector sensitive total  volume is " << fTotalVolume->GetName();

    InitShapeSize();

  }

  //......................................................................
  geo::Point_t AuxDetSensitiveGeo::GetCenter(double localz /* = 0.0 */) const
    { return toWorldCoords(LocalPoint_t{ 0.0, 0.0, localz }); }

  //......................................................................
  void AuxDetSensitiveGeo::GetCenter(double* xyz, double localz) const {
    auto const& center = GetCenter(localz);
    xyz[0] = center.X();
    xyz[1] = center.Y();
    xyz[2] = center.Z();
  } // AuxDetSensitiveGeo::GetCenter(double*)

  //......................................................................

  // Return the unit normal vector (0,0,1) in local coordinates to global coordinates
  geo::Vector_t AuxDetSensitiveGeo::GetNormalVector() const
    { return toWorldCoords(geo::Zaxis<LocalVector_t>()); }

  //......................................................................

  // Return the unit normal vector (0,0,1) in local coordinates to global coordinates
  void AuxDetSensitiveGeo::GetNormalVector(double* xyzDir) const {
    auto const& norm = GetNormalVector();
    xyzDir[0] = norm.X();
    xyzDir[1] = norm.Y();
    xyzDir[2] = norm.Z();
  } // AuxDetSensitiveGeo::GetNormalVector(double*)


  //......................................................................
  geo::Length_t AuxDetSensitiveGeo::DistanceToPoint(double const* point) const
    { return DistanceToPoint(geo::vect::makePointFromCoords(point)); }


  //......................................................................
  std::string AuxDetSensitiveGeo::AuxDetInfo
    (std::string indent /* = "" */, unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintAuxDetInfo(sstr, indent, verbosity);
    return sstr.str();
  } // AuxDetSensitiveGeo::AuxDetInfo()


  //......................................................................
  void AuxDetSensitiveGeo::InitShapeSize() {
    // set the ends depending on whether the shape is a box or trapezoid
    std::string volName(fTotalVolume->GetName());
    if( volName.find("Trap") != std::string::npos ) {

      //       Small Width
      //          ____          Height is the thickness
      //         /    \     T     of the trapezoid
      //        /      \    |
      //       /        \   | Length
      //      /__________\  _
      //         Width
      fHalfHeight      =     ((TGeoTrd2*)fTotalVolume->GetShape())->GetDy1(); // same as Dy2()
      fLength          = 2.0*((TGeoTrd2*)fTotalVolume->GetShape())->GetDz();
      fHalfWidth1      =     ((TGeoTrd2*)fTotalVolume->GetShape())->GetDx1(); // at -Dz
      fHalfWidth2      =     ((TGeoTrd2*)fTotalVolume->GetShape())->GetDx2(); // at +Dz
    }
    else {
      fHalfWidth1      =     ((TGeoBBox*)fTotalVolume->GetShape())->GetDX();
      fHalfHeight      =     ((TGeoBBox*)fTotalVolume->GetShape())->GetDY();
      fLength          = 2.0*((TGeoBBox*)fTotalVolume->GetShape())->GetDZ();
      fHalfWidth2      = fHalfWidth1;
    }
  } // AuxDetSensitiveGeo::InitShapeSize()
}
////////////////////////////////////////////////////////////////////////
