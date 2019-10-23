////////////////////////////////////////////////////////////////////////
/// \file  AuxDetGeo.cxx
/// \brief Encapsulate the geometry of an auxilary detector
///
/// \author  miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/AuxDetGeo.h"

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect namespace
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/CoreUtils/SortByPointers.h"

// ROOT
#include "TGeoTrd2.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C/C++ libraries
#include <sstream> // std::ostringstream
#include <string>
#include <limits>

namespace geo{

  //-----------------------------------------
  AuxDetGeo::AuxDetGeo(
    TGeoNode const& node, geo::TransformationMatrix&& trans,
    AuxDetSensitiveList_t&& sensitive
    )
    : fTotalVolume(node.GetVolume())
    , fTrans(std::move(trans))
    , fSensitive(std::move(sensitive))
  {

    if(!fTotalVolume)
      throw cet::exception("AuxDetGeo") << "cannot find AuxDet volume\n";

    MF_LOG_DEBUG("Geometry") << "detector total  volume is " << fTotalVolume->GetName();

    // if there are no sensitive volumes then this aux det
    // could be from an older gdml file than the introduction of AuxDetSensitiveGeo
    // in that case assume the full AuxDetGeo is sensitive and copy its information
    // into a single AuxDetSensitive
    if (fSensitive.empty())
      fSensitive.emplace_back(node, geo::TransformationMatrix(fTrans.Matrix()));

    InitShapeSize();

  }

  //......................................................................
  geo::Point_t AuxDetGeo::GetCenter(double localz /* = 0.0 */) const
    { return toWorldCoords(LocalPoint_t{ 0.0, 0.0, localz }); }

  //......................................................................
  void AuxDetGeo::GetCenter(double* xyz, double localz) const {
    auto const& center = GetCenter(localz);
    xyz[0] = center.X();
    xyz[1] = center.Y();
    xyz[2] = center.Z();
  } // AuxDetGeo::GetCenter(double*)

  //......................................................................

  // Return the unit normal vector (0,0,1) in local coordinates to global coordinates
  geo::Vector_t AuxDetGeo::GetNormalVector() const
    { return toWorldCoords(geo::Zaxis<LocalVector_t>()); }

  //......................................................................

  // Return the unit normal vector (0,0,1) in local coordinates to global coordinates
  void AuxDetGeo::GetNormalVector(double* xyzDir) const {
    auto const& norm = GetNormalVector();
    xyzDir[0] = norm.X();
    xyzDir[1] = norm.Y();
    xyzDir[2] = norm.Z();
  } // AuxDetGeo::GetNormalVector(double*)


  //......................................................................
  geo::Length_t AuxDetGeo::DistanceToPoint(double const* point) const
    { return DistanceToPoint(geo::vect::makePointFromCoords(point)); }

  //......................................................................
  std::size_t AuxDetGeo::FindSensitiveVolume(geo::Point_t const& point) const
  {
    for(std::size_t a = 0; a < fSensitive.size(); ++a) {
      auto const& sensVol = SensitiveVolume(a);

      auto const local = sensVol.toLocalCoords(point);

      double const HalfCenterWidth = sensVol.HalfCenterWidth();

      double const deltaWidth
        = local.Z()*(HalfCenterWidth-sensVol.HalfWidth2())/sensVol.HalfLength();

      if (local.Z() >= - sensVol.HalfLength()       &&
          local.Z() <=   sensVol.HalfLength()       &&
          local.Y() >= - sensVol.HalfHeight()     &&
          local.Y() <=   sensVol.HalfHeight()     &&
          // if SensitiveVolume a is a box, then HalfSmallWidth = HalfWidth
          local.X() >= - HalfCenterWidth + deltaWidth &&
          local.X() <=   HalfCenterWidth - deltaWidth
        )  return a;

    }// for loop over AuxDetSensitive a

    throw cet::exception("AuxDetGeo")
      << "Can't find AuxDetSensitive for position " << point << "\n";
    // the following is not very useful right now...
    return std::numeric_limits<std::size_t>::max();
  } // AuxDetGeo::FindSensitiveVolume(geo::Point_t)

  //......................................................................
  std::size_t AuxDetGeo::FindSensitiveVolume(double const worldPos[3]) const
    { return FindSensitiveVolume(geo::vect::makePointFromCoords(worldPos)); }


  //......................................................................
  AuxDetSensitiveGeo const& AuxDetGeo::PositionToSensitiveVolume
    (geo::Point_t const& point, size_t& sv) const
  {
    sv = FindSensitiveVolume(point);
    if (sv == std::numeric_limits<std::size_t>::max()) {
      throw cet::exception("AuxDetGeo")
        << "Can't find AuxDetSensitiveGeo for position " << point << "\n";
    }
    return SensitiveVolume(sv);
  }

  //......................................................................
  AuxDetSensitiveGeo const& AuxDetGeo::PositionToSensitiveVolume
    (double const worldLoc[3], size_t& sv) const
    { return PositionToSensitiveVolume(geo::vect::makePointFromCoords(worldLoc), sv); }

  //......................................................................
  void AuxDetGeo::SortSubVolumes(GeoObjectSorter const& sorter)
  {
    // the sorter interface requires a vector of pointers;
    // sorting is faster, but some gymnastics is required:
    util::SortByPointers
      (fSensitive, [&sorter](auto& coll){ sorter.SortAuxDetSensitive(coll); });
  }

  //......................................................................
  std::string AuxDetGeo::AuxDetInfo
    (std::string indent /* = "" */, unsigned int verbosity /* = 1 */) const
  {
    std::ostringstream sstr;
    PrintAuxDetInfo(sstr, indent, verbosity);
    return sstr.str();
  } // AuxDetGeo::AuxDetInfo()


  //......................................................................
  void AuxDetGeo::InitShapeSize() {
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
  } // AuxDetGeo::InitShapeSize()

  //......................................................................

}
////////////////////////////////////////////////////////////////////////
