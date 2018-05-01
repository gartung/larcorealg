////////////////////////////////////////////////////////////////////////
/// \file larcorealg/Geometry/CryostatGeo.cxx
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/CryostatGeo.h"

// LArSoft includes
#include "larcorealg/CoreUtils/SortByPointers.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// ROOT includes
#include "TGeoNode.h"
#include "TGeoBBox.h"
#include "TClass.h"

// C++ standard libraries
#include <limits> // std::numeric_limits<>
#include <algorithm> // std::sort()


namespace geo{


  //......................................................................
  // Define sort order for detector tpcs.
  static bool opdet_sort(const OpDetGeo* t1, const OpDetGeo* t2) 
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    if(xyz1[2]!=xyz2[2])
      return xyz1[2]>xyz2[2];
    else if(xyz1[1]!=xyz2[1])
      return xyz1[1]>xyz2[1];
    else
      return xyz1[0]>xyz2[0];
  }

  // DUNE specific sorting originally intended for 10kt
  // executed when there are 600+ opdets. -talion
  ///\todo: move dune opdet sorting to appropriate place in dunetpc 
  static bool DUNE_opdet_sort(const OpDetGeo* t1, const OpDetGeo* t2)
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    if(xyz1[0]!=xyz2[0])
      return xyz1[0]>xyz2[0];
    else if(xyz1[2]!=xyz2[2])
      return xyz1[2]>xyz2[2];
    else
    return xyz1[1]>xyz2[1];
  }

  //......................................................................
  CryostatGeo::CryostatGeo(std::vector<const TGeoNode*>& path, int depth)
    : fTrans(path, depth)
    , fVolume(nullptr)
  {
    
    // all planes are going to be contained in the volume named volCryostat
    // now get the total volume of the Cryostat
    fVolume = path[depth]->GetVolume();
    if(!fVolume)
      throw cet::exception("CryostatGeo") << "cannot find cryostat outline volume\n";
    
    LOG_DEBUG("Geometry") << "cryostat  volume is " << fVolume->GetName();

    // set the bounding box
    InitCryoBoundaries();
    
    // find the tpcs for the cryostat so that you can use them later
    this->FindTPC(path, depth);


    // Set OpDetName;
    fOpDetGeoName = "volOpDetSensitive";
    
    // find the opdets for the cryostat so that you can use them later
    this->FindOpDet(path, depth);
    
    // sort the OpDets according to xyz position
    // 600 intended to separate dune10kt geometry from others when sorting
    ///\todo: remove the hard-coded 600 in favor of selecting sorting the same way as in ChannelMapAlgs
    /// (LArSoft issue #16812)
    auto sorter = (fOpDets.size() != 600)? opdet_sort: DUNE_opdet_sort;
    util::SortByPointers(fOpDets,
      [&sorter](auto& coll){ std::sort(coll.begin(), coll.end(), sorter); }
      );
    
  }


  //......................................................................
  void CryostatGeo::FindTPC(std::vector<const TGeoNode*>& path,
			    unsigned int depth) 
  {

    const char* nm = path[depth]->GetName();
    if( (strncmp(nm, "volTPC", 6) == 0) ){
      this->MakeTPC(path,depth);
      return;
    }

    //explore the next layer down
    unsigned int deeper = depth+1;
    if(deeper >= path.size()){
      throw cet::exception("BadTGeoNode") << "exceeded maximum TGeoNode depth\n";
    }

    const TGeoVolume *v = path[depth]->GetVolume();
    int nd = v->GetNdaughters();
    for(int i = 0; i < nd; ++i){
      path[deeper] = v->GetNode(i);
      this->FindTPC(path, deeper);
    }
  
  }

  //......................................................................
  void CryostatGeo::MakeTPC(std::vector<const TGeoNode*>& path, int depth) 
  {
    fTPCs.emplace_back(path, depth);
  }

  
  //......................................................................
  // sort the TPCGeo objects, and the PlaneGeo objects inside
  void CryostatGeo::SortSubVolumes(geo::GeoObjectSorter const& sorter)
  {
    //
    // TPCs
    //
    util::SortByPointers
      (fTPCs, [&sorter](auto& coll){ sorter.SortTPCs(coll); });
    
    for (geo::TPCGeo& TPC: fTPCs) { 
      TPC.SortSubVolumes(sorter);
    } // for TPCs
    
    //
    // optical detectors
    //
    
    // sorting of optical detectors happens elsewhere
    
  } // CryostatGeo::SortSubVolumes()


  //......................................................................
  void CryostatGeo::UpdateAfterSorting(geo::CryostatID cryoid) {
    
    // update the cryostat ID
    fID = cryoid;
    
    // trigger all the TPCs to update as well
    for (unsigned int tpc = 0; tpc < NTPC(); ++tpc)
      fTPCs[tpc].UpdateAfterSorting(geo::TPCID(fID, tpc));
    
  } // CryostatGeo::UpdateAfterSorting()
  
  
  //......................................................................
  const TPCGeo& CryostatGeo::TPC(unsigned int itpc) const
  {
    TPCGeo const* pTPC = TPCPtr(itpc);
    if(!pTPC){
      throw cet::exception("TPCOutOfRange") << "Request for non-existant TPC "
                                            << itpc << "\n";
    }

    return *pTPC;
  }



  //......................................................................
  void CryostatGeo::FindOpDet(std::vector<const TGeoNode*>& path,
			      unsigned int depth) 
  {

    const char* nm = path[depth]->GetName();
    if( (strncmp(nm, OpDetGeoName().c_str(), 6) == 0) ){
      this->MakeOpDet(path,depth);
      return;
    }

    //explore the next layer down
    unsigned int deeper = depth+1;
    if(deeper >= path.size()){
      throw cet::exception("BadTGeoNode") << "exceeded maximum TGeoNode depth\n";
    }

    const TGeoVolume *v = path[depth]->GetVolume();
    int nd = v->GetNdaughters();
    for(int i = 0; i < nd; ++i){
      path[deeper] = v->GetNode(i);
      this->FindOpDet(path, deeper);
    }
  
  }

  //......................................................................
  void CryostatGeo::MakeOpDet(std::vector<const TGeoNode*>& path, int depth) 
  {
    fOpDets.emplace_back(path, depth);
  }

  //......................................................................
  const OpDetGeo& CryostatGeo::OpDet(unsigned int iopdet) const
  {
    if(iopdet >= fOpDets.size()){
      throw cet::exception("OpDetOutOfRange") << "Request for non-existant OpDet " 
					      << iopdet;
    }

    return fOpDets[iopdet];
  }
  
  
  
  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the 
  // passed in world loc relative to the boundaries.
  geo::TPCID::TPCID_t CryostatGeo::FindTPCAtPosition
    (double const worldLoc[3], double wiggle) const
  {
    geo::TPCID tpcid
      = PositionToTPCID(geo::vect::makePointFromCoords(worldLoc), wiggle);
    return tpcid? tpcid.TPC: geo::TPCID::InvalidID;
  } // CryostatGeo::FindTPCAtPosition()

  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the 
  // passed in world loc relative to the boundaries.
  geo::TPCID CryostatGeo::PositionToTPCID
    (geo::Point_t const& point, double wiggle) const
  { 
    geo::TPCGeo const* tpc = PositionToTPCptr(point, wiggle);
    return tpc? tpc->ID(): geo::TPCID{};
  }

  //......................................................................
  // wiggle is 1+a small number to allow for rounding errors on the 
  // passed in world loc relative to the boundaries.
  TPCGeo const& CryostatGeo::PositionToTPC
    (geo::Point_t const& point, double wiggle) const
  {
    geo::TPCGeo const* tpc = PositionToTPCptr(point, wiggle);
    if (!tpc) {
      throw cet::exception("CryostatGeo")
        << "Can't find any TPC for position " << point << " within " << ID()
        << "\n";
    }
    return *tpc;
  }
  
  //......................................................................
  geo::TPCGeo const* CryostatGeo::PositionToTPCptr
    (geo::Point_t const& point, double wiggle) const
  {
    for (auto const& tpc: TPCs())
      if (tpc.ContainsPosition(point, wiggle)) return &tpc;
    return nullptr;
  } // CryostatGeo::PositionToTPCptr()
  
  
  //......................................................................
  unsigned int CryostatGeo::MaxPlanes() const {
    unsigned int maxPlanes = 0;
    for (geo::TPCGeo const& TPC: fTPCs) {
      unsigned int maxPlanesInTPC = TPC.Nplanes();
      if (maxPlanesInTPC > maxPlanes) maxPlanes = maxPlanesInTPC;
    } // for
    return maxPlanes;
  } // CryostatGeo::MaxPlanes()
  
  //......................................................................
  unsigned int CryostatGeo::MaxWires() const {
    unsigned int maxWires = 0;
    for (geo::TPCGeo const& TPC: fTPCs) {
      unsigned int maxWiresInTPC = TPC.MaxWires();
      if (maxWiresInTPC > maxWires) maxWires = maxWiresInTPC;
    } // for
    return maxWires;
  } // CryostatGeo::MaxWires()
  
  //......................................................................
  double CryostatGeo::HalfWidth()  const 
  {
    return static_cast<TGeoBBox const*>(fVolume->GetShape())->GetDX();
  }

  //......................................................................
  double CryostatGeo::HalfHeight() const 
  {
    return static_cast<TGeoBBox const*>(fVolume->GetShape())->GetDY();
  }

  //......................................................................
  double CryostatGeo::HalfLength() const
  { 
    return static_cast<TGeoBBox const*>(fVolume->GetShape())->GetDZ();
  }

  //......................................................................
  void CryostatGeo::Boundaries(double* boundaries) const {
    boundaries[0] = MinX();
    boundaries[1] = MaxX();
    boundaries[2] = MinY();
    boundaries[3] = MaxY();
    boundaries[4] = MinZ();
    boundaries[5] = MaxZ();
  } // CryostatGeo::CryostatBoundaries(double*)
  
  
  //......................................................................
  // Find the nearest opdet to point in this cryostat

  geo::OpDetGeo const* CryostatGeo::GetClosestOpDetPtr
    (geo::Point_t const& point) const
  {
    unsigned int iOpDet = GetClosestOpDet(point);
    return
      (iOpDet == std::numeric_limits<double>::max())? nullptr: &OpDet(iOpDet);
  }
  
  //......................................................................
  unsigned int CryostatGeo::GetClosestOpDet(geo::Point_t const& point) const {
    unsigned int ClosestDet = std::numeric_limits<unsigned int>::max();
    double ClosestDist = std::numeric_limits<double>::max();
    
    for(unsigned int o = 0U; o < NOpDet(); ++o) {
      double const ThisDist = OpDet(o).DistanceToPoint(point); 
      if(ThisDist < ClosestDist) {
        ClosestDist = ThisDist;
        ClosestDet  = o;
      }
    } // for
    return ClosestDet;
  } // CryostatGeo::GetClosestOpDet(geo::Point_t)
  
  //......................................................................
  unsigned int CryostatGeo::GetClosestOpDet(double const* point) const
    { return GetClosestOpDet(geo::vect::makePointFromCoords(point)); }
  
  //......................................................................
  void CryostatGeo::InitCryoBoundaries() {
    
    // check that this is indeed a box
    if (!dynamic_cast<TGeoBBox*>(Volume()->GetShape())) {
      // at initialisation time we don't know yet our real ID
      throw cet::exception("CryostatGeo") << "Cryostat is not a box! (it is a "
        << Volume()->GetShape()->IsA()->GetName() << ")\n";
    }
    
    // get the half width, height, etc of the cryostat
    const double halflength = HalfLength();
    const double halfwidth  = HalfWidth();
    const double halfheight = HalfHeight();
    
    SetBoundaries(
      toWorldCoords(LocalPoint_t{ -halfwidth, -halfheight, -halflength }),
      toWorldCoords(LocalPoint_t{ +halfwidth, +halfheight, +halflength })
      );
    
  } // CryostatGeo::InitCryoBoundaries()
  
  
  //......................................................................

} // namespace geo
////////////////////////////////////////////////////////////////////////
