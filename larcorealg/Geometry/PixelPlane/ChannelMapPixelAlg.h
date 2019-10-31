/**
 * @file   larcorealg/Geometry/PixelPlane/ChannelMapPixelAlg.h
 * @brief  Simple channel mapping for pixel readout.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   October 24, 2019
 * @see    `larcorealg/Geometry/PixelPlane/ChannelMapPixelAlg.cxx`
 */


#ifndef LARCOREALG_GEOMETRY_PIXELPLANE_CHANNELMAPPIXELALG_H
#define LARCOREALG_GEOMETRY_PIXELPLANE_CHANNELMAPPIXELALG_H


// LArSoft libraries
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/GeometryDataContainers.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h" // readout::TPCsetID, ...
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// support libraries
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Comment.h"

// C/C++ standard libraries
#include <vector>
#include <set>
#include <optional>
#include <algorithm> // std::lower_bound()
#include <utility> // std::pair


//------------------------------------------------------------------------------
namespace geo {
  class GeoObjectSorter;
  class ChannelMapPixelAlg;
} // namespace geo

/**
 * @brief A simple channel mapping for "3D" TPC readout detectors.
 * 
 * This channel mapping is as simple as the standard one
 * (`geo::ChannelMapStandardAlg`), but even simpler.
 * 
 * The mapping of TPC sets to TPC's is one to one, as it is the mapping of
 * readout planes with sensitive planes. Also the mapping between readout
 * channels and sensitive elements (pixels) is one to one.
 * 
 * The channel numbering is a sequence following the same order as the
 * `geo::WireID` associated with each sensitive element, starting from `0`.
 * 
 */
class geo::ChannelMapPixelAlg: public geo::ChannelMapAlg {
  
    public:
  
  /// Configuration of the algorithm.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    // no configuration so far
    
  }; // struct Config
  
  
  // --- BEGIN -- object state management --------------------------------------
  
  /// Constructor: reads the configuration.
  ChannelMapPixelAlg(Config const& config);
  
  /// Extracts and caches necessary information from geometry structure.
  virtual void Initialize(geo::GeometryData_t const& geodata) override;
  
  /// Deconfiguration: prepare for a following call of Initialize().
  virtual void Uninitialize() override;
  
  // --- END -- object state management ----------------------------------------
  
  
  // --- BEGIN -- channel information ------------------------------------------
  
  /// Returns a list of elements connected to the specified readout channel ID.
  /// @throws cet::exception (category: "Geometry") if non-existent channel
  std::vector<geo::WireID> ChannelToWire
    (raw::ChannelID_t channel) const override;
  
  /// Returns the total number of channels present (they are contiguous).
  virtual unsigned int Nchannels() const override;
  
  /// @brief Returns the number of channels in the specified readout plane.
  /// @return number of channels in the specified ROP, 0 if non-existent
  virtual unsigned int Nchannels(readout::ROPID const& ropid) const override;
  
  // --- END -- channel information --------------------------------------------
  
  
  // --- BEGIN -- channel topology ---------------------------------------------
  //@{
  /// @brief `WireCoordinate()` functions are not supported.
  virtual double WireCoordinate
    (double YPos, double ZPos, geo::PlaneID const& planeID) const override;
  virtual double WireCoordinate(double YPos, double ZPos,
                                unsigned int PlaneNo,
                                unsigned int TPCNo,
                                unsigned int cstat) const override
    { return WireCoordinate(YPos, ZPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
  //@}

  //@{
  /// @brief `NearestWireID()` functions are not supported.
  virtual geo::WireID NearestWireID
    (const TVector3& worldPos, geo::PlaneID const& planeID) const override;
  virtual geo::WireID NearestWireID(const TVector3& worldPos,
                                unsigned int    PlaneNo,
                                unsigned int    TPCNo,
                                unsigned int    cstat) const override
    { return NearestWireID(worldPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); }
  //@}
  
  // --- END -- channel topology -----------------------------------------------
  
  
  // --- BEGIN -- channel mapping ----------------------------------------------
  
  //@{
  /**
   * @brief Returns the channel ID a sensitive element is connected to.
   * @param wireID ID of the element
   * @return the ID of the channel
   * @see `PlaneWireToChannel(geo::WireID const&) const`
   *
   * Behaviour on an invalid or not present wires is undefined.
   *
   * @deprecated Use the version with `geo::WireID`
   */
  virtual raw::ChannelID_t PlaneWireToChannel
    (geo::WireID const& wireID) const override;
  virtual raw::ChannelID_t PlaneWireToChannel(
    unsigned int plane, unsigned int wire,
    unsigned int tpc, unsigned int cstat
    ) const override
    { return PlaneWireToChannel({ cstat, tpc, plane, wire }); }
  //@}

  /// @brief `PlaneIDs()` function is not supported.
  virtual std::set<PlaneID> const& PlaneIDs() const override;
  
  // --- END -- channel mapping ------------------------------------------------


  // --- BEGIN -- TPC set interface --------------------------------------------
  /// @name TPC set mapping
  /// @{
  /**
   * @brief Returns the total number of TPC sets in the specified cryostat.
   * @param cryoid cryostat ID
   * @return number of TPC sets in the cryostat, or 0 if no cryostat found
   *
   * In this mapping, TPCs have independent readout and there is one TPC in
   * each TPC set and one TPC set for each TPC.
   */
  virtual unsigned int NTPCsets
    (readout::CryostatID const& cryoid) const override;

  /// Returns the largest number of TPC sets any cryostat in the detector has.
  virtual unsigned int MaxTPCsets() const override;

  /// Returns whether we have the specified TPC set.
  /// @return whether the TPC set is valid and exists
  virtual bool HasTPCset(readout::TPCsetID const& tpcsetid) const override;

  /**
   * @brief Returns the ID of the TPC set the specified TPC belongs to.
   * @param tpcid ID of the TPC
   * @return the ID of the corresponding TPC set, or invalid ID when tpcid is
   *
   * In this mapping, TPC sets and TPCs are mapped one-to-one.
   * The returned value mirrors the TPC ID in the readout space.
   * If the TPC ID is not valid, an invalid TPC set ID is returned.
   * Note that this check is performed on the validity of the TPC ID, that
   * does not necessarily imply that the TPC specified by the ID actually
   * exists.
   */
  virtual readout::TPCsetID TPCtoTPCset
    (geo::TPCID const& tpcid) const override;

  /**
   * @brief Returns a list of ID of TPCs belonging to the specified TPC set.
   * @param tpcsetid ID of the TPC set to convert into TPC IDs
   * @return the list of TPCs, empty if TPC set is invalid
   *
   * In this mapping, TPC sets and TPCs are mapped one-to-one.
   * The returned list contains always one entry, unless the specified TPC
   * set ID is invalid, in which case the list is empty.
   * Note that the check is performed on the validity of the TPC set ID, that
   * does not necessarily imply that the TPC set specified by the ID actually
   * exists. Check the existence of the TPC set first (HasTPCset()).
   * Behaviour on valid, non-existent TPC set IDs is undefined.
   */
  virtual std::vector<geo::TPCID> TPCsetToTPCs
    (readout::TPCsetID const& tpcsetid) const override;

  /// Returns the ID of the first TPC belonging to the specified TPC set.
  virtual geo::TPCID FirstTPCinTPCset
    (readout::TPCsetID const& tpcsetid) const override;

  /// @}
  // --- END -- TPC set interface ----------------------------------------------


  // --- BEGIN -- Readout plane interface --------------------------------------
  /// @name Readout plane mapping
  /// @{

  /**
   * @brief Returns the total number of ROPs in the specified TPC set.
   * @param tpcsetid TPC set ID
   * @return number of readout planes in the TPC sets, or 0 if ID is invalid
   *
   * Note that this methods explicitly check the existence of the TPC set.
   *
   * In this mapping, planes have independent readout and there is one wire
   * plane in each readout plane and one readout plane for each wire plane.
   */
  virtual unsigned int NROPs
    (readout::TPCsetID const& tpcsetid) const override;

  /// Returns the largest number of ROPs a TPC set in the detector has.
  virtual unsigned int MaxROPs() const override;

  /// Returns whether we have the specified ROP.
  /// @return whether the readout plane is valid and exists
  virtual bool HasROP(readout::ROPID const& ropid) const override;

  /**
   * @brief Returns the ID of the ROP planeid belongs to, or invalid if none.
   * @param planeid ID of the plane
   * @return the ID of the corresponding ROP, or invalid ID when planeid is
   *
   * In this mapping, readout planes and wire planes are mapped one-to-one.
   * The returned value mirrors the plane ID in the readout space.
   * If the plane ID is not valid, an invalid readout plane ID is returned.
   * Note that this check is performed on the validity of the plane ID, that
   * does not necessarily imply that the plane specified by the ID actually
   * exists.
   */
  virtual readout::ROPID WirePlaneToROP
    (geo::PlaneID const& planeid) const override;

  /**
   * @brief Returns a list of ID of planes belonging to the specified ROP.
   * @param ropid ID of the readout plane to convert into wire planes
   * @return the list of wire plane IDs, empty if readout plane ID is invalid
   *
   * In this mapping, readout planes and wire planes are mapped one-to-one.
   * The returned list contains always one entry, unless the specified readout
   * plane ID is invalid, in which case the list is empty.
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists.
   */
  virtual std::vector<geo::PlaneID> ROPtoWirePlanes
    (readout::ROPID const& ropid) const override;

  /**
   * @brief Returns a list of ID of TPCs the specified ROP spans.
   * @param ropid ID of the readout plane
   * @return the list of TPC IDs, empty if readout plane ID is invalid
   *
   * In this mapping, readout planes and wire planes are mapped one-to-one.
   * The returned list contains always one entry, unless the specified readout
   * plane ID is invalid, in which case the list is empty.
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour on non-existing readout planes is undefined.
   */
  virtual std::vector<geo::TPCID> ROPtoTPCs
    (readout::ROPID const& ropid) const override;

  /// Returns the ID of the ROP the channel belongs to (invalid if none).
  virtual readout::ROPID ChannelToROP
    (raw::ChannelID_t channel) const override;

  /**
   * @brief Returns the ID of the first channel in the specified readout plane.
   * @param ropid ID of the readout plane
   * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
   *
   * Note that this check is performed on the validity of the readout plane
   * ID, that does not necessarily imply that the readout plane specified by
   * the ID actually exists. Check if the ROP exists with HasROP().
   * The behaviour for non-existing readout planes is undefined.
   */
  virtual raw::ChannelID_t FirstChannelInROP
    (readout::ROPID const& ropid) const override;

  /// Returns the ID of the first plane belonging to the specified ROP.
  virtual geo::PlaneID FirstWirePlaneInROP
    (readout::ROPID const& ropid) const override;

  /// @}
  // --- END -- Readout plane interface ----------------------------------------

  
  /// Sorter not supported.
  /// @throws geo::ChannelMapAlg::NoSorter
  [[noreturn]] virtual geo::GeoObjectSorter const& Sorter() const override;
  
  
    private:

  /// Information about the size of all geometry objects.
  struct GeometrySizes_t {
    
    GeometrySizes_t
      (unsigned int nCryostats, unsigned int maxTPCs, unsigned int maxPlanes)
      : nCryostats{ nCryostats }
      , nTPCs     { nCryostats, 0U }
      , nPlanes   { nCryostats, maxTPCs, 0U }
      , nWires    { nCryostats, maxTPCs, maxPlanes, 0U }
      , TPCmax { maxTPCs }
      , PlaneMax{ maxPlanes }
      {}
    
    /// Number of cryostats in the detector.
    unsigned int nCryostats = 0U;
    
    /// Number of TPCs in each cryostat.
    std::vector<unsigned int> nTPCs;
    
    /// Number of sensitive planes in in each TPC.
    geo::TPCDataContainer<unsigned int> nPlanes;
    
    /// Number of sensitive elements in each sensitive plane.
    geo::PlaneDataContainer<unsigned int> nWires;
    
#if 0
    // FIXME when latest geometry data containers land here
//     unsigned int maxTPCs() const { return nWires.dimSize<1U>(); }
    
    /// Returns the maximum number of sensitive planes per TPC.
    unsigned int maxPlanes() const { return nWires.dimSize<2U>(); }
#else
    unsigned int TPCmax = 0U;
    unsigned int PlaneMax = 0U;
      
    /// Returns the maximum number of TPCs per cryostat.
    unsigned int maxTPCs() const { return TPCmax; }
    
    /// Returns the maximum number of sensitive planes per TPC.
    unsigned int maxPlanes() const { return PlaneMax; }
#endif // 0?
    
    /// Returns whether any of the information is missing.
    bool empty() const
      {
        return (
          (nCryostats == 0U) || nTPCs.empty()
          || nPlanes.empty() || nWires.empty()
          );
      }
    
  }; // struct GeometrySizes_t
  
  
  /// A wrapper associating a datum to a channel. The datum is owned.
  template <typename T>
  struct ChannelOrderedDataStruct: std::pair<raw::ChannelID_t const, T> {
    
    using Data_t = T; ///< Type of data begin carried.
    
    using Base_t = std::pair<raw::ChannelID_t const, Data_t>; ///< Type of base class.
    
    /// Functor to compare two objects by channel.
    struct Comparer {
      
      /// Comparison between two objects.
      template <typename A, typename B>
      bool operator() (A const& a, B const& b) const { return cmp(a, b); }
      
      /// Returns the channel of... a channel.
      static raw::ChannelID_t channelOf(raw::ChannelID_t channel)
        { return channel; }
      
      /// Returns the channel of a `ChannelOrderedDataStruct` object.
      template <typename U>
      static raw::ChannelID_t channelOf(ChannelOrderedDataStruct<U> const& data)
        { return data.channel(); }
      
      /// Comparison between two objects (static).
      template <typename A, typename B>
      static bool cmp(A const& a, B const& b)
        { return channelOf(a) < channelOf(b); }
      
    }; // struct Comparer
    
    // Inherit constructors.
    using Base_t::Base_t;
    
    // @{
    /// Returns the channel used as index.
    raw::ChannelID_t channel() const { return Base_t::first; }
    raw::ChannelID_t& channel() { return Base_t::first; }
    // @}
    
    //@{
    /// Returns the data payload.
    Data_t& data() { return Base_t::second; }
    Data_t const& data() const { return Base_t::second; }
    //@}
    
    //@{
    /// Direct access to the data.
    Data_t* operator->() { return &data(); }
    Data_t const* operator->() const { return &data(); }
    //@}
    
    /// Simple comparison with another `ChannelOrderedDataStruct`.
    template <typename U>
    bool operator< (ChannelOrderedDataStruct<U> const& other) const
      { return Comparer::cmp(*this, other); }
    
  }; // struct ChannelOrderedDataStruct
  
  template <typename T>
  class ChannelDataContainer {
    
      public:
    
    using Data_t = T; /// Type of the contained data.
        
    /// Data associated with channel range start.
    using ChannelAndData_t = ChannelOrderedDataStruct<Data_t>;
    
    
    /// Returns a pointer to the datum including the specified channel.
    /// @return a pointer to the requested datum, or `nullptr` if not found
    ChannelAndData_t* find(raw::ChannelID_t channel);
    
    /// Returns a pointer to the datum including the specified channel.
    /// @return a pointer to the requested datum, or `nullptr` if not found
    ChannelAndData_t const* find(raw::ChannelID_t channel) const;
    
    
    /// Returns a pointer to the datum including the specified channel.
    /// @return a pointer to the requested datum, or `nullptr` if not found
    Data_t* findData(raw::ChannelID_t channel);
    
    /// Returns a pointer to the datum including the specified channel.
    /// @return a pointer to the requested datum, or `nullptr` if not found
    Data_t const* findData(raw::ChannelID_t channel) const;
    
    /// Returns the ID of the first invalid channel.
    raw::ChannelID_t endChannel() const { return fEndChannel; }
    
    /// Returns whether this container is empty (no channel range is covered).
    bool empty() const { return endChannel() == raw::ChannelID_t{ 0 }; }
    
    
    /// Sets the ID of first invalid channel (end of the covered channel range).
    void setEndChannel(raw::ChannelID_t channel) { fEndChannel = channel; }
    
    /// Adds a `datum` covering `nChannels` starting from `endChannel()`.
    void append(Data_t const& datum, unsigned int nChannels = 1U);
    
    /// Adds a `datum` for the channel range starting with `channel`.
    void append(Data_t&& datum, unsigned int nChannels = 1U);
    
    /// Removes all data.
    void clear();
    
      private:
    
    std::vector<ChannelAndData_t> fData; ///< Contained data.
    raw::ChannelID_t fEndChannel { 0 }; ///< First invalid channel.
    
  }; // struct ChannelDataContainer
  
  /// Data associated to the first channel in a channel range (sensitive plane).
  struct FirstChannelInfo_t {
    
    geo::PlaneID plane; ///< Wire corresponding to the first channel.
    
  }; // struct FirstChannelInfo_t
  
  /// Information per sensitive plane.
  struct PlaneInfo_t {
    
    /// First channel in the plane.
    raw::ChannelID_t firstChannel = raw::InvalidChannelID;
    
  }; // struct PlaneInfo_t
  
  
  /// Information of first channel of each sensitive plane.
  ChannelDataContainer<FirstChannelInfo_t> fFirstChannelInfo;
  
  // optional to be able to delay its initialization
  /// Number of subelements in each geometry component.
  std::optional<GeometrySizes_t> fGeometrySizes;
  
  // optional to allow delayed initialization
  /// Contains the information per plane.
  std::optional<geo::PlaneDataContainer<PlaneInfo_t>> fPlaneInfo;
  
  /// Returns the number of cryostats in the detector.
  unsigned int NCryostats() const
    { assert(fGeometrySizes); return fGeometrySizes->nCryostats; }
  
  /// Returns the number of TPCs in the specified cryostat.
  unsigned int NTPCs(geo::CryostatID const& cid) const
    { assert(fGeometrySizes); return fGeometrySizes->nTPCs[cid.Cryostat]; }
  
  /// Returns the number of sensitive planes in the specified TPC.
  unsigned int NPlanes(geo::TPCID const& tid) const
    { assert(fGeometrySizes); return fGeometrySizes->nPlanes[tid]; }
  
  /// Returns the number of sensitive elements in the specified sensitive plane.
  unsigned int NWires(geo::PlaneID const& pid) const
    { assert(fGeometrySizes); return fGeometrySizes->nWires[pid]; }
  
  /* // belongs to the interface
  /// Returns the number of TPC sets in the specified cryostat.
  unsigned int NTPCsets(geo::CryostatID const& cid) const
    { return NTPCs(cid); }
  
  /// Returns the number of readout planes in the specified TPC set.
  unsigned int NROPs(readout::TPCsetID const& sid) const
    { return NPlanes(convertTPCsetToTPC(sid)); }
  */
  
  /// Returns the number of channels in the specified readout plane.
  unsigned int NChannels(readout::ROPID const& rid) const
    { return NWires(convertROPtoPlane(rid)); }
  
  /// Returns the number of TPC channels in the detector.
  unsigned int NChannels() const
    {
      assert(!fFirstChannelInfo.empty());
      return fFirstChannelInfo.endChannel() - raw::ChannelID_t{ 0U };
    }
  
  /// Returns the maximum number of TPCs per cryostat.
  unsigned int MaxTPCs() const
    { assert(fGeometrySizes); return fGeometrySizes->maxTPCs(); }
  
  /// Returns the maximum number of sensitive planes per TPC.
  unsigned int MaxPlanes() const
    { assert(fGeometrySizes); return fGeometrySizes->maxPlanes(); }
  
  /// Returns whether the specified cryostat is present.
  bool hasCryostat(geo::CryostatID const& cid) const
    { return cid.Cryostat < NCryostats(); }
  
  /// Returns whether the specified TPC is present.
  bool hasTPC(geo::TPCID const& tid) const
    { return hasCryostat(tid) && tid.TPC < NTPCs(tid); }
  
  /// Returns whether the specified sensitive plane is present.
  bool hasPlane(geo::PlaneID const& pid) const
    { return hasTPC(pid) && pid.Plane < NPlanes(pid); }
  
  /// Returns whether the specified sensitive element is present.
  bool hasWire(geo::WireID const& wid) const
    { return hasPlane(wid) && wid.Wire < NWires(wid); }
  
  /// Returns whether the specified TPC set is present.
  bool hasTPCset(readout::TPCsetID const& sid) const
    { return hasTPC(convertTPCsetToTPC(sid)); }
  
  /// Returns whether the specified readout plane is present.
  bool hasROP(readout::ROPID const& rid) const
    { return hasPlane(convertROPtoPlane(rid)); }
  
  /// Returns whether the specified readout plane is present.
  bool hasChannel(raw::ChannelID_t const channel) const
    { return channel < NChannels(); }
  
  /// Returns the information for the specified plane.
  PlaneInfo_t const& planeInfo(geo::PlaneID const& pid) const
    { assert(fPlaneInfo); return fPlaneInfo.value()[pid]; }
  
  
  /**
    * @brief Return the signal type of the specified channel.
    * @param channel ID of the channel
    * @return signal type of the channel, or geo::kMysteryType if not known
    *
    * A valid sensitive element has signal type `geo::kCollection` (these are
    * supposed to be pixels).
    * On any type of error (e.g., invalid or unknown channel ID),
    * `geo::kMysteryType` is returned.
    */
  virtual geo::SigType_t SignalTypeForChannelImpl
    (raw::ChannelID_t const channel) const override;
  
  
  /// Extracts the first channel information from the geometry.
  void fillFirstChannelInfo
    (geo::GeometryData_t::CryostatList_t const& cryostats);
  
  /// Extracts the number of subelements of each detector element.
  void fillGeometrySizes(geo::GeometryData_t::CryostatList_t const& cryostats);
  
  /// Fills the map from plane to first channel.
  void makePlaneChannelMap();
  
  
  // --- BEGIN --- Conversion between IDs for trivial mapping ------------------
  /// @name Conversion between IDs for trivial readout/geometry element mapping
  /// @{
  
  /// Convert a TPC ID into a TPC set ID in a one-to-one mapping.
  static readout::TPCsetID convertTPCtoTPCset(geo::TPCID const& tpcid)
    {
      return {
        tpcid.asCryostatID(),
        static_cast<readout::TPCsetID::TPCsetID_t>(tpcid.TPC)
        };
    } // TPCtoTPCset()
  
  /// Convert a sensitive plane ID into a readout plane ID
  /// in a one-to-one mapping.
  static readout::ROPID convertPlaneToROP(geo::PlaneID const& planeid)
    {
      return {
        convertTPCtoTPCset(planeid),
        static_cast<readout::ROPID::ROPID_t>(planeid.Plane)
        };
    }
  
  /// Convert a TPC set ID into a TPC ID in a one-to-one mapping.
  static geo::TPCID convertTPCsetToTPC(readout::TPCsetID const& tpcsetid)
    {
      return {
        tpcsetid.asCryostatID(),
        static_cast<geo::TPCID::TPCID_t>(tpcsetid.TPCset)
        };
    } // TPCsetToTPC()
  
  /// Convert a readout plane ID into a sensitive plane ID
  /// in a one-to-one mapping.
  static geo::PlaneID convertROPtoPlane(readout::ROPID const& ropid)
    {
      return {
        convertTPCsetToTPC(ropid),
        static_cast<geo::PlaneID::PlaneID_t>(ropid.ROP)
        };
    }
  
  /// @}
  // --- END --- Conversion between IDs for trivial mapping --------------------
  
  
}; // class geo::ChannelMapPixelAlg


// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
template <typename T>
auto geo::ChannelMapPixelAlg::ChannelDataContainer<T>::find
  (raw::ChannelID_t channel) -> ChannelAndData_t*
{
  assert(!empty() || fData.empty());
  if (channel >= endChannel()) return nullptr;
  auto const iDatum = std::lower_bound(
    fData.begin(), fData.end(), channel, typename ChannelAndData_t::Comparer()
    );
  assert(iDatum != fData.end());
  return &*iDatum;
} // geo::ChannelMapPixelAlg::ChannelDataContainer<T>::find()


// -----------------------------------------------------------------------------
template <typename T>
auto geo::ChannelMapPixelAlg::ChannelDataContainer<T>::find
  (raw::ChannelID_t channel) const -> ChannelAndData_t const*
{
  assert(!empty() || fData.empty());
  if (channel >= endChannel()) return nullptr;
  auto const iDatum = std::lower_bound(
    fData.begin(), fData.end(), channel, typename ChannelAndData_t::Comparer{}
    );
  assert(iDatum != fData.end());
  return &*iDatum;
} // geo::ChannelMapPixelAlg::ChannelDataContainer<T>::find() const


// -----------------------------------------------------------------------------
template <typename T>
auto geo::ChannelMapPixelAlg::ChannelDataContainer<T>::findData
  (raw::ChannelID_t channel) -> Data_t*
{
  ChannelAndData_t* channelAndData = find(channel);
  return channelAndData? &(channelAndData->data()): nullptr;
} // geo::ChannelMapPixelAlg::ChannelDataContainer<T>::findData()


// -----------------------------------------------------------------------------
template <typename T>
auto geo::ChannelMapPixelAlg::ChannelDataContainer<T>::findData
  (raw::ChannelID_t channel) const -> Data_t const*
{
  ChannelAndData_t const* channelAndData = find(channel);
  return channelAndData? &(channelAndData.data()): nullptr;
} // geo::ChannelMapPixelAlg::ChannelDataContainer<T>::findData() const


// -----------------------------------------------------------------------------
template <typename T>
void geo::ChannelMapPixelAlg::ChannelDataContainer<T>::append
  (Data_t const& datum, unsigned int nChannels /* = 1U */)
{
  fData.emplace_back(endChannel(), datum);
  setEndChannel(endChannel() + nChannels);
} // geo::ChannelMapPixelAlg::ChannelDataContainer<T>::append()


// -----------------------------------------------------------------------------
template <typename T>
void geo::ChannelMapPixelAlg::ChannelDataContainer<T>::append
  (Data_t&& datum, unsigned int nChannels /* = 1U */)
{
  fData.emplace_back(endChannel(), datum);
  setEndChannel(endChannel() + nChannels);
} // geo::ChannelMapPixelAlg::ChannelDataContainer<T>::append()


// -----------------------------------------------------------------------------
template <typename T>
void geo::ChannelMapPixelAlg::ChannelDataContainer<T>::clear() {
  setEndChannel(raw::ChannelID_t{ 0 });
  fData.clear();
} // geo::ChannelMapPixelAlg::ChannelDataContainer<T>::clear()


// -----------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_PIXELPLANE_CHANNELMAPPIXELALG_H
