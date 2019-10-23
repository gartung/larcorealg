/**
 * @file   DriftPartitions.h
 * @brief  Data structures and algorithms to partition a cryostat volume.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 13, 2017
 * @see    DriftPartitions.cxx
 * @ingroup Geometry
 *
 */

#ifndef LARCOREALG_GEOMETRY_DRIFTPARTITIONS_H
#define LARCOREALG_GEOMETRY_DRIFTPARTITIONS_H

// LArSoft libraries
#include "larcorealg/Geometry/Partitions.h"
#include "larcorealg/Geometry/Decomposer.h"
#include "larcorealg/Geometry/SimpleGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "Math/GenVector/Cartesian2D.h"
#include "Math/GenVector/DisplacementVector2D.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <functional> // std::less()
#include <algorithm> // std::upper_bound()
#include <utility> // std::forward()
#include <memory> // std::unique_ptr()
#include <type_traits> // std::decay_t<>, std::declval()


namespace geo {

  class CryostatGeo;

  namespace part {

    //--------------------------------------------------------------------------
    /// Describes a `geo::TPCGeo` object for `Partition::describe()`.
    template <>
    struct PartitionDataDescriber<geo::TPCGeo> {
      template <typename Stream>
      PartitionDataDescriber(
        Stream&& out, geo::TPCGeo const* pTPC,
        std::string indent = "", std::string firstIndent = ""
        );
    }; // PartitionDataDescriber(TPCGeo)

    //--------------------------------------------------------------------------

  } // namespace part


  namespace details {

    //--------------------------------------------------------------------------
    /// Function translation of `std::less`.
    template <typename T>
    auto static_less(T a, T b)
      { return std::less<T>()(a, b); }

    //--------------------------------------------------------------------------
    /// Class managing comparisons between `T` objects via a `Key` key.
    template <
      typename T, typename Key,
      Key KeyExtractor(T const&),
      bool KeyComparer(Key, Key) = static_less<Key>
      >
    struct Comparer;

    //--------------------------------------------------------------------------

  } // namespace details


  // --- BEGIN -----------------------------------------------------------------
  /// @ingroup Geometry
  /// @{
  /**
  * @brief Set of drift volumes.
  *
  * A drift volume is a set of TPCs whose readout planes lie on the same
  * geometric plane.
  *
  */
  class DriftPartitions {

      public:
    /// Type of TPC collection for the partition of a single drift volume.
    using TPCPartition_t = geo::part::Partition<geo::TPCGeo const>;

    /// Type for description of drift range.
    using Range_t = lar::util::simple_geo::Range<double>;

    /// Data associated to a single drift volume.
    struct DriftVolume_t {
      /// A partition of the volume in width and depth.
      std::unique_ptr<TPCPartition_t> partition;
      /// Interval of drift direction covered by this drift volume.
      Range_t driftCoverage;

      /// Constructor: imports the specified partition and drift coverage range.
      DriftVolume_t
        (std::unique_ptr<TPCPartition_t>&& part, Range_t const& cover)
        : partition(std::move(part)), driftCoverage(cover) {}

      /// Returns whether this drift volume covers specified drift coordinate.
      bool coversDrift(double drift) const
        { return driftCoverage.contains(drift); }

      /// Returns the drift coordinate of the specified partition.
      static double Position(DriftVolume_t const& part)
        { return part.driftCoverage.lower; }

      /// Type of static object to compare `DriftVolume_t` objects.
      using Comparer_t
        = details::Comparer<DriftVolume_t, double, DriftVolume_t::Position>;

    }; // DriftPartitions::DriftVolume_t


    /// Type representing a position in 3D space.
    using Position_t = geo::Point_t;

    /// Type representing a direction in 3D space (norm is not constrained).
    using Direction_t = geo::Vector_t;

    /// Type representing a position in the 2D space.
    using Projection_t
      = ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<double>>;

    /// Type representing the drift direction (assumed to have norm 1).
    using DriftDir_t = Direction_t;

    /// Object used to compute projections on drift volume readout plane.
    using Decomposer_t = geo::Decomposer<Direction_t, Position_t, Projection_t>;


    /// All drift volumes, sorted by position.
    std::vector<DriftVolume_t> volumes;
    Decomposer_t decomposer; ///< Decomposition on drift, width and depth axes.

    /// Constructor: no partition, but sets the main "drift" direction.
    explicit DriftPartitions(Decomposer_t const& decomp): decomposer(decomp) {}

    /// @{
    /// @name Drift volume lookup.

    /// Returns drift coordinate (in the drift-volume-specific frame) of `pos`.
    double driftCoord(Position_t const& pos) const
      { return decomposer.PointNormalComponent(pos); }

    /// Returns which partition contains the specified position.
    /// @return volume containing the specified position (`nullptr` if none)
    DriftVolume_t const* driftVolumeAt(Position_t const& pos) const
      { return driftVolumeAt(driftCoord(pos)); }

    /// Returns which volume contains the specified drift (`nullptr` if none).
    DriftVolume_t const* driftVolumeAt(double drift) const;

    /// Returns which TPC contains the specified position (`nullptr` if none).
    geo::TPCGeo const* TPCat(Position_t const& pos) const;

    /// @}

    /// Printout of the drift volume information.
    template <typename Stream>
    void print(Stream&& out) const;

    /// Adds the specified partition as a new drift volume.
    void addPartition(std::unique_ptr<TPCPartition_t>&& part);


      private:

    /// Returns an iterator to the drift volume starting after `pos`.
    std::vector<DriftVolume_t>::iterator volumeAfter(double pos);

    /// Returns an iterator to the drift volume starting after `pos`.
    std::vector<DriftVolume_t>::const_iterator volumeAfter(double pos) const;

    /// Computes the coverage of the specified partition in the drift direction.
    Range_t computeCoverage(TPCPartition_t const& TPCpart) const;

  }; // class DriftPartitions


  //----------------------------------------------------------------------------
  /**
   * @brief Creates a `DriftPartitions` object from the TPCs in a cryostat.
   * @param cryo the cryostat with the TPCs to be partitioned.
   *
   * The algorithm groups all TPCs with the same drift direction and readout
   * planes with similar drift coordinates.
   * A "global drift direction" is chosen. The drift directions of all TPCs in
   * the cryostat are required to be parallel (at most with different verse).
   * Projections of the TPC wire planes on this direction determine TPC
   * grouping. TPCs with planes within a small range (typically, 10 plane
   * pitches) of this projected coordinate are grouped together.
   * TPCs with opposite drift directions are still kept in different groups.
   *
   * The information of each drift volume is kept in a `DriftVolume_t` class
   * which contains a hierarchical structure of type `geo::part::Partition`
   * (with data `geo::TPCGeo const` coming directly from the geometry). This
   * structure describes the topology of the TPCs within the drift volume.
   *
   */
  DriftPartitions buildDriftVolumes(geo::CryostatGeo const& cryo);


  /// @}
  // --- END -------------------------------------------------------------------

} // namespace geo


//------------------------------------------------------------------------------
//---  inline and template implementation
//---
//------------------------------------------------------------------------------
namespace geo {
  namespace details {

    //--------------------------------------------------------------------------
    // TODO use SorterByKey instead?
    template <
      typename T, typename Key,
      Key KeyExtractor(T const&),
      bool KeyComparer(Key, Key) /* = static_less<Key> */
      >
    struct Comparer {
      using Key_t = Key; ///< Type of comparison key.
      using Object_t = T; ///< Type of object to be compared.

      bool operator() (Key_t a, Key_t b) const { return key_comp(a, b); }
      bool operator() (Object_t const& a, Key_t b) const
        { return this->operator() (key(a), b); }
      bool operator() (Key_t a, Object_t const& b) const
        { return this->operator() (a, key(b)); }
      bool operator() (Object_t const& a, Object_t const& b) const
        { return this->operator() (key(a), b); }

        private:
      static auto key(T const& v) { return KeyExtractor(v); }
      static auto key_comp(Key_t a, Key_t b) { return KeyComparer(a, b); }

    }; // struct Comparer<>

    //--------------------------------------------------------------------------

  } // namespace details
} // namespace geo


//------------------------------------------------------------------------------
template <typename Stream>
geo::part::PartitionDataDescriber<geo::TPCGeo>::PartitionDataDescriber(
  Stream&& out, geo::TPCGeo const* pTPC,
  std::string indent /* = "" */, std::string firstIndent /* = "" */
  )
{
  out << firstIndent;
  if (!pTPC) {
    out << "no TPC";
    return;
  }
  // verbosity: 2 => ID, size, drift information;
  // 5: also number of wires and active volume
  pTPC->PrintTPCInfo(std::forward<Stream>(out), indent, 2U);
} // geo::part::PartitionDataDescriber<TPCGeo>::PartitionDataDescriber()


//------------------------------------------------------------------------------
//--- geo::DriftPartitions
//---
inline auto geo::DriftPartitions::volumeAfter(double pos)
  -> std::vector<DriftVolume_t>::iterator
{
  return std::upper_bound
    (volumes.begin(), volumes.end(), pos, DriftVolume_t::Comparer_t());
} // geo::DriftPartitions::volumeAfter()

inline auto geo::DriftPartitions::volumeAfter(double pos) const
  -> std::vector<DriftVolume_t>::const_iterator
{
  return std::upper_bound
    (volumes.begin(), volumes.end(), pos, DriftVolume_t::Comparer_t());
} // geo::DriftPartitions::volumeAfter()


//------------------------------------------------------------------------------
template <typename Stream>
void geo::DriftPartitions::print(Stream&& out) const {

  out << volumes.size() << " drift volume partitions:";
  for (auto const& driftVol: volumes) {
    out << "\n[" << driftVol.driftCoverage.lower
      << " -- " << driftVol.driftCoverage.upper << "]: "
      << driftVol.partition->describe("  ", "");
  } // for
  out << "\n";
} // geo::DriftPartitions::print()


//------------------------------------------------------------------------------

#endif // LARCOREALG_GEOMETRY_DRIFTPARTITIONS_H
