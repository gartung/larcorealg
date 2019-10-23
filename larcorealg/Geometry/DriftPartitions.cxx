/**
 * @file   DriftPartitions.cxx
 * @brief  Classes describing partition of cryostat volume.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 13, 2017
 * @see    DriftPartitions.h
 *
 */

#include "larcorealg/Geometry/DriftPartitions.h"


// LArSoft libraries
#include "larcorealg/Geometry/PlaneGeo.h"                  // for PlaneGeo::...
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/Decomposer.h"
#include "larcorealg/CoreUtils/DumpUtils.h" // lar::dump::vector3D(), ...
#include "larcorealg/CoreUtils/RealComparisons.h"

// utility libraries
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <vector>
#include <iterator> // std::back_inserter(), std::distance(), ...
#include <algorithm> // std::transform(), std::sort(), ...
#include <utility> // std::pair<>, std::move()
#include <memory> // std::addressof(), std::unique_ptr<>, ...
#include <type_traits> // std::declval(), ...
#include <cassert>
#include <cmath> // std::abs()
#include <cstdlib> // std::size_t

//------------------------------------------------------------------------------
namespace {

  //----------------------------------------------------------------------------
  /// Appends a copy of the content of `src` to the end of `dest` vector.
  template <typename T, typename CONT>
  void append(std::vector<T>& dest, CONT const& src) {
    using std::cbegin;
    using std::cend;
    dest.insert(dest.end(), cbegin(src), cend(src));
  } // append()

  /// Copies or moves the content of `src` vector to the end of `dest` vector.
  template <typename T>
  void append(std::vector<T>& dest, std::vector<T>&& src) {
    if (dest.empty()) dest = std::move(src);
    else append(dest, src);
  } // append(vector<T>&&)


  //----------------------------------------------------------------------------
  /// Returns an iterable through all keys in `m` (map or pair collection).
  template <typename Cont>
  auto keys(Cont const& m) {
    // TODO use range library instead (this temporary implementation copies keys)
    std::vector<typename Cont::value_type::first_type> result;
    result.reserve(m.size());
    std::transform(
      m.cbegin(), m.cend(),
      std::back_inserter(result), [](auto const& p){ return p.first; }
      );
    return result;
  } // keys()


  //----------------------------------------------------------------------------

} // local namespace


//------------------------------------------------------------------------------
//---  geo::DriftPartitions
//------------------------------------------------------------------------------
geo::DriftPartitions::DriftVolume_t const* geo::DriftPartitions::driftVolumeAt
  (double drift) const
{
  auto iVol = volumeAfter(drift); // points to partition after the good one
  if (iVol == volumes.cbegin()) return nullptr;
  if (!(--iVol)->coversDrift(drift)) return nullptr; // maybe not that good?
  return &*iVol;
} // geo::DriftPartitions::driftVolumeAt()


//------------------------------------------------------------------------------
geo::TPCGeo const* geo::DriftPartitions::TPCat(Position_t const& pos) const {
  auto comp = decomposer.DecomposePoint(pos);
  auto volume = driftVolumeAt(comp.distance);
  return (volume && volume->partition)
    ? volume->partition->atPoint(comp.projection.X(), comp.projection.Y())
    : nullptr;
} // geo::DriftPartitions::TPCat()


//------------------------------------------------------------------------------
void geo::DriftPartitions::addPartition(std::unique_ptr<TPCPartition_t>&& part)
{
  auto const range = computeCoverage(*part);
  auto iVol = volumeAfter(range.lower);
  volumes.emplace(iVol, std::move(part), range);
} // geo::DriftPartitions::addPartition()


//------------------------------------------------------------------------------
auto geo::DriftPartitions::computeCoverage
  (TPCPartition_t const& TPCpart) const -> Range_t
{
  /*
   * Computes the range of drift covered by the TPCs in the specified partition.
   * The range is currently defined to include any drift distance covered by
   * any single TPC.
   * The drift distance is computed in "absolute" coordinates, meaning that the
   * origin of the drift is at the origin of the global coordinate system.
   * The drift direction is the one stored in this object.
   */

  struct CoverageExtractor {
    Range_t coverage;

    CoverageExtractor(DriftPartitions::Decomposer_t const& decomp)
      : decomp(decomp) {}

    void operator() (TPCPartition_t const& TPCpart)
      {
        geo::TPCGeo const* TPC = TPCpart.data();
        if (!TPC) return;
        includePoint(TPC->GetCathodeCenter<Position_t>());
        includePoint(TPC->LastPlane().GetCenter<Position_t>());
      }

      private:
    DriftPartitions::Decomposer_t const& decomp;

    double driftCoord(Position_t const& pos) const
      { return decomp.PointNormalComponent(pos); }
    void includePoint(Position_t const& pos)
      { coverage.extendToInclude(driftCoord(pos)); }
  }; // struct CoverageExtractor

  CoverageExtractor extractor(decomposer);
  TPCpart.walk(extractor);
  return extractor.coverage;

} // DriftPartitions::computeCoverage()


//******************************************************************************
//***  Algorithms for building the drift volume partition.
//******************************************************************************
//------------------------------------------------------------------------------
// Keeps a TPC and its position as a single coordinate on a axis.
using TPCandPos_t = std::pair<double, geo::TPCGeo const*>;

// Data structure collecting all TPCs in a drift volume.
struct TPCgroup_t {
  double pos; // Common coordinate of the group.
  std::vector<geo::TPCGeo const*> TPCs; // All TPCs in the group.

  // Explicit constructor.
  TPCgroup_t(double pos, std::vector<geo::TPCGeo const*>&& TPCs)
    : pos(pos), TPCs(std::move(TPCs)) {}

  // Returns the position of a TPC group.
  static double Position(TPCgroup_t const& tpcg) { return tpcg.pos; }

  // Comparison object useful to sort and find TPC groups.
  using Comparer_t
    = geo::details::Comparer<TPCgroup_t, double, TPCgroup_t::Position>;

}; // struct TPCgroup_t


//------------------------------------------------------------------------------
std::vector
  <std::pair<geo::DriftPartitions::DriftDir_t, std::vector<geo::TPCGeo const*>>>
groupTPCsByDriftDir(geo::CryostatGeo const& cryo)
{
  /*
   * TPCs are grouped by their drift direction, allowing for a small rounding
   * error.
   * The result is a collection with one element for each group of TPCs sharing
   * the same drift direction. Each element is a pair with that drift direction
   * first, and a collection of all TPCs with that drift direction.
   * Elements are in no particular order.
   */
  std::vector<std::pair
    <geo::DriftPartitions::DriftDir_t, std::vector<geo::TPCGeo const*>>
    >
    result;

  lar::util::RealComparisons<double> coordIs(1e-4);
  auto vectorIs = lar::util::makeVector3DComparison(coordIs);

  auto const nTPCs = cryo.NTPC();
  for (unsigned int iTPC = 0; iTPC < nTPCs; ++iTPC) {
    geo::TPCGeo const& TPC = cryo.TPC(iTPC);

    decltype(auto) driftDir = TPC.DriftDir();

    std::size_t iGroup = 0;
    for (; iGroup < result.size(); ++iGroup) {
      if (vectorIs.nonEqual(driftDir, result[iGroup].first)) continue;
      result[iGroup].second.push_back(&TPC);
      break;
    } // for

    // if we did not find a group yet, make a new one
    if (iGroup == result.size()) {
      result.emplace_back(
        geo::vect::rounded01(driftDir, coordIs.threshold),
        std::vector<geo::TPCGeo const*>{ &TPC }
        );
    } // if

  } // for

  return result;
} // groupTPCsByDriftDir()


//------------------------------------------------------------------------------
std::vector<TPCandPos_t> sortTPCsByDriftCoord(
  std::vector<geo::TPCGeo const*> const& TPCs,
  geo::DriftPartitions::Decomposer_t const& decomp
) {
  /*
   * TPCs in the argument are sorted by their drift coordinate.
   * The drift coordinate is defined as the coordinate specified by the normal
   * direction of the decomposer in argument.
   * The result is a collection of data structures containing each a TPC and
   * the drift coordinate that was used to sort it.
   *
   * Sorting happens by the drift coordinate of the first wire plane;
   * the absolute value of the drift coordinate value is not relevant nor it is
   * well defined.
   * The result preserves that coordinate for further processing (grouping).
   */
  auto const driftCoord = [&decomp](geo::TPCGeo const& TPC)
    { return decomp.PointNormalComponent(geo::vect::convertTo<geo::DriftPartitions::Position_t>(TPC.FirstPlane().GetCenter())); };

  std::vector<TPCandPos_t> result;
  result.reserve(TPCs.size());
  std::transform(TPCs.cbegin(), TPCs.cend(), std::back_inserter(result),
    [&driftCoord](geo::TPCGeo const* pTPC)
      { return TPCandPos_t(driftCoord(*pTPC), pTPC); }
    );
  // std::pair sorts by first key first, and second key on par
  // (on par, which may happen often, we don't have means to decide here...)
  std::sort(result.begin(), result.end());
  return result;
} // sortTPCsByDriftCoord()


//------------------------------------------------------------------------------
std::vector<TPCgroup_t> groupByDriftCoordinate
  (std::vector<TPCandPos_t> const& TPCs)
{
  /*
   * Produces a list of TPC groups. Within each group, the TPCs have similar
   * drift coordinate.
   * Similar is defined arbitrarily as ten times the plane pitch (of the first
   * planes of the TPC).
   */
  if (TPCs.empty()) return {};

  geo::TPCGeo const& firstTPC = *(TPCs.front().second);
  // arbitrary 5 cm if the first TPC has only one plane (pixel readout?);
  // protect against the case where planes have the same position
  // (e.g. dual phase)
  double const groupThickness = 10.0
    * std::min(((firstTPC.Nplanes() > 1)? firstTPC.Plane0Pitch(1): 0.5), 0.1);

  auto iFirstTPC = TPCs.cbegin(), tend = TPCs.cend();

  std::vector<TPCgroup_t> result;
  while (iFirstTPC != tend) {
    double const posEnd = iFirstTPC->first + groupThickness; // not beyond here
    double sumPos = 0.0;
    std::vector<geo::TPCGeo const*> TPCs;
    auto iEndGroup = iFirstTPC;
    do {
      TPCs.push_back(iEndGroup->second);
      sumPos += iEndGroup->first;
      ++iEndGroup;
    } while ((iEndGroup != tend) && (iEndGroup->first < posEnd));

    double const averagePos = sumPos / TPCs.size();
    result.emplace_back(averagePos, std::move(TPCs));

    iFirstTPC = iEndGroup;
  } // while (outer)

  return result;
} // groupByDriftCoordinate()


//------------------------------------------------------------------------------
unsigned int checkTPCcoords(std::vector<geo::TPCGeo const*> const& TPCs) {
  /*
   * Verify coordinate system consistency between TPCs:
   *   * need to have the same drift direction
   *   * need to have the same drift coordinate
   *
   * On error, it prints information on the error stream ("GeometryPartitions").
   * It returns the number of errors found.
   */

  auto iTPC = TPCs.cbegin(), tend = TPCs.cend();
  if (iTPC == tend) {
    mf::LogProblem("GeometryPartitions")
      << "checkTPCcoords() got an empty partition.";
    return 0;
  }

  geo::TPCGeo const& refTPC = **iTPC;
  decltype(auto) refDriftDir = refTPC.DriftDir();

  auto driftCoord = [&refDriftDir](geo::TPCGeo const& TPC)
    { return geo::vect::dot(TPC.FirstPlane().GetCenter(), refDriftDir); };

  auto const refDriftPos = driftCoord(refTPC);

  lar::util::RealComparisons<double> coordIs(1e-5);
  auto vectorIs = lar::util::makeVector3DComparison(coordIs);

  unsigned int nErrors = 0U;
  while (++iTPC != tend) {
    geo::TPCGeo const& TPC = **iTPC;

    if (vectorIs.nonEqual(TPC.DriftDir(), refDriftDir)) {
      mf::LogProblem("GeometryPartitions")
        << "Incompatible drift directions between " << TPC.ID()
        << " " << lar::dump::vector3D(TPC.DriftDir()) << " and " << refTPC.ID()
        << " " << lar::dump::vector3D(refTPC.DriftDir());
      ++nErrors;
    }
    auto const driftPos = driftCoord(TPC);
    if (coordIs.nonEqual(driftPos, refDriftPos)) {
      mf::LogProblem("GeometryPartitions")
        << "Incompatible drift coordinate between " << TPC.ID()
        << " (" << driftPos << "( and " << refTPC.ID() << " ("
        << refDriftPos << ")";
      ++nErrors;
    }
  } // while
  return nErrors;
} // checkTPCcoords()


//------------------------------------------------------------------------------
template <typename Range>
geo::DriftPartitions::DriftDir_t detectGlobalDriftDir(Range&& directions) {
  /*
   * Returns the first of the specified directions (possibly flipped);
   * throws an exception if any of them is not parallel to that one.
   */
  using std::cbegin;
  using std::cend;
  auto iDir = cbegin(directions);
  auto dend = cend(directions);
  if (!(iDir != dend)) {
    throw cet::exception("buildDriftVolumes")
      << "detectGlobalDriftDir(): no TPCs provided!\n";
  }

  lar::util::RealComparisons<double> comp(1e-5);
  auto compatibleDir = [comp](auto const& a, auto const& b)
    { return comp.equal(std::abs(geo::vect::dot(a, b)), +1.0); };

  auto const dir = *(iDir++);
  for (; iDir != dend; ++iDir) {
    if (compatibleDir(dir, *iDir)) continue;
    throw cet::exception("buildDriftVolumes")
      << "Found drift directions not compatible: " << lar::dump::vector3D(dir)
      << " and " << lar::dump::vector3D(*iDir) << "\n";
  } // for

  // mildly prefer positive directions
  return ((dir.X() <= 0.0) && (dir.Y() <= 0.0) && (dir.Z() <= 0.0))? -dir: dir;
} // detectGlobalDriftDir()


//------------------------------------------------------------------------------
// A TPC and the area it covers in the partition.
struct TPCwithArea_t: public geo::part::AreaOwner {
  using Area_t = geo::part::AreaOwner::Area_t;

  geo::TPCGeo const* TPC = nullptr;

  TPCwithArea_t(Area_t area, geo::TPCGeo const* TPC)
    : geo::part::AreaOwner(area), TPC(TPC) {}

}; // TPCwithArea_t


//------------------------------------------------------------------------------
geo::part::AreaOwner::Area_t TPCarea
  (geo::TPCGeo const& TPC, geo::DriftPartitions::Decomposer_t const& decomposer)
{
  /*
   * Returns the "area" of the TPC.
   * The area is delimited by TPC bounding box
   */

  geo::DriftPartitions::Position_t const lower
    { TPC.MinX(), TPC.MinY(), TPC.MinZ() };
  geo::DriftPartitions::Position_t const upper
    { TPC.MaxX(), TPC.MaxY(), TPC.MaxZ() };

  auto const lowerProj = decomposer.ProjectPointOnPlane(lower);
  auto const upperProj = decomposer.ProjectPointOnPlane(upper);

  // we ask to sort the ranges, since the reference base may be flipped
  return {
    { lowerProj.X(), upperProj.X(), true },
    { lowerProj.Y(), upperProj.Y(), true }
    };

} // TPCarea()


std::vector<TPCwithArea_t> addAreaToTPCs(
  std::vector<geo::TPCGeo const*> const& TPCs,
  geo::DriftPartitions::Decomposer_t const& decomposer
) {
  /*
   * Transforms a collection of TPCs into a collection of TPCs with area.
   */
  std::vector<TPCwithArea_t> result;
  result.reserve(TPCs.size());

  for (auto const& TPC: TPCs)
    result.emplace_back(TPCarea(*TPC, decomposer), TPC);

  return result;
} // addAreaToTPCs()


//------------------------------------------------------------------------------
template <typename BeginIter, typename EndIter>
geo::part::AreaOwner::Area_t computeTotalArea
  (BeginIter TPCbegin, EndIter TPCend)
{
  /*
   * Computes an area covering the areas from all TPCs delimited by the
   * iterators in argument.
   * Each iterator points to a AreaOwner pointer.
   */
  auto iTPC = TPCbegin;
  geo::part::AreaOwner::Area_t totalArea((*iTPC)->area());
  while (++iTPC != TPCend) totalArea.extendToInclude((*iTPC)->area());
  return totalArea;
} // computeTotalArea()


//------------------------------------------------------------------------------
// Class applying a comparison function to keys from the arguments.
template <typename Key, typename ExtractKey, typename Comparer = std::less<Key>>
struct SorterByKey {
  // differently from geo::details::Comparer, this class is not focused on a
  // specific type and its ExtractKey is not required to be a function operating
  // on that object.
  using Key_t = Key;

  static bool sortKey(Key_t a, Key_t b) { return Comparer()(a, b); }
  static Key_t key(Key_t k) { return k; }
  template <typename Data>
  static Key_t key(Data const& obj) { return ExtractKey()(obj); }

  template <typename A, typename B>
  bool operator() (A const& a, B const& b) const
    { return sortKey(key(a), key(b)); }

}; // struct SorterByKey


//------------------------------------------------------------------------------
// Class sorting any datum by a TPC area range lower boundary.
template<geo::part::AreaOwner::AreaRangeMember_t Range>
using SortTPCareaByAreaRangeLower
  = SorterByKey<double, geo::part::details::RangeLowerBoundExtractor<Range>>;

using SortTPCwithAreaByWidth
  = SortTPCareaByAreaRangeLower<&geo::part::AreaOwner::Area_t::width>;
using SortTPCwithAreaByDepth
  = SortTPCareaByAreaRangeLower<&geo::part::AreaOwner::Area_t::depth>;


//------------------------------------------------------------------------------
template<
  geo::part::AreaOwner::AreaRangeMember_t sortingRange,
  typename BeginIter, typename EndIter
  >
std::vector<std::vector<TPCwithArea_t const*>::const_iterator>
groupTPCsByRangeCoord
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  /*
   * Groups each TPC with all the following ones overlapping in the selected
   * sorting range.
   * The range of TPCs must be already sorted by lower sorting range coordinate.
   * The result is a list of iterators like the ones in input (in fact, the
   * first is always beginTPCwithArea).
   * The iterators are expected to be valid also after this function has
   * returned (lest the result be unusable).
   */

  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> groupStart;

  // tolerate 1mm overlaps; this is way forgiving, but apparently there are
  // geometries around (DUNE 35t) with overlaps of that size.
  lar::util::RealComparisons<double> coordIs(0.1);

  auto gbegin = beginTPCwithArea;
  while (gbegin != endTPCwithArea) {

    groupStart.push_back(gbegin);

    //
    // detect the end of this group
    //
    auto range = (*gbegin)->area().*sortingRange;
    auto gend = gbegin;
    while (++gend != endTPCwithArea) {
      //
      // check if the sorting range of this TPC (gend) overlaps the accumulated
      // one; since TPCs are sorted by lower end of the range, gend has that one
      // larger than the accumulated one, and overlap happens only if that lower
      // bound is smaller than the upper bound of the accumulated range;
      // we need to avoid rounding errors: close borders are decided as
      // non-overlapping
      //
      auto const& TPCrange = (*gend)->area().*sortingRange;
      if (coordIs.nonSmaller(TPCrange.lower, range.upper))
        break;
      range.extendToInclude(TPCrange);
    } // while (inner)

    // prepare for the next group
    gbegin = gend;
  } // while (outer)

  return groupStart;

} // groupTPCsByRangeCoord<>()


//------------------------------------------------------------------------------
// Returns the TPCs sorted by range coordinate, and the group limits.
template<
  geo::part::AreaOwner::AreaRangeMember_t sortingRange,
  typename BeginIter, typename EndIter
  >
std::pair<
  std::vector<TPCwithArea_t const*>,
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator>
  >
sortAndGroupTPCsByRangeCoord
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  //
  // sort by coordinate; work with pointers for convenience
  //
  std::vector<TPCwithArea_t const*> TPCs(beginTPCwithArea, endTPCwithArea);
  if (TPCs.size() <= 1) return {}; // with only one TPC, refuse to operate
  std::sort
    (TPCs.begin(), TPCs.end(), SortTPCareaByAreaRangeLower<sortingRange>());

  //
  // group
  //
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> TPCgroups
    = groupTPCsByRangeCoord<sortingRange>(TPCs.cbegin(), TPCs.cend());
  assert(!TPCgroups.empty());

  return { std::move(TPCs), std::move(TPCgroups) };

} // sortAndGroupTPCsByRangeCoord()


//------------------------------------------------------------------------------
template <
  typename BeginIter, typename EndIter, typename TPCendIter,
  typename SubpartMaker
  >
geo::part::Partition<geo::TPCGeo const>::Subpartitions_t createSubpartitions(
  BeginIter itTPCbegin, EndIter itTPCend, TPCendIter TPCend,
  SubpartMaker subpartMaker
) {
  /*
   * Internal helper to create a sequence of partitions (geo::part::Partition)
   * from groups of TPCs. Each TPC is specified as a pointer to TPCwithArea_t
   * object.
   *
   * The groups are specified in a way that is more or less convenient after
   * calling groupTPCsByRangeCoord(): the uber-iterators itTPCbegin and itTPCend
   * delimit a collection of iterators pointing to the first TPC of a group
   * (the pointed TPCs are on a different collection). This can define all
   * groups, each group delimited by the TPC pointed by an iterator pointing at
   * the beginning of that group and the TPC pointed by the next iterator,
   * pointing at the beginning of the next group. The exception is the last
   * group, for which there is no "next iterator pointing to the beginning of
   * the next group". That information is provided separately by the third
   * iterator argument, that is an iterator of a different type than the other
   * two (which are "uber-iterators" whose values are iterators), and which
   * points to the last available TPC, that is also the end iterator for the
   * TPCs of the last group.
   *
   */
  geo::part::Partition<geo::TPCGeo const>::Subpartitions_t subparts;

  // TPCgroups has an iterator to the starting TPC of each group.
  // Iterators refer to the `TPCs` collection.
  // The end iterator of the group is the starting iterator of the next one;
  // the last group includes all the remaining TPCs and its end iterator is
  // the end iterator of `TPCs`.
  auto igbegin = itTPCbegin;
  while (igbegin != itTPCend) {
    auto const gbegin = *igbegin;
    auto const gend = (++igbegin == itTPCend)? TPCend: *igbegin;

    //
    // create a partition from the new group
    //
    if (std::distance(gbegin, gend) == 1) {
      subparts.emplace_back(makeTPCPartitionElement(**gbegin));
    }
    else {
      auto subpart = subpartMaker(gbegin, gend);
      if (!subpart) return {}; // failure!!

      subparts.emplace_back(std::move(subpart));
    }
  } // while
  return subparts;
} // createSubpartitions()


//------------------------------------------------------------------------------
auto makeTPCPartitionElement(TPCwithArea_t const& TPCinfo)
{
  return std::make_unique<geo::part::PartitionElement<geo::TPCGeo const>>
    (TPCinfo.area(), TPCinfo.TPC);
}


//------------------------------------------------------------------------------
// prototypes for makePartition-type functions: they call each other.
template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makeWidthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea);
template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makeDepthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea);
template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makeGridPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea);
template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makePartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea);


//------------------------------------------------------------------------------
template <
  typename TPCPartitionResultType,
  geo::part::AreaOwner::AreaRangeMember_t Range,
  typename BeginIter, typename EndIter,
  typename SubpartMaker
  >
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makeSortedPartition(
  BeginIter beginTPCwithArea, EndIter endTPCwithArea,
  SubpartMaker subpartMaker
) {
  /*
   * TPCs in input are arranged into a partition split on width direction.
   * In case of failure, a null pointer is returned.
   * Do not use this function for a single TPC.
   * The algorithm is as follows:
   *
   * 1. sort the TPCs by width coordinate
   * 2. for each one, group it with all the following ones overlapping in width
   * 3. for each group with:
   *     .1 more than one element:
   *         .1 let the group be partitioned on depth
   *         .2 if that fails, we fail
   *     .2 just one element: create a partition element
   * 4. the resulting partition is the list of one-TPC "groups" and depth-based
   *    subpartitions of many-TPC groups
   * 5. if the result is a single partition, it must be a depth partition;
   *    then, there is no room for a partition along width, and we declare
   *    failure
   */

  //
  // sort by coordinate and group TPCs; work with pointers for convenience
  //
  auto const TPCgroupInfo
    = sortAndGroupTPCsByRangeCoord<Range>(beginTPCwithArea, endTPCwithArea);
  std::vector<TPCwithArea_t const*> const& TPCs = TPCgroupInfo.first;
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
    TPCgroups = TPCgroupInfo.second;

  if (TPCs.empty()) return {}; // failure?

  //
  // for each group, create a subpartition
  //
  auto subparts = createSubpartitions
    (TPCgroups.cbegin(), TPCgroups.cend(), TPCs.cend(), subpartMaker);

  // if we have grouped everything in a single unit, we have not done any good
  if (subparts.size() == 1) return {};

  //
  // compute the total area (it might have been merged in a previous loop...)
  //
  auto totalArea = computeTotalArea(TPCs.cbegin(), TPCs.cend());

  //
  // construct and return the final partition
  //
  return std::make_unique<TPCPartitionResultType>
    (totalArea, std::move(subparts));

} // makeSortedPartition()


//------------------------------------------------------------------------------
template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makeWidthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  return makeSortedPartition<
    geo::part::WidthPartition<geo::TPCGeo const>,
    &geo::part::AreaOwner::Area_t::width
    >
    (beginTPCwithArea, endTPCwithArea, &makeDepthPartition<BeginIter, EndIter>);
} // makeWidthPartition()


template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makeDepthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  return makeSortedPartition<
    geo::part::DepthPartition<geo::TPCGeo const>,
    &geo::part::AreaOwner::Area_t::depth
    >
    (beginTPCwithArea, endTPCwithArea, &makeWidthPartition<BeginIter, EndIter>);
} // makeDepthPartition()


//------------------------------------------------------------------------------
template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makeGridPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  /*
   * Requires at least 4 input TPCs (otherwise, do not use GridPartition).
   *
   * 1. attempt a partition on width
   *    1. if failed, return failure
   * 2. attempt to partition the first subpartition on depth
   *    1. if failed, return failure
   * 3. extend each depth partition to the other width partitions
   *    1. if the extension of one partition line fails, discard it
   *    2. if no partition line survives, return failure
   * 4. run makePartition() on each of the cells
   * 5. create and return a GridPartition object from the subpartitions so
   *    created
   *
   * This algorithm could use some factorization...
   */
  using Area_t = geo::part::AreaOwner::Area_t;

  //
  // sort by width coordinate; work with pointers for convenience
  //
  auto const TPCgroupInfo = sortAndGroupTPCsByRangeCoord<&Area_t::width>
    (beginTPCwithArea, endTPCwithArea);
  std::vector<TPCwithArea_t const*> const& TPCs = TPCgroupInfo.first;
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
    TPCgroups = TPCgroupInfo.second;

  if (TPCs.empty()) return {}; // failure?
  // with only one TPC, then makeTPCPartitionElement() should be used instead!
  if (TPCs.size() < 4) return {};

  unsigned int const nWidthParts = TPCgroups.size();
  if (nWidthParts <= 1) return {}; // only one group ain't no good

  //
  // sort TPCs in the first width partition by depth coordinate
  //
  auto const FirstColGroupInfo
    = sortAndGroupTPCsByRangeCoord<&Area_t::depth>(TPCgroups[0], TPCgroups[1]);
  std::vector<TPCwithArea_t const*> const& FirstColTPCs
    = FirstColGroupInfo.first;
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
    FirstColGroups = FirstColGroupInfo.second;

  if (FirstColTPCs.empty()) return {}; // failure?
  if (FirstColGroups.size() <= 1 ) return {}; // only one row ain't good either

  //
  // collect all candidate separation ranges
  //
  // First depth partition has no lower limit, last one has no upper limit
  // (they include all TPCs with depth lower than the upper limit in the first
  // case, all TPCs with depth higher of the lower limit in the last case).
  // Checks need to be done in the gaps between depth partitions.
  // So we start by skipping the first border.

  lar::util::RealComparisons<double> coordIs(0.1);
  std::vector<Area_t::Range_t> depthGaps; // candidate gaps
  auto icnext = FirstColGroups.cbegin(), icprev = icnext,
    icend = FirstColGroups.cend();
  while (++icnext != icend) {
    //
    // establish the range of the group in the candidate depth partition
    // from the group in the first width partition
    //
    auto const cprev = *icprev;
    auto const cnext = *icnext;

    depthGaps.emplace_back
      ((*cprev)->area().depth.upper, (*cnext)->area().depth.lower);

    icprev = icnext;
  } // while
  assert(!depthGaps.empty());

  //
  // see that for every other width partition separations hold
  //
  auto igbegin = TPCgroups.cbegin();
  while (++igbegin != TPCgroups.cend()) {
    //
    // prepare the TPC groups within this width partition
    //
    auto igend = std::next(igbegin);
    auto gbegin = *igbegin;
    auto gend = (igend == TPCgroups.cend())? TPCs.cend(): *igend;

    auto const ColGroupInfo
      = sortAndGroupTPCsByRangeCoord<&Area_t::depth>(gbegin, gend);
    std::vector<TPCwithArea_t const*> const& ColTPCs = ColGroupInfo.first;
    std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
      ColGroups = ColGroupInfo.second;

    // failure to partition a single column means total failure
    if (ColTPCs.empty()) return {};
    if (ColGroups.size() <= 1) return {}; // only one row ain't good either

    //
    // compute the coverage of each of the depth groups
    //
    std::vector<TPCwithArea_t::Area_t::Range_t> groupDepths(ColGroups.size());
    auto iGDepth = groupDepths.begin();
    for (auto icgstart = ColGroups.cbegin(); icgstart != ColGroups.cend();
      ++icgstart, ++iGDepth)
    {
      auto const icgend = std::next(icgstart);
      auto ictpc = *icgstart;
      auto const ictend = (icgend == ColGroups.cend())? ColTPCs.cend(): *icgend;
      while (ictpc != ictend)
        iGDepth->extendToInclude((*(ictpc++))->area().depth);
    } // for

    //
    // check each of the remaining candidate gaps
    //
    auto iGap = depthGaps.begin();
    while (iGap != depthGaps.end()) {
      Area_t::Range_t& gap = *iGap;

      //
      // check that the gap holds
      //
      bool bGoodGap = false;
      // first TPC starting after the gap (even immediately after):
      auto iCGroup = std::lower_bound(
        groupDepths.cbegin(), groupDepths.cend(), gap.upper,
        SortTPCwithAreaByDepth()
        );

      // any TPCs before/after this gap?
      if ((iCGroup != groupDepths.begin()) && (iCGroup != groupDepths.end())) {
        Area_t::Range_t const& before = *(std::prev(iCGroup));
        Area_t::Range_t const& after = *iCGroup;
        Area_t::Range_t const TPCgap{ before.upper, after.lower };

        // correct the gap
        if (coordIs.strictlySmaller(iGap->lower, TPCgap.lower))
          iGap->lower = TPCgap.lower;
        if (coordIs.strictlyGreater(iGap->upper, TPCgap.upper))
          iGap->upper = TPCgap.upper;

        // if nothing is left, gap is gone
        bGoodGap = coordIs.nonSmaller(iGap->upper, iGap->lower);
      } // if TPCs around the gap

      //
      // if the gap has been flagged as bad, remove it
      //
      if (bGoodGap) ++iGap;
      else iGap = depthGaps.erase(iGap);

    } // while (separation)

    if (depthGaps.empty()) return {}; // no surviving gaps means failure

  } // while (width partition)

  //
  // turn the gaps into separators
  //
  std::vector<double> depthSep;
  std::transform(
    depthGaps.cbegin(), depthGaps.cend(), std::back_inserter(depthSep),
    [](auto const& r){ return (r.lower + r.upper) / 2.0; }
    );
  unsigned int const nDepthParts = depthSep.size() + 1;

  //
  // fill the groups with TPCs, and create subpartitions from each of them
  //
  geo::part::Partition<geo::TPCGeo const>::Subpartitions_t subparts
    (nWidthParts * nDepthParts);
  Area_t totalArea;

  unsigned int iWidth = 0;
  for (auto igbegin = TPCgroups.cbegin(); igbegin != TPCgroups.cend();
    ++igbegin, ++iWidth
  ) {

    // sort TPCs in this group (yes, again; this time we don't group just yet)
    auto igend = std::next(igbegin);
    auto gbegin = *igbegin;
    auto gend = (igend == TPCgroups.cend())? TPCs.cend(): *igend;
    std::vector<TPCwithArea_t const*> ColTPCs(gbegin, gend);
    std::sort(ColTPCs.begin(), ColTPCs.end(),
      SortTPCareaByAreaRangeLower<&Area_t::depth>());

    unsigned int iDepth = 0;
    auto cgstart = ColTPCs.cbegin(), TPCend = ColTPCs.cend();
    for (double sep: depthSep) {

      //
      // collect all TPCs for this partition
      //
      // the first TPC that starts *after* the separator:
      auto cgend
        = std::upper_bound(cgstart, TPCend, sep, SortTPCwithAreaByDepth());
      // if we cut out TPCs that were included because of some tolerance,
      // recover them now
      while (cgend != cgstart) {
        auto cglast = std::prev(cgend);
        if (coordIs.strictlySmaller((*cglast)->area().depth.lower, sep)) break;
        cgend = cglast;
      } // while
      assert(cgstart != cgend); // separator selection should guarantee this

      //
      // create and register the partition
      //
      auto part = makePartition(cgstart, cgend);
      if (!part) return {}; // late failure!
      totalArea.extendToInclude(part->area());
      subparts[iDepth * nWidthParts + iWidth] = std::move(part);

      ++iDepth;
      cgstart = cgend;
    } // for all depth separators

    //
    // collect all the TPCs after the last separator
    //
    auto part = makePartition(cgstart, TPCend);
    if (!part) return {}; // super-late failure!
    totalArea.extendToInclude(part->area());
    subparts[iDepth * nWidthParts + iWidth] = std::move(part);

  } // for all width partitions

  return std::make_unique<geo::part::GridPartition<geo::TPCGeo const>>
    (totalArea, std::move(subparts), nWidthParts, nDepthParts);

} // makeGridPartition()


//------------------------------------------------------------------------------
template <typename BeginIter, typename EndIter>
std::unique_ptr<geo::part::Partition<geo::TPCGeo const>> makePartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  /*
   * Organizes a list of TPCs in a hierarchical partition.
   * Three main elements are used:
   * - single element partition objects: that's the single TPC end point
   * - TPC groups organised along width
   * - TPC groups organised along depth
   *
   * The procedure is recursively analysing a set of TPCs:
   * - if the set is actually one TPC only, use a PartitionElement
   * - attempt partitioning o a grid; if fails:
   * - attempt partitioning along width:
   *    * determine overlapping groups: a set of TPCs which share part of the
   *      width range
   *    * recurse on each overlapping group with more than one TPC,
   *      attempting a depth partition
   *    * if that fails, bail out since we don't have code to deal with a layout
   *      with areas overlapping on both directions at the same time
   *    * add the single elements and the overlapping groups to the width
   *      partition
   * - attempt partitioning along height:
   *    * same algorithm as for width
   * - pick the partitioning with less elements
   */
  using value_type = std::remove_reference_t<decltype(*beginTPCwithArea)>;
  static_assert(
    std::is_pointer<value_type>()
    && std::is_same
      <std::decay_t<std::remove_pointer_t<value_type>>, TPCwithArea_t>(),
    "Iterators must point to TPCwithArea_t pointers."
    );


  auto const size = std::distance(beginTPCwithArea, endTPCwithArea);
  if (size == 1) {
    return makeTPCPartitionElement(**beginTPCwithArea);
  }

  auto gPart = makeGridPartition(beginTPCwithArea, endTPCwithArea);
  if (gPart) return gPart;

  auto wPart = makeWidthPartition(beginTPCwithArea, endTPCwithArea);
  auto dPart = makeDepthPartition(beginTPCwithArea, endTPCwithArea);

  if (wPart) {

    if (dPart) { // wPart && dPart
      if (wPart->nParts() < dPart->nParts()) return wPart;
      else                                   return dPart; // slight preference
    }
    else { // wPart && !dPart
      return wPart; // easy choice
    }

  }
  else {

    if (dPart) { // !wPart && dPart
      return dPart; // easy choice
    }
    else { // !wPart && !dPart
      return {}; // failure!!
    }

  }

} // makePartition(Iter)


//------------------------------------------------------------------------------
template <typename BeginIter, typename EndIter>
auto makeCPointerVector(BeginIter b, EndIter e) {
  using value_type = typename BeginIter::value_type;
  std::vector<value_type const*> result;
  result.reserve(std::distance(b, e));
  std::transform(b, e, std::back_inserter(result),
    [](auto& obj){ return std::addressof(obj); });
  return result;
} // makeCPointerVector()

template <typename T>
auto makeCPointerVector(std::vector<T> const& v)
  { return makeCPointerVector(v.cbegin(), v.cend()); }


//------------------------------------------------------------------------------
std::unique_ptr<geo::DriftPartitions::TPCPartition_t> makePartition
  (std::vector<TPCwithArea_t> const& TPCs)
{
  // TODO use range library instead:
//  auto TPCptrs = TPCs | ranges::view::transform(std::addressof);
  auto TPCptrs = makeCPointerVector(TPCs);
  using std::cbegin;
  using std::cend;
  return makePartition(cbegin(TPCptrs), cend(TPCptrs));
} // makePartition(coll)


//------------------------------------------------------------------------------
geo::DriftPartitions geo::buildDriftVolumes(geo::CryostatGeo const& cryo) {

  //
  // group TPCs by drift direction
  //
  auto TPCsByDriftDir = groupTPCsByDriftDir(cryo);

  //
  // determine the cryostat-wide drift direction (arbitrary but consistent)
  // and the decomposition base (using the same for all drift partitions);
  //

  // In practice we use the coordinate system from the first TPC;
  // we still check that all drift directions are compatible,
  // but the result of detection is ignored.
  /* auto globalDriftDir = */ detectGlobalDriftDir(keys(TPCsByDriftDir));
  using Direction_t = geo::DriftPartitions::Direction_t;
  geo::TPCGeo const& firstTPC = cryo.TPC(0);
  geo::DriftPartitions::Decomposer_t decomposer
    ({ cryo.Center(), firstTPC.RefWidthDir<Direction_t>(), firstTPC.RefDepthDir<Direction_t>() });

  //
  // further group TPCs by plane position in drift direction
  //
  std::vector<TPCgroup_t> TPCgroups;
  for (auto const& TPCsOnDriftDir: TPCsByDriftDir) {
    auto TPCs = sortTPCsByDriftCoord(TPCsOnDriftDir.second, decomposer);
    append(TPCgroups, groupByDriftCoordinate(TPCs));
  } // for

  //
  // verify coordinate system consistency between TPCs
  //
  for (auto const& TPCgroup: TPCgroups) {
    unsigned int errors = checkTPCcoords(TPCgroup.TPCs);
    if (errors > 0) {
      throw cet::exception("buildDriftVolumes")
        << "TPCs in partition have different drift directions ("
        << errors << " errors found in " << TPCgroup.TPCs.size() << " TPCs).\n";
    } // if
  } // for

  //
  // partition each group
  //
  geo::DriftPartitions partitions(decomposer);
  for (auto const& TPCgroup: TPCgroups) {
    auto TPCs = addAreaToTPCs(TPCgroup.TPCs, decomposer);
    auto part = makePartition(TPCs);
    if (!part) {
      cet::exception e("buildDriftVolumes");
      e << "Failed to construct partition out of " << TPCs.size() << " TPCs:";
      for (auto const& TPCinfo: TPCs) {
        e << "\n at " << TPCinfo.area() << " TPC ";
        TPCinfo.TPC->PrintTPCInfo(e, "   ", 5U);
      } // for
      throw e;
    } // if error
    partitions.addPartition(std::move(part));
  } // for

  return partitions;
} // geo::buildDriftVolumes()


//------------------------------------------------------------------------------
