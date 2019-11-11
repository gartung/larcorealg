/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeoBase.cxx
 * @brief   Base class for readout plane with pixels as sensitive elements.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 23, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeoBase.h`
 * @ingroup GeometryPixel
 */

// class header
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeoBase.h"

// LArSoft includes
#include "larcorealg/Geometry/GeoMetadataParser.h"
#include "larcorealg/Geometry/GeoNodePath.h"
#include "larcorealg/Geometry/LocalTransformation.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::middlePoint()..
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/zip.h"
#include "larcorealg/CoreUtils/counter.h"
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// support libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoExtension.h"
#include "TMap.h"

// C/C++ standard library
#include <ostream>
#include <vector>
#include <array>
#include <algorithm> // std::swap()
#include <regex>
#include <optional>
#include <limits> // std::numeric_limits<>
#include <cmath> // std::round()
#include <cassert>


// -----------------------------------------------------------------------------
namespace {
  
  // ---------------------------------------------------------------------------
  class PixelCenterFinder {
    
    using Path_t = geo::GeoNodePath; ///< Type of node path.
    
    std::regex const fPattern; ///< The pattern to recognize a pixel volume.
    unsigned int fMaxDepth; ///< Do not descend below this deepness.
    
    /// Fills `centers` descending into `path`.
    template <typename Point>
    void fillPixelCenters(std::vector<Point>& centers, Path_t& path) const;
    
      public:
    
    /// Constructor: initializes the regular expression to recognize a pixel.
    PixelCenterFinder
      (std::regex const& namePattern, unsigned int maxDepth = 10U)
      : fPattern(namePattern), fMaxDepth(maxDepth) {}
    
    /// Returns whether the current node of the `path` is a pixel.
    bool isPixel(TGeoNode const& node) const
      { return std::regex_match(node.GetName(), fPattern); }
    
    /// Returns a list of points starting with the specified node.
    template <typename Point = geo::Point_t>
    std::vector<Point> findPixelCenters(TGeoNode const& startNode) const;
    
  }; // class PixelCenterFinder
  
  // ---------------------------------------------------------------------------
  template <typename MetaDataObj>
  class GeoMetadataCollector {
      public:
    using Metadata_t = MetaDataObj;
    using MetadataPtr_t = Metadata_t const*;
    using MetadataColl_t = std::vector<MetadataPtr_t>;
    
    MetadataColl_t collectFrom(TGeoNode const& node) const
      { MetadataColl_t coll; collectMetadata(coll, node); return coll; }
    MetadataColl_t operator() (TGeoNode const& node) const
      { return collectFrom(node); }
    
    MetadataColl_t collectFrom(TGeoVolume const& node) const
      { MetadataColl_t coll; collectMetadata(coll, node); return coll; }
    MetadataColl_t operator() (TGeoVolume const& volume) const
      { return collectFrom(volume); }
    
      private:
    
    void collectMetadata
      (MetadataColl_t& metadata, TGeoVolume const& startVolume) const
      {
        TMap const* nodeMetadata = getMetadata(startVolume);
        if (nodeMetadata) metadata.push_back(nodeMetadata);
        
        for (auto iNode: util::counter(startVolume.GetNdaughters())) {
          TGeoNode const* node = startVolume.GetNode(iNode);
          if (node) collectMetadata(metadata, *node);
        }
      } // collectMetadata(TGeoVolume)

    void collectMetadata
      (MetadataColl_t& metadata, TGeoNode const& startNode) const
      {
        TMap const* nodeMetadata = getMetadata(startNode);
        if (nodeMetadata) metadata.push_back(nodeMetadata);
        TGeoVolume const* volume = startNode.GetVolume();
        if (volume) collectMetadata(metadata, *volume);
      } // collectMetadata(TGeoNode)
    
    template <typename NodeOrVolume>
    static MetadataPtr_t getMetadata(NodeOrVolume const& node)
      {
        //
        // we expect the node/volume to have a user extension of type 
        // TGeoRCExtension which holds a TMap. That's what we're after.
        //
        auto pExt = dynamic_cast<TGeoRCExtension const*>(node.GetUserExtension());
        return pExt? dynamic_cast<TMap const*>(pExt->GetUserObject()): nullptr;
      } // getMetadata()

  }; // class GeoMetadataCollector<>
  
  
  // ---------------------------------------------------------------------------
  template <typename T>
  bool sameOptional(std::optional<T> const& a, std::optional<T> const& b)
    { return a? (a.value() == b.value()): !b; }


  // ---------------------------------------------------------------------------
  
} // local namespace



// -----------------------------------------------------------------------------
template <typename Point /* = geo::Point_t */>
std::vector<Point> PixelCenterFinder::findPixelCenters
  (TGeoNode const& node) const
{
  std::vector<Point> centers;
  Path_t path { &node };
  fillPixelCenters(centers, path);
  return centers;
} // PixelCenterFinder::findPixelCenters()


// -----------------------------------------------------------------------------
template <typename Point>
void PixelCenterFinder::fillPixelCenters
  (std::vector<Point>& centers, Path_t& path) const
{
  //
  // if this is a target object, we are set
  //
  if (isPixel(path.current())) {
    centers.push_back(
      geo::LocalTransformation
        (path.currentTransformation<geo::TransformationMatrix>())
        .LocalToWorld(geo::origin<Point>())
      );
    return;
  }

  //
  // descend into the next layer down, concatenate the results and return them
  //
  if (path.depth() >= fMaxDepth) return;

  TGeoVolume const& volume = *(path.current().GetVolume());
  for (auto i: util::counter<int>(volume.GetNdaughters())) {
    path.append(*(volume.GetNode(i)));
    fillPixelCenters(centers, path);
    path.pop();
  } // for
  
} // PixelCenterFinder::fillPixelCenters()


// -----------------------------------------------------------------------------
// --- geo::PixelPlaneGeoBase::RectPixelGeometry_t
// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoBase::RectPixelGeometry_t::AxisInfo_t::isComplete() const
{
  return dir && length && nPixels && pitch;
} // geo::PixelPlaneGeoBase::RectPixelGeometry_t::AxisInfo_t::isComplete()


// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoBase::RectPixelGeometry_t::AxisInfo_t::clear() {
  dir.reset();
  length.reset();
  nPixels.reset();
  pitch.reset();
} // geo::PixelPlaneGeoBase::RectPixelGeometry_t::AxisInfo_t::clear()


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoBase::RectPixelGeometry_t::isComplete() const {
  
  for (auto const& side: sides) if (!side.isComplete()) return false;
  return center.has_value();
  
} // geo::PixelPlaneGeoBase::RectPixelGeometry_t::isComplete()


// -----------------------------------------------------------------------------
bool geo::PixelPlaneGeoBase::RectPixelGeometry_t::checkAxes() const {
  /*
   * If both axis directions are not set, their content must be the same.
   */
  if (sides[0U].dir || sides[1U].dir) return true;
  
  // neither direction is set: content of the two axes must match
  return sameOptional(sides[0U].length, sides[1U].length)
    && sameOptional(sides[0U].nPixels, sides[1U].nPixels)
    && sameOptional(sides[0U].pitch, sides[1U].pitch)
    ;
  
} // geo::PixelPlaneGeoBase::RectPixelGeometry_t::checkAxes()


// -----------------------------------------------------------------------------
std::ostream& geo::operator<<
  (std::ostream& out, geo::PixelPlaneGeoBase::RectPixelGeometry_t const& info)
  { info.print(out); return out; }


// -----------------------------------------------------------------------------
// --- geo::PixelPlaneGeoBase
// -----------------------------------------------------------------------------
#if 0
std::regex const geo::PixelPlaneGeoBase::DefaultPixelPattern {
  ".*pixel.*",
  std::regex::basic | std::regex::icase | std::regex::optimize
  };


// -----------------------------------------------------------------------------
#endif // 0
// -----------------------------------------------------------------------------
void geo::PixelPlaneGeoBase::InitializePixelGeometry
  (RectPixelGeometry_t const& pixelGeometry)
{
  
  /*
   * Assigns the origin, pitch, number and direction of the pixels and grid.
   * This is a elaboration of the raw information from the geometry, passed
   * via `pixelGeometry` argument.
   * 
   * The center and directions from `pixelGeometry` are only partially
   * respected.
   * The directions are used to associate the pitches to the right axis, but
   * the axes are imposed to be the ones of the plane frame.
   * The center vector is supposed to point to the center of the pixelated area.
   * The origin will be moved from there to the center of the first pixel,
   * which is going to be located at one of the four corners of the pixelated
   * area.
   * 
   * Requires:
   *  * local-world transformations be available
   *  * frame geometry to be set (exact origin depth coordinate does not matter)
   * 
   */
  
  /*
   * NOTE this is part of the pixel plane initialization procedure documented
   *      in `geo::PixelPlaneGeoInterface` class. If changing *what* is being
   *      initialized or its *order*, please also update that documentation at
   *      the top of `geo::PixelPlaneGeoInterface` class Doxygen documentation (in
   *      `larcorealg/Geometry/PixelPlane/PixelPlaneGeoInterface.h`, section
   *      "Initialization steps").
   */
  
  /*
   * REMINDER: do not rely on virtual methods from derived classes here, as
   *           they might not be available yet (this method is called in class
   *           constructor, when thee derived class constructor hasn't been run
   *           yet)
   */
  
  constexpr lar::util::RealComparisons cmp(1e-2); // 0.1 mm tolerance
  constexpr auto vcmp = lar::util::makeVector3DComparison(cmp);
  
  assert(vcmp.nonZero(fDecompFrame.MainDir()));
  assert(vcmp.nonZero(fDecompFrame.SecondaryDir()));
  
  using namespace geo::pixel;
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "PixelPlaneGeoInterface::InitializePixelGeometry(): "
    << "Information received from plane geometry:"
    << "\n" << pixelGeometry;
  
  
  //
  // 1. assign the correct axis information to each direction;
  //    no nonsense here: axes must follow the frame one way or the other way
  //    (just two ways are allowed: either axis aligned with `WidthDir()`)
  //
  struct ExtendedAxisInfo_t {
    RectPixelGeometry_t::AxisInfo_t basicInfo;
    geo::Vector_t axisDir; ///< Direction and size in the world frame.
  }; // struct ExtendedAxisInfo_t
  
  assert(pixelGeometry.sides.size() == NCoords);
  
  std::array<ExtendedAxisInfo_t, NCoords> axes;
  for (auto&& [ side, axis ]: util::zip(pixelGeometry.sides, axes)) {
    if (!side.dir || vcmp.zero(side.dir.value())) { // ideally, norm should be 1
      throw cet::exception("PixelPlaneGeoInterface")
        << "Specification of pixel plane axis was received as "
        << side.dir.value()
        << ".\n";
    }
    axis = { side, toWorldCoords(side.dir.value()) };
  } // for
  
  //
  // make sure that the data on `ixSec` is aligned with `WidthDir()`:
  //
  if (vcmp.parallel(axes[ixSec].axisDir, WidthDir())) {
    // all good already
  }
  else if (vcmp.parallel(axes[ixSec].axisDir, DepthDir())) {
    std::swap(axes[ixSec], axes[ixMain]);
  }
  else {
    throw cet::exception("PixelPlaneGeoInterface")
      << "InitializePixelGeometry(): pixel axis system { "
      << axes[ixMain].axisDir << " x " << axes[ixSec].axisDir
      << " } is misaligned with the plane frame system { "
      << WidthDir<geo::Vector_t>() << " x " << DepthDir<geo::Vector_t>() << " }"
      "\n";
  }
  
  fDecompPixel.SetMainDir
    (geo::vect::rounded01(axes[ixMain].axisDir.Unit(), 1e-5));
  fDecompPixel.SetSecondaryDir
    (geo::vect::rounded01(axes[ixSec].axisDir.Unit(), 1e-5));
  
  //
  // 2. set the pixel geometry on each axis
  //
  
  for (auto&& [ iAxis, axisInfo, pitch, nPixels ]
    : util::enumerate(axes, fPitches, fNPixels ))
  {
    // purpose of each loop is to set a value for `pitch` and `nPixels`
    
    DirIndex_t const ixDir { iAxis };
    RectPixelGeometry_t::AxisInfo_t const& rawInfo = axisInfo.basicInfo;
    
    //
    // We have as input side length, number of pixels and wire pitch:
    // we may be underconstrained or overconstrained here...
    // In case of overconstraining, we trust number of pixels, side length and
    // pitch in decreasing order.
    //
    if (rawInfo.nPixels) { // if given number of pixels...
      nPixels = rawInfo.nPixels.value();
      if (nPixels == 0U) {
        throw cet::exception("PixelPlaneGeoInterface")
          << "Number of pixel is provided for "
          << getDirectionName(ixDir) << " axis[" << rawInfo.dir.value()
          << "], and it is provided wrong! (" << nPixels << ")\n"
          << pixelGeometry << "\n"
          ;
      } // if number of pixels is 0
      
      if (rawInfo.length) { // if given number of pixels, total length...
        pitch = rawInfo.length.value() / nPixels; // should we round here?
        if (rawInfo.pitch) { // if given pixels, total length and pitch
          // just check for consistency, but eventually use the provided value
          if (
            lar::util::RealComparisons(1e-2) // 10 um tolerance
              .nonEqual(pitch, rawInfo.pitch.value())
            )
          {
            MF_LOG_WARNING("Geometry") << "Information from raw geometry for "
              << getDirectionName(ixDir) << " [" << rawInfo.dir.value()
              << "] pixel axis is inconsistent, yielding a pitch of " << pitch
              << " cm (retained) while expecting it to be "
              << rawInfo.pitch.value() << " cm.\nFull geometry information:\n"
              << pixelGeometry << "\n";
          } // if inconsistency detected
          pitch = rawInfo.pitch.value();
        } // if overconstrained
      } // if
      else if (rawInfo.pitch) { // if given number of pixels and pitch
        // well, since we do not save total length anyway we are pretty much set
        pitch = rawInfo.pitch.value();
      }
      else { // we are given only pixel number!
        // should we supply with the direction size here instead?
        throw cet::exception("PixelPlaneGeoInterface")
          << "We are not given enough information on "
          << getDirectionName(ixDir) << " axis[" << rawInfo.dir.value()
          << "]: only the number of pixels!\n"
          << pixelGeometry << "\n"
          ;
      } // if ... else
      
    }
    else if (rawInfo.length && rawInfo.pitch) { // if given length and pitch
      
      pitch = rawInfo.pitch.value();
      
      // we *truncate* so that the specified length is not exceeded
      // (with some tolerance)
      
      nPixels
        = static_cast<unsigned int>(std::round(rawInfo.length.value() / pitch));
      double const length = nPixels * pitch;
      if (cmp.strictlyGreater(length, rawInfo.length.value())) {
        --nPixels;
      }
      assert(cmp.nonGreater(nPixels * pitch, rawInfo.length.value()));
    }
    else {
      throw cet::exception("PixelPlaneGeoInterface")
        << "We are not given enough information on "
        << getDirectionName(ixDir) << " axis [" << rawInfo.dir.value() << "]:\n"
        << pixelGeometry << "\n"
        ;
    }
    
    // final sanity checks
    if ((nPixels == 0U) || cmp.nonPositive(pitch)) {
      throw cet::exception("PixelPlaneGeoInterface")
        << "Failure in pixel geometry determination algorithm for "
        << getDirectionName(ixDir) << " axis [" << rawInfo.dir.value()
        << "] yielded "
        << nPixels << " pixels with " << pitch << " cm pitch. Geometry info:\n"
        << pixelGeometry << "\n";
    }
    
    MF_LOG_TRACE("Geometry")
      << "PixelPlaneGeoInterface::InitializePixelGeometry(): "
      << getDirectionName(ixDir) << " grid axis: "
      << getSensElemDir(ixDir) << " axis covering "
      << getSensElemDirSize(ixDir) << " cm with "
      << getNsensElem(ixDir) << " pixels of side "
      << getSensElemPitch(ixDir) << " cm each"
      ;
    
  } // for
  
  //
  // 3. now set the position of the pixel plane;
  //    that is driven by the origin of the pixel decomposition frame,
  //    which corresponds to the position of the first pixel:
  //    that is what we are pursuing now.
  //    We do not bother with the position along thickness direction,
  //    which stays the same as from the geometry
  //    (which is not that bad a choice);
  //    note that `fromCenterToFirstPixel()` requires some of the quantities
  //    just set (pretty much all of them, in fact)
  //
  LocalPoint_t const center
    = pixelGeometry.center.value_or(geo::origin<LocalPoint_t>());
  fDecompPixel.SetReferencePoint
    (fromCenterToFirstPixel(toWorldCoords(center)));
  
  MF_LOG_TRACE("PixelPlaneGeoInterface")
    << "geo::PixelPlaneGeoInterface::InitializePixelGeometry():"
      " directions set to:"
    << "\n * main: " << fDecompPixel.MainDir() << " after geometry dir "
      << axes[ixMain].axisDir << " with " << fNPixels[ixMain]
      << " pixels with " << fPitches[ixMain] << " cm pitch"
    << "\n * secondary: " << fDecompPixel.SecondaryDir()
      << " after geometry dir " << axes[ixSec].axisDir << " with "
      << fNPixels[ixSec] << " pixels with " << fPitches[ixSec] << " cm pitch"
    << "\n * origin: " << fDecompPixel.ReferencePoint()
      << " (center of pixelized area: " << toWorldCoords(center)
      << ")"
    ;
  
  
} // geo::PixelPlaneGeoInterface::InitializePixelGeometry()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoBase::ReadPixelGeometryFromMetadata(
  TGeoNode const& startNode,
  RectPixelGeometry_t const& startValues /* = {} */
  )
  -> RectPixelGeometry_t
{
  using namespace std::string_literals;
  
  // BUG the double brace syntax is required to work around clang bug 21629
  static constexpr std::array<const char*, RectPixelGeometry_t::NSides> SideTags
//    = { "A", "B" };
    = {{ "A", "B" }};
  
  std::vector<TMap const*> metadataColl
    = GeoMetadataCollector<TMap>().collectFrom(startNode);
  
  MF_LOG_TRACE("geo::PixelPlaneGeoBase")
    << "ReadPixelGeometryFromMetadata(): collected " << metadataColl.size()
    << " metadata sets from node [" << startNode.GetName() << "]";
  
  RectPixelGeometry_t info = startValues;
  for (TMap const* metadata: metadataColl) {
    
    geo::GeoMetadataParser const fetcher(*metadata);
    
    for (auto&& [ side, sideTag ]: util::zip(info.sides, SideTags)) {
      
      /*
       * Welcome to The Corner of Amazing C++.
       * So the structure binding above does not introduce new variables,
       * meaning that "side" and "sideTag" are not new variables.
       * Also, lambdas can capture only variables. So, "sideTag" can't be
       * captured! Which is something C++ people acknowledge to be undesirable.
       * Well, we can create a new lambda-local variable, "sideTag", and assign
       * to it as initialization value the one of the non-variable "sideTag".
       * Then the latter is not captured, but its value is used for a different
       * thing that is not exactly a capture but in the end looks very much like
       * one. If C++ will ever solve this issue, we can just then capture
       * "sideTag" directly and be happy with it.
       * Note that GCC 7.3 bit it, although 8.0 seems to reject it, as does
       * Clang 5.
       */
      // prefixed keys become something like "pixelAdirection" etc.
      auto prefix
        = [sideTag=sideTag](std::string const& key)
        { return "pixel"s + sideTag + key; };
      
      //
      // "pixelDdirection"
      //
      if (!side.dir) {
        side.dir = fetcher.getVector<LocalVector_t>(prefix("direction"));
        if (side.dir) {
          MF_LOG_TRACE("geo::PixelPlaneGeoBase")
            << "ReadPixelGeometryFromMetadata(): found side " << sideTag
            << " direction: " << side.dir.value();
        } // dir
      } // if we haven't this information yet
      
      //
      // "pixelDpitch"
      //
      if (!side.pitch) {
        side.pitch = fetcher.getValue(prefix("pitch"), "m");
        if (side.pitch) {
          side.pitch.value() *= 100.0; // [m] => [cm]
          MF_LOG_TRACE("geo::PixelPlaneGeoBase")
            << "ReadPixelGeometryFromMetadata(): found side " << sideTag
            << " pitch: " << (side.pitch.value() * 10.0) << " mm";
        } // pitch
      } // if we haven't this information yet
      
      //
      // "pixelDnumber"
      //
      if (!side.nPixels) {
        side.nPixels = fetcher.getValue<unsigned int>(prefix("number"));
        if (side.nPixels) {
          MF_LOG_TRACE("geo::PixelPlaneGeoBase")
            << "ReadPixelGeometryFromMetadata(): found side " << sideTag
            << " pixels: " << side.nPixels.value();
        } // nPixels
      } // if we haven't this information yet
      
      //
      // "pixelDsideLength"
      //
      if (!side.length) {
        side.length = fetcher.getValue(prefix("sideLength"), "m");
        if (side.length) {
          side.length.value() *= 100.0; // [m] => [cm]
          MF_LOG_TRACE("geo::PixelPlaneGeoBase")
            << "ReadPixelGeometryFromMetadata(): found side " << sideTag
            << " length: " << side.length.value() << " cm";
        } // length
      } // if we haven't this information yet
      
    } // for sides
    
    //
    // "pixelActiveCenter"
    //
    if (!info.center) {
      info.center = fetcher.getVector<LocalPoint_t>("pixelActiveCenter", "m");
      if (info.center) {
        info.center.value() *= 100.0; // [m] => [cm]
        MF_LOG_TRACE("geo::PixelPlaneGeoBase")
          << "ReadPixelGeometryFromMetadata(): found pixel plane center: "
          << info.center.value() << " cm";
      } // center
    } // if we haven't this information yet
    
    
  } // for metadata objects
  
  // 
  // parameter check: do not allow information without its direction
  // 
  for (auto&& [ side, sideTag ]: util::zip(info.sides, SideTags)) {
    if ((side.pitch || side.nPixels || side.length) && !side.dir) {
      throw cet::exception("geo::PixelPlaneGeoBase")
        << "Metadata is incomplete on side " << sideTag
        << ", which misses its direction:\n" << info << "\n";
    } // if has info
  } // for sides
  
  return info;
} // geo::PixelPlaneGeoBase::ReadPixelGeometryFromMetadata()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoBase::ExtractPixelGeometry(
  TGeoNode const& startNode, std::regex const& pixelNamePattern,
  RectPixelGeometry_t const& startValues /* = {} */
  )
  -> RectPixelGeometry_t
{
  
  /*
   * We assume that pixels are accommodated into a grid lying on a geometric
   * 2D plane which is "parallel" to the wire plane object, that is, the normal
   * of the former and the normal of `geo::PlaneGeo`s frame base match.
   * We also assume pixels not to have gaps in between, and to be evenly spaced
   * in each of the two grid axes, although the spacing may be different between
   * the two axes.
   * 
   * Goals are:
   *  * identify the two axes
   *  * determine the number of pixels on each axis
   *  * determine the extension of the grid on each axis
   *  * determine the center of the grid
   * 
   * 1. collect all pixel centers in the local plane frame
   * 2. identify the four corners
   * 3. identify the axes; for each axis:
   * 3.1. identify the number of pixels
   * 3.2. find the pixel pitch
   * 3.3. quantify the extent of the plane
   * 4. determine the center of the pixel grid
   * 
   * Requirements:
   *  * local-to-world transformation already defined
   *  * GDML/ROOT plane volume and descendants available
   *  * frame base already defined: `Width()` `Depth()`, and `WidthDir()` and
   *    `DepthDir()` (their verse is not relevant)
   * 
   */
  
  constexpr std::size_t NAxes = 2U;
  
  constexpr auto vcmp = lar::util::makeVector3DComparison(1e-4); // 1 um tol.
  
  if (!startValues.checkAxes()) {
    throw cet::exception("geo::PixelPlaneGeoBase")
      << "ExtractPixelGeometry(): start values are not consistent.\n"
      << startValues << "\n";
  }
  
  // shortcut: if we already have all the information, we won't change it anyway
  if (startValues.isComplete()) return startValues;
  
  // to extract any of the information, we do need the complete set of pixels
  // (actually, at least the four corners, but then the pitch and number of
  //  pixel should better be already there)
  
  //
  // 1. collect all pixel centers in the local plane frame
  //
  // We descend into subvolumes.
  //
  
  std::vector<LocalPoint_t> const centers
    = findPixelCenters(startNode, pixelNamePattern);
  if (centers.size() < 4U) {
    throw cet::exception("geo::PixelPlaneGeoBase")
      << "Only " << centers.size() << " pixels found under node '"
      << startNode.GetName() << "' (at least 4 are required)\n";
  } // if
  
  //
  // 2. identify the four corners
  //
  
  std::array<LocalPoint_t, 4U> const corners = findCorners(centers);
  
  //
  // 3. identify the axes
  //
  // These axes do not need to have to do with the plane frame base axes.
  //
  
  LocalPoint_t const& refPoint = corners[0U];
  // BUG the double brace syntax is required to work around clang bug 21629
//   std::array<LocalVector_t, NAxes> const localAxes = { // rounding: 10 um
  std::array<LocalVector_t, NAxes> const localAxes = {{ // rounding: 10 um
    geo::vect::rounded01(corners[1U] - refPoint, 1e-3),
    geo::vect::rounded01(corners[3U] - refPoint, 1e-3)
//   };
  }};
  
  
  unsigned int TotalPixels = 1U;
  bool bNpixelAutodetection = true;
  // we deal with starting values of axes separately later:
  RectPixelGeometry_t info = startValues;
  for (auto& axisInfo: info.sides) axisInfo.clear();
  
  // C++17 note: the tuple { axis, axisInfo } is constant, but the elements
  // are not necessarily so.
  for (auto const& [ axis, axisInfo ]: util::zip(localAxes, info.sides)) {
    
    LocalVector_t const axisDir = axis.Unit();
    
    // if we have start values,
    // we need to figure out which are the ones for this axis
    for (auto const& startInfo: startValues.sides) {
      //
      // if this startInfo has a direction that is ours, or if it has none and
      // therefore it's some default, then we use it;
      // also, if this startInfo has a direction that is ours, we are done
      //
      
      if (startInfo.dir) { // we have info to consider
        
        // this is not information for our axis: ignore it
        if (!vcmp.parallel(startInfo.dir.value(), axisDir)) continue;
        
        // full direction match here, we are good to go
        axisInfo = startInfo;
        break;
        
      }
      
      // this looks like default info: good for us, pending better luck
      axisInfo = startInfo; // we keep looking for a match though
    } // for
    
    // if we ended up with some default info, direction is not set and we do now
    if (!axisInfo.dir) axisInfo.dir = axisDir;
    
    //
    // 4. identify the number of pixels
    //
    
    if (axisInfo.nPixels) {
      bNpixelAutodetection = false;
    }
    else {
      std::vector<LocalPoint_t> const pointsOnAxis
        = findPointsOnAxis(refPoint, axisDir, centers, 0.01); // tenth of mm
      
      if (pointsOnAxis.size() < 2U) {
        throw cet::exception("geo::PixelPlaneGeoBase")
          << "Found only " << pointsOnAxis.size()
          << " on direction " << axis << " (we require at least 2)\n"
          ;
      } // if
      
      axisInfo.nPixels = pointsOnAxis.size();
      
    } // if extract number of pixels
    unsigned int const nPixels = axisInfo.nPixels.value();
    
    TotalPixels *= nPixels;
    
    // 
    // 3.2. find the pixel pitch
    // 
    
    if (!axisInfo.pitch) {
      // all rounding issues are left unchecked so far...
      axisInfo.pitch = axis.R() / (nPixels - 1U);
    }
    
    // 
    // 3.3. quantify the extent of the plane
    // 
    // We use the corner-to-corner size, plus one pitch.
    // Maybe round?
    // 
    
    if (!axisInfo.length) {
      // we leave the axis length untouched at this time
    //  axisInfo.length = axisInfo.pitch.value() * nPixels;
    }
    
  } // for axes
  
  if (TotalPixels != centers.size()) {
    if (bNpixelAutodetection) {
      throw cet::exception("geo::PixelPlaneGeoBase")
        << "Logic error: deduced " << info.sides[0U].nPixels.value() << "x"
        << info.sides[1U].nPixels.value() << " = " << TotalPixels
        << " pixels (" << centers.size() << " expected)\n";
    }
    else {
      //
      // This warning means that we have found a number of pixel volumes
      // in the geometry description different from what the number of pixels
      // we thing by now to be on the plane. This may be fine if the number
      // of pixel was given as input instead of us figuring it out of the
      // pixel volumes. For example, if the geometry contains the minimal number
      // of pixel volumes (one per corner) and then specifies how many pixels
      // there are via metadata, the `TotalPixels` we find is `4` but the real
      // number is from the metadata.
      //
      MF_LOG_DEBUG("geo::PixelPlaneGeoBase")
        << "We deduced " << info.sides[0U].nPixels.value() << "x"
        << info.sides[1U].nPixels.value() << " = " << TotalPixels
        << " pixels but expected " << centers.size() << " from initial requests"
        ;
    } // see if we care...
  } // if pixels are not as many as expected
  
  // 
  // 4. determine the center of the pixel grid
  // 
  // This is easily achieved by averaging the four corners.
  //
  if (!info.center) {
    info.center = geo::vect::middlePoint(corners.begin(), corners.end());
  }
  
  return info;
  
} // geo::PixelPlaneGeoBase::ExtractPixelGeometry()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoBase::findPixelCenters
  (TGeoNode const& startNode, std::regex const& pixelNamePattern)
  -> std::vector<LocalPoint_t>
{
  PixelCenterFinder pixelFinder(pixelNamePattern);
  
  std::vector<LocalPoint_t> centers
    = pixelFinder.findPixelCenters<LocalPoint_t>(startNode);
  
  MF_LOG_TRACE("geo::PixelPlaneGeoBase")
    << "Found " << centers.size() << " pixel volumes in '"
    << startNode.GetName() << "'"
    ;
  
  return centers;
} // geo::PixelPlaneGeoBase::findPixelCenters()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoBase::findCorners
  (std::vector<LocalPoint_t> const& points)
  -> std::array<LocalPoint_t, 4U>
{
  // 
  // returned points are in clockwise order
  // 
  
  lar::util::RealComparisons cmp(1e-4); // this is 1 um
  
  assert(points.size() >= 3);
  geo::AffinePlaneBase<LocalVector_t, LocalPoint_t> const decompBase = makeBase(
    points.front(),
    points[1U] - points.front(),
    points.back() - points.front()
    );
  if (cmp.nonPositive(decompBase.NormalDir().Mag2())) {
    throw cet::exception("geo::PixelPlaneGeoBase")
      << "Logic error: could not make a base out of points #0 "
      << points.front() << ", #1 " << points[1U] << " and #"
      << (points.size() - 1U) << " " << points.front() << "\n";
  } // if
  
  using limits = std::numeric_limits<double>;
  struct RecordPoint_t {
    double m = 0.0;
    double s = 0.0;
    LocalPoint_t point;
  }; // struct RecordPoint_t
  
  // we need four extremes, and we also need to decide what to do when multiple
  // points have the same extreme in one coordinate; and we have to choose that
  // so that all four corners are found; only two ways available... we pick one:
  RecordPoint_t max_m{ limits::lowest(), limits::max   () }; // max m (or min s)
  RecordPoint_t max_s{ limits::lowest(), limits::lowest() }; // max s (or max m)
  RecordPoint_t min_m{ limits::max   (), limits::lowest() }; // min m (or max s)
  RecordPoint_t min_s{ limits::max   (), limits::max   () }; // min s (or min m)
  
  for (LocalPoint_t const& point: points) {
    
    double const m = decompBase.MainDir().Dot(decompBase.ToVector(point));
    double const s = decompBase.SecondaryDir().Dot(decompBase.ToVector(point));
    
    // rounded comparisons matter only in the main extreme coordinate;
    // if that coordinate is close to the extreme and the other one is close
    // too, then the two points are matching within errors
    // (which we might decide to check against)
    
    // max m (or min s)
    if (
      cmp.strictlyGreater(m, max_m.m)
      || (cmp.equal(m, max_m.m) && (s < max_m.s))
      )
    {
      max_m = { m, s, point };
    }
    
    // max s (or max m)
    if (
      cmp.strictlyGreater(s, max_s.s)
      || (cmp.equal(s, max_s.s) && (m > max_s.m))
      )
    {
      max_s = { m, s, point };
    }
    
    // min m (or max s)
    if (
      cmp.strictlySmaller(m, min_m.m)
      || (cmp.equal(m, min_m.m) && (s > min_m.s))
      )
    {
      min_m = { m, s, point };
    }
    
    // min s (or min m)
    if (
      cmp.strictlySmaller(s, min_s.s)
      || (cmp.equal(s, min_s.s) && (m < min_s.m))
      )
    {
      min_s = { m, s, point };
    }
    
  } // for
  
  // BUG the double brace syntax is required to work around clang bug 21629
//  return { max_m.point, max_s.point, min_m.point, min_s.point };
  return {{ max_m.point, max_s.point, min_m.point, min_s.point }};
} // geo::PixelPlaneGeoBase::findCorners()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoBase::findPointsOnAxis(
  LocalPoint_t const& ref, LocalVector_t const& axis,
  std::vector<LocalPoint_t> const& points,
  double const tol /* = 0.01 */
  )
  -> std::vector<LocalPoint_t>
{
  lar::util::RealComparisons cmp(tol*tol); // will compare with squares
  assert(cmp.equal(axis.Mag2(), 1.0));
  
  std::vector<LocalPoint_t> onAxis;
  for (LocalPoint_t const& point: points) {
    
    if (cmp.strictlyPositive((point - ref).Cross(axis).Mag2()))
      continue; // transverse component too large
    
    onAxis.push_back(point);
  } // for
  
  return onAxis;
} // geo::PixelPlaneGeoBase::findPointsOnAxis()


// -----------------------------------------------------------------------------
auto geo::PixelPlaneGeoBase::makeBase
  (LocalPoint_t const& o, LocalVector_t const& a, LocalVector_t const& b)
  -> geo::AffinePlaneBase<LocalVector_t, LocalPoint_t>
{
  
  LocalVector_t const n = a.Cross(b);
  
  LocalPoint_t const origin = o;
  LocalVector_t const main = a;
  LocalVector_t const sec = n.Cross(a);
  
  return { origin, main.Unit(), sec.Unit() };
  
} // geo::PixelPlaneGeoBase::makeBase()


// -----------------------------------------------------------------------------
