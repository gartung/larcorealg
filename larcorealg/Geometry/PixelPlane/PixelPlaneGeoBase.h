/**
 * @file    larcorealg/Geometry/PixelPlane/PixelPlaneGeoBase.h
 * @brief   Base class for readout plane with pixels as sensitive elements.
 * @author  Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    October 23, 2019
 * @see     `larcorealg/Geometry/PixelPlane/PixelPlaneGeoBase.cxx`
 * @ingroup GeometryPixel
 */

#ifndef LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEOBASE_H
#define LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEOBASE_H

// LArSoft libraries
#include "larcorealg/Geometry/PixelPlane/PixelPlaneGeoInterface.h"
#include "larcorealg/Geometry/Decomposer.h" // geo::AffinePlaneBase

// ROOT libraries
#include <TGeoNode.h>

// C/C++ standard libraries
#include <vector>
#include <array>
#include <regex>
#include <utility> // std::move()
#include <optional>
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace geo { class PixelPlaneGeoBase; }

/**
 * @brief Base class for geometry description of a pixel readout plane.
 * 
 * This class implements the `geo:PlaneGeo` interface via
 * `geo::PixelPlaneGeoBase` base class, which should be referred to for any
 * detail on the class behaviour.
 * 
 * This class is an implementation of `geo::PixelPlaneGeoBase` which only adds
 * a specific initialization pattern to it.
 * 
 * 
 * Initialization
 * ===============
 * 
 * This class uses the standard pixel geometry initialization algorithm
 * `geo::PixelPlaneGeoBase::InitializePixelGeometry()`.
 * The initialization of the pixel geometry happens immediately at construction,
 * just after the base class initialization.
 * 
 * The pixel geometry information needed by the algorithm is obtained from the
 * GDML/ROOT description of the plane geometry.
 * The details are documented in `ExtractPixelGeometry()` method.
 * 
 */
class geo::PixelPlaneGeoBase: public geo::PixelPlaneGeoInterface {
  
    protected:
  
  /// Pixel initialization helper data structure.
  struct RectPixelGeometry_t {
    
    /// Information on a single pixel plane axis.
    struct AxisInfo_t {
      /// Direction of this side in local (GDML) plane coordinates;
      /// modulus is ignored.
      std::optional<LocalVector_t> dir;
      /// Length of the area covered by pixel along this direction [cm]
      std::optional<double> length;
      std::optional<unsigned int> nPixels; ///< Number of pixels.
      /// Size of the side of each pixel (i.e. the pitch) [cm]
      std::optional<double> pitch;
      
      /// Returns whether *all* the information is set.
      bool isComplete() const;
      
      /// Removes all information, reverting it to unset.
      void clear();

    }; // struct AxisInfo_t
    
    static constexpr std::size_t NSides = 2U;
    
    /// Information about one side of the pixels.
    std::array<AxisInfo_t, NSides> sides;
    
    std::optional<LocalPoint_t> center; ///< Center of the pixelized area.
    
    /// Returns whether all information is present.
    bool isComplete() const;
    
    /// Returns whether information on the axes is acceptable.
    bool checkAxes() const;
    
    /// Prints the information into a stream (no indent).
    template <typename Stream>
    void print(Stream&& out) const;

  }; // struct RectPixelGeometry_t
  
  
  /**
   * @brief Constructor: extracts pixel information from geometry description.
   * @param node GDML/ROOT object describing the plane itself
   * @param trans transformation describing position and orientation of the
   *              plane in the world (based on local-to-world transformation)
   * 
   * Pass-through constructor.
   */
  PixelPlaneGeoBase(TGeoNode const& node, geo::TransformationMatrix&& trans)
    : geo::PixelPlaneGeoInterface(node, std::move(trans)) {}
  
  
  // --- BEGIN --- Initialization customization for derived classes ------------
  /**
   * @name Initialization customization for derived classes
   *
   * The methods in this section provide some options for the initialization of
   * a pixel plane object.
   * None of this is used by default, and the derived classes will have to make
   * their choice. An example of such a choice is in `geo::PixelPlaneGeo`.
   * 
   * The method `InitializePixelGeometry()` fulfils the set up step of the
   * `geo::PixelPlaneGeoInterface` object that is required before the set up
   * after sorting (see `geo::PlaneGeo` for the details on the setup steps).
   * It does this part of set up based on the pixel geometry information
   * passed to it by argument (data structure of type `RectPixelGeometry_t`).
   * 
   * Two methods are also provided that can fill such pixel geometry
   * information: `ExtractPixelGeometry()` extracts it from the pixel volumes
   * in the GDML/ROOT detector description, while
   * `ReadPixelGeometryFromMetadata()` reads them from the auxiliary fields
   * in the GDML/ROOT geometry nodes.
   * 
   */
  /// @{
  
  /**
   * @brief Reads the pixel information from auxiliary data in the geometry.
   * @param startNode node of the plane geometry object
   * @param startValues information already present that will be used directly
   * @return the pixel geometry information as complete as possible
   * 
   * This method produces information about the placement and geometry of the
   * pixels on the plane, to be taken by
   * `geo::PixelPlaneGeoBase::InitializePixelGeometry()` and turned into actual
   * plane set up.
   * 
   * The algorithm attempts to extract only the information that is not already
   * set in `startValues`.
   * 
   * This algorithm looks for metadata describing the pixel grid structure.
   * The metadata is expected to be present in any of the volumes under the
   * plane volume. If multiple versions of the metadata are found, the
   * information read the earliest survives.
   * 
   * 
   * Supported metadata
   * -------------------
   * 
   * The following metadata elements are sought and interpreted:
   * 
   * * `pixelAdirection`, `pixelBdirection` (3D vectors in the frame of the
   *   anode plane volume) define the direction the rest of the
   *   information is assigned to; for lack of imagination, we call these two
   *   directions, or sides, "A" and "B"
   * * `pixelApitch`, `pixelBpitch` (length; unit required to be meter-based)
   *   the distance between consecutive pixels along one of the directions
   * * `pixelAnumber`, `pixelBnumber` (positive integral number)
   *   the number of pixels along one of the directions
   * * `pixelAsideLength`, `pixelBsideLength` (length; unit required to be
   *   meter-based) the length covered by the pixels along one of the directions
   * * `pixelActiveCenter` (3D point in the frame of the anode plane) fixes
   *   the _center_ of the rectangular area covered by pixels.
   * 
   * An example of a complete set of information:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.xml}
   * <auxiliary auxtype="pixelAdirection"   auxvalue="( 0, 0, 1 )"              />
   * <auxiliary auxtype="pixelApitch"       auxvalue="4.8421875"   auxunit="mm" />
   * <auxiliary auxtype="pixelAnumber"      auxvalue="64"                       />
   * <auxiliary auxtype="pixelAsideLength"  auxvalue="30.99"       auxunit="cm" />
   * <auxiliary auxtype="pixelBdirection"   auxvalue="( 0, 1, 0 )"              />
   * <auxiliary auxtype="pixelBpitch"       auxvalue="4.8421875"   auxunit="mm" />
   * <auxiliary auxtype="pixelBnumber"      auxvalue="64"                       />
   * <auxiliary auxtype="pixelBsideLength"  auxvalue="30.99"       auxunit="cm" />
   * <auxiliary auxtype="pixelActiveCenter" auxvalue="( 0, 0, 0 )" auxunit="cm" />
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * 
   * Adding metadata to the GDML description
   * ----------------------------------------
   * 
   * The way to add metadata to the GDML detector description so that LArSoft
   * can parse and utilize it are documented with the `geo::GeoMetadataParser`
   * class.
   * 
   */
  static RectPixelGeometry_t ReadPixelGeometryFromMetadata(
    TGeoNode const& startNode,
    RectPixelGeometry_t const& startValues = RectPixelGeometry_t{}
    );
  
  
  /**
   * @brief Extracts the pixel information.
   * @param startNode node of the plane geometry object
   * @param pixelNamePattern pattern used to recognize pixel nodes
   * @param startValues information already present that will be used directly
   * @return the pixel geometry information as complete as possible
   * 
   * This method produces information about the placement and geometry of the
   * pixels on the plane, to be taken by
   * `geo::PixelPlaneGeoBase::InitializePixelGeometry()` and turned into actual
   * plane set up.
   * 
   * The algorithm attempts to extract only the information that is not already
   * set in `startValues`.
   * 
   * The algorithm expects the pixels to be lying on a grid with orthogonal
   * directions, and that they are all at the same pitch on each of the
   * directions.
   * The position (center) of the pixels is extracted from the GDML/ROOT
   * geometry description: pixels are expected to be volumes somewhere under
   * the `startNode` and whose node names match the `pixelNamePattern` pattern
   * (`std::regex_search()`).
   * 
   * The algorithm will identify the geometric plane all the pixel lie upon,
   * then identify the pixels at the four corners of the covered rectangular
   * area. From the corners, the two grid directions are identified, together
   * with a reference corner pixel for convenience.
   * For each direction, the number of pixels is determined checking how many
   * pixels are aligned to the reference one along each of the directions.
   * The distance between pixels (pitch) is taken dividing the distance between
   * the corner pixels by the number of pixels in between them.
   * 
   * Finally, the center of the pixel area is taken as the middle point between
   * the four corners.
   */
  static RectPixelGeometry_t ExtractPixelGeometry(
    TGeoNode const& startNode, std::regex const& pixelNamePattern,
    RectPixelGeometry_t const& startValues = RectPixelGeometry_t{}
    );
  
  
  /**
   * @brief Uses `pixelGeometry` information for "wire" frame initialization.
   * @param pixelGeometry the pixel geometry information as complete as possible
   * 
   * This method is provided as standard initialization of pixel geometry as
   * part of `geo::PixelPlaneGeoInterface`. It is _not_ directly called by
   * `geo::PixelPlaneGeoInterface` initialization, and it is provided as an
   * option for the implementers of derived classes.
   * For some patterns of pixel geometry initialization, see the
   * @ref geo_PixelPlaneGeoInterface_Derivation "derived class initialization"
   * documentation.
   * 
   * 
   * Default implementation
   * -----------------------
   * 
   * The default implementation of this method does the following:
   *  * the pixel grid directions are set to match the plane frame width
   *    and depth direction
   *  * each of the two side information records in `pixelGeometry` is assigned
   *    to one of the two grid pixel directions, expecting the directions in
   *    `pixelGeometry` and the ones of the plane frame to be matching
   *  * for each side, the needed pixel geometry is extracted, with the
   *    information that is provided in `pixelGeometry` and trying to deduce
   *    any information that is not provided:
   *      * if the number of pixels is not provided, it is computed from the
   *        length of the side and the pitch, which are mandatory;
   *      * if the pitch is not provided, it is computed from the length of the
   *        side, which is mandatory; if the pitch and the length are both
   *        provided, the consistency of the pitch is checked but the provided
   *        value is used nevertheless
   *  * the origin of the decomposition frame is set exactly at the coordinate
   *    specified in `RectPixelGeometry_t::center`; if `pixelGeometry` does not
   *    specify the `center`, the origin is set so that the center of the area
   *    covered by the pixels is at the exact center of the plane box (including
   *    the drift coordinate being in the middle of that box). The origin is
   *    set to the center of the pixel with index #0, as prescribed by the
   *    protocol.
   * 
   * This function expects the information of the pixel geometry to be provided.
   * Two algorithms are provided in `geo::PixelPlaneGeoInterface` extract that
   * information: `ReadPixelGeometryFromMetadata()` reads the information from
   * the metadata in the geometry description, while `ExtractPixelGeometry()`
   * extracts it from the pixel volumes in the geometry description.
   * 
   * Whatever the mean, pixel geometry information needs to be provided in
   * `RectPixelGeometry_t` format, and the following information is *required*:
   *  * for each side, the direction of the side
   *  * for each side, at least two among the number of required pixels, their
   *    pitch or the length of the side of the pixel grid in their direction
   * 
   * All lengths are in centimeters.
   * 
   */
  void InitializePixelGeometry(RectPixelGeometry_t const& pixelGeometry);
  
  /// @}
  // --- END -- Initialization customization for derived classes ---------------
  
  
  
  // --- BEGIN --- Initialization procedure implementation ---------------------
  /// @name Initialization procedure implementation
  /// @{
  
  /// Returns an unsorted list of center of all pixels under `startNode`.
  static std::vector<LocalPoint_t> findPixelCenters
    (TGeoNode const& startNode, std::regex const& pixelNamePattern);
  
  /// Returns the most extreme points among the specified ones (clockwise).
  static std::array<LocalPoint_t, 4U> findCorners
    (std::vector<LocalPoint_t> const& points);
  
  /// Returns copies of the `points` which have component transverse to axis`
  /// no larger than `tol`.
  static std::vector<LocalPoint_t> findPointsOnAxis(
    LocalPoint_t const& ref, LocalVector_t const& axis,
    std::vector<LocalPoint_t> const& points,
    double const tol = 0.01
    );
  
  /// Creates an orthonormal base with origin at `o`.
  static geo::AffinePlaneBase<LocalVector_t, LocalPoint_t> makeBase
    (LocalPoint_t const& o, LocalVector_t const& a, LocalVector_t const& b);
  
  /// @}
  // --- END --- Initialization procedure implementation -----------------------
  
  /// Returns a string with human-readable translation of the content of `info`.
  static std::string DumpPixelGeometry(RectPixelGeometry_t const& info);
  
  
}; // class geo::PixelPlaneGeoBase


namespace geo {
  
  //----------------------------------------------------------------------------
  template <typename Stream>
  decltype(auto) operator<<
    (Stream&& out, geo::PixelPlaneGeoBase::RectPixelGeometry_t const& info)
    { info.print(out); return out; }
  
  //----------------------------------------------------------------------------
  
} // namespace geo


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Stream>
void geo::PixelPlaneGeoBase::RectPixelGeometry_t::print(Stream&& out) const
{
  
  out << " * center of the pixelized area: ";
  if (center) out << center.value();
  else        out << "n/a";
  
  out << "\n * sides:";
  for (auto const& side: sides) {
    
    out << "\n    - direction ";
    if (side.dir) out << side.dir->Unit();
    else          out << "unknown";
    out << ": ";
    
    if (side.nPixels) out << side.nPixels.value() << " ";
    
    if (side.pitch) out << side.pitch.value() << "-cm pixels";
    else            out << "pixels of unspecified size";
    
    out << " covering";
    if (side.length) out << " " << side.length.value() << " cm";
    else             out << " an unspecified length";
    
  } // for
  
} // geo::PixelPlaneGeoBase::RectPixelGeometry_t::print()


//------------------------------------------------------------------------------




#endif // LARCOREALG_GEOMETRY_PIXELPLANE_PIXELPLANEGEOBASE_H
