/**
 * @file   larcorealg/Geometry/GeoNodePath.h
 * @brief  Class representing a path in ROOT geometry.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/Geometry/GeometryBuilder.h`,
 *         `larcorealg/Geometry/GeoNodePath.cxx`
 */

#ifndef LARCOREALG_GEOMETRY_GEONODEPATH_H
#define LARCOREALG_GEOMETRY_GEONODEPATH_H

// LArSoft libraries
#include "larcorealg/Geometry/LocalTransformation.h"

// ROOT libraries
#include "TGeoNode.h"

// C++ standard library
#include <vector>
#include <string>
#include <initializer_list>
#include <cstddef> // std::size_t


namespace geo {

  /**
   * @brief Representation of a node and its ancestry.
   *
   * A `GeoNodePath` contains a sequence of nodes, from the `root()` node down
   * to a `current()` one.
   *
   * It behaves like a `stack` in that it inserts and removes elements at the
   * "top", which is also what defines the current node.
   *
   */
  class GeoNodePath {

      public:

    // --- BEGIN Data types ----------------------------------------------------
    /// Type of node object.
    using Node_t = TGeoNode const;

    /// Type of list of nodes.
    using Nodes_t = std::vector<Node_t const*>;

    /// Type used to represent the depth of the path.
    using Depth_t = std::size_t;

    // --- END Data types ------------------------------------------------------

    // --- BEGIN Constructors and destructor -----------------------------------
    /// Default constructor: an empty path.
    GeoNodePath() = default;

    /// Sets all the the specified nodes into the current path.
    GeoNodePath(std::initializer_list<TGeoNode const*> nodes)
      : fNodes(nodes)
      {}

    /// Sets the nodes from `begin` to `end` as the path content.
    template <typename Iter>
    GeoNodePath(Iter begin, Iter end): fNodes(begin, end) {}

    // --- END Constructors and destructor -------------------------------------


    // --- BEGIN Query and access ----------------------------------------------
    /// Returns whether there is a current node.
    bool empty() const { return fNodes.empty(); }

    /// Returns the depth of the path (elements including up to the current).
    Depth_t depth() const { return fNodes.size(); }

    /// Returns the root node. Undefined if the path is empty.
    Node_t const& root() const { return *(fNodes.front()); }

    /// Returns the current node. Undefined if the path is empty.
    Node_t const& current() const { return *(fNodes.back()); }

    // --- END Query and access ------------------------------------------------


    // --- BEGIN Content management --------------------------------------------
    /// Adds a node to the current path.
    void append(Node_t const& node) { fNodes.push_back(&node); }

    /// Removes the current node from the path, moving the current one up.
    void pop() { fNodes.pop_back(); }
    // --- END Content management ----------------------------------------------


    /// Returns the total transformation to the current node, as a `Matrix`.
    template <typename Matrix = TGeoHMatrix>
    Matrix currentTransformation() const;

    /// Prints the full path (as node names) into a string.
    operator std::string() const;


      private:

    Nodes_t fNodes; ///< Local path of pointers to ROOT geometry nodes.

  }; // class GeoNodePath


} // namespace geo


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
template <typename Matrix /* = TGeoHMatrix */>
Matrix geo::GeoNodePath::currentTransformation() const {
  return geo::LocalTransformation<Matrix>::transformationFromPath
    (fNodes.begin(), fNodes.end());
} // geo::GeoNodePath::currentTransformation()


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_GEONODEPATH_H
