/**
 * @file   NameSelector.h
 * @brief  A class providing a selection list
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 24, 2015
 * @see    NameSelector.cxx
 */

#ifndef TEST_GEOMETRY_NAMESELECTOR_H
#define TEST_GEOMETRY_NAMESELECTOR_H 1

// C/C++ standard library
#include <initializer_list>
#include <iosfwd>
#include <string>
#include <set>
#include <map>
#include <stddef.h>

// ... and more in the template implementation section below

namespace testing {

  /**
   * @brief Manages a set of names
   *
   * This class contains a set of names to be "accepted" and a set of names
   * to be "rejected". Queries on unknown names will generate a default
   * answer that can be acceptance, rejection or an exception.
   *
   * The class is initialized by a list of strings following the convention
   * detailed in ProcessItem().
   * The default answer is set on construction, but it can be overridden to
   * a "reject by default" using a specific identifier.
   *
   */
  class NameSelector {

      public:
    using Name_t = std::string; ///< type representing a name
    using Names_t = std::set<Name_t>; ///< list of names

    /// Possible responses
    typedef enum {
      rsRejected,   ///< rejected
      rsAccepted,   ///< accepted
      rsThrow,      ///< throw art::Exception (art::errors::Configuration)
      rsDefault = rsAccepted
    } Response_t;


    using NameList = std::initializer_list<Name_t>;


    /// Constructor: an empty selector with a default answer for unknown names
    NameSelector(Response_t default_answer = rsDefault)
      { SetDefaultResponse(default_answer); }

    /// Constructor: Parse()s the specified items
    template <typename LIST>
    NameSelector(LIST const& items, Response_t default_answer = rsDefault)
      : NameSelector(default_answer)
      { Parse(items); }

    /// Sets the default answer for names that are not registered
    void SetDefaultResponse(Response_t default_answer)
      { known_names[DefaultName] = { default_answer }; }

    /**
     * @brief Parses the names in the list and adds them to the selector
     * @tparam LIST type of list of strings
     * @param items the definitions to be used
     *
     * Specified items are processed with ProcessItem() and the current set is
     * modified.
     */
    template <typename LIST>
    void Parse(LIST const& items) { BuildNameSet(known_names, items); }


    /**
     * @brief Parses a name and adds it to the selector
     * @param name name specification to be parsed
     *
     * The name is processed with ProcessItem() and the current set is modified.
     */
    void ParseName(Name_t name);

    /**
     * @brief Parses a list of names and adds them to the selector
     * @tparam NAMES Name_t-compatible types
     * @param names the definitions to be used
     *
     * Specified items are processed with ProcessItem() and modify the current
     * set.
     */
    template <typename... NAMES>
    void ParseNames(NAMES... names) { AddFirstName(known_names, names...); }

    /// Erases all the names in the selector (default answer is unchanged)
    void Clear() { ClearNameSet(known_names); }

    /**
     * @brief Defines a set
     * @tparam LIST type of list of strings
     * @param set_name name of the set to be defined
     * @param items the definitions to be used
     *
     * Specified items are processed with BuildSet(); the result is stored as
     * a set with the specified name.
     * If a set with that name already exists, it's replaced and the old one
     * is lost.
     */
    template <typename LIST>
    void Define(std::string set_name, LIST const& items);

    /**
     * @brief Parses a list of names and add them to the specified definition
     * @tparam NAMES Name_t-compatible types
     * @param set_name name of the set to be defined
     * @param names the definitions to be used
     *
     * Specified items are processed with ProcessItem(); the result is added to
     * the set with the specified name.
     * If a set with that name does not exist, it's created empty.
     */
    template <typename... NAMES>
    void AddToDefinition(std::string set_name, NAMES... names)
      { AddFirstName(definitions[set_name], names...);  }

    /// Returns the response for the specified name (does not throw)
    Response_t Query(Name_t name) const;

    /// Returns whether the name is accepted as good
    bool Accepted(Name_t name) const;

    /// Returns whether the name is rejected as bad
    bool Rejected(Name_t name) const { return !Accepted(name); }

    /// Returns the default answer for names that are not registered
    Response_t DefaultResponse() const
      { return known_names.find(DefaultName)->second.response; }

    /// Returns whether the name is accepted as good (alias for accept())
    bool operator() (Name_t name) const { return Accepted(name); }

    /// Prints the configuration into a stream
    void PrintConfiguration(std::ostream&) const;

    /// Returns a list of names that were accepted
    Names_t AcceptedNames() const { return QueriedWithStatus(rsAccepted); }

    /// Returns a list of names that were rejected
    Names_t RejectedNames() const { return QueriedWithStatus(rsRejected); }

    /// Reests the query registry
    void ClearQueryRegistry() { query_registry.clear(); }

    //@{
    /// Checks that no known element with valid response was left unqueried
    bool CheckQueryRegistry() const { return DoCheckQueryRegistry(); }
    bool CheckQueryRegistry(std::ostream& out) const
      { return DoCheckQueryRegistry(&out); }
    //@}


    static Name_t const DefaultName; ///< name representing the default
    static Name_t const ClearAllName; ///< name instructing to delete all names

      protected:

    /// A data structure containing how to react to a name
    typedef struct {
      Response_t response; ///< the response
    } NameResponse_t;

    /// Information about known names
    using KnownNames_t = std::map<Name_t, NameResponse_t>;

    /// Type of list of definitions
    using Definitions_t = std::map<Name_t, KnownNames_t>;

    /// Type of query counters
    using QueryRegistry_t = std::map<Name_t, size_t>;

    KnownNames_t known_names; ///< list of known names, with category

    Definitions_t definitions; ///< a set of definitions

    mutable QueryRegistry_t query_registry; ///< record of all the queries

    /// Returns the response for the specified name (does not register query)
    Response_t LookupResponse(Name_t name) const;

    /**
     * @brief Fills name_set with a list of items
     * @tparam LIST type of container of name directives
     * @param name_set set to be modified
     * @param items name directives
     *
     * A name set is modified according to the instructions in each of the
     * items. The items are parsed sequentially, and their order matters
     * for the final result.
     * Each item is processed through ProcessItem().
     */
    template <typename LIST>
    void BuildNameSet(KnownNames_t& name_set, LIST const& items) const;

    /// Parses the first of the provided names, and recurs
    template <typename... NAMES>
    void AddFirstName
      (KnownNames_t& name_set, Name_t name, NAMES... other_names);

    /// Adds an item to the name set, working in specified mode
    void InsertItem
      (KnownNames_t& name_set, Name_t item, Response_t response) const;

    /// Adds an item with response to the name set, working in specified mode
    void InsertItem(
      KnownNames_t& name_set, KnownNames_t::value_type item,
      Response_t response
      ) const;

    /**
     * @brief Fills name_set with an item
     * @param name_set set to be modified
     * @param item name directive
     *
     * A name set is modified according to the instruction in the items.
     * An identifier may represent either a single literal item or a set name;
     * it can appear in two forms:
     * - "\@identifier": always denotes a set named "identifier"
     * - "identifier": if a set named "identifier" exists, the name
     *   represents that set; otherwise, it represents a literal item
     *   named "identifier"
     *
     * The directives in items are:
     * - "identifier": the element or elements represented by the identifier
     *   replace the content of the current set
     * - "+identifier": the element or elements represented by the identifier
     *   are added to the current set; if the identifier is a predefined set,
     *   the elements accepted in the set are added as accepted here, the ones
     *   rejected are added as rejected here
     * - "-identifier": the element or elements represented by the identifier;
     *   if the identifier is a predefined set, the elements accepted in the
     *   set are added as rejected here, the ones rejected are added as
     *   accepted here
     *
     */
    void ProcessItem(KnownNames_t& name_set, Name_t item) const;

    /// Strips set specifier and returns iterator to the definition, or end()
    Definitions_t::const_iterator FindDefinition(Name_t& item) const;

    /// Erases all the names in the selector (default answer is unchanged)
    void ClearNameSet(KnownNames_t& name_set) const;

    /// Returns the list of queried names whose response is answer
    Names_t QueriedWithStatus(Response_t answer) const;

    // Performs the actual registry check, optionally printing a message
    bool DoCheckQueryRegistry(std::ostream* out = nullptr) const;

    /// Strips the mode specifier from item and returns the insertion mode
    static Response_t ParseMode
      (Name_t& item, Response_t default_answer = rsAccepted);

  }; // class NameSelector

} // namespace testing


//------------------------------------------------------------------------------
//--- Template implementation
//---

// C/C++ standard library
#include <utility> // std::move()

//------------------------------------------------------------------------------
template <typename LIST>
void testing::NameSelector::Define(std::string set_name, LIST const& items) {
  KnownNames_t name_set;
  BuildNameSet(name_set, items);
  definitions[set_name] = std::move(name_set);
} // testing::NameSelector::Define()


//------------------------------------------------------------------------------
template <typename LIST>
void testing::NameSelector::BuildNameSet
  (KnownNames_t& name_set, LIST const& items) const
{
  for (Name_t item: items) ProcessItem(name_set, item);
} // testing::NameSelector::BuildNameSet()


//------------------------------------------------------------------------------
namespace testing {
  // forward declaration of specialization
  template <>
  void NameSelector::AddFirstName<>(KnownNames_t& name_set, Name_t name);
} // namespace testing

template <typename... NAMES>
void testing::NameSelector::AddFirstName
  (KnownNames_t& name_set, Name_t name, NAMES... other_names)
{
  AddFirstName(name_set, name);
  AddFirstName(name_set, other_names...);
} // testing::NameSelector::AddFirstName()


//------------------------------------------------------------------------------

#endif // TEST_GEOMETRY_NAMESELECTOR_H
