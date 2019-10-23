/**
 * @file   NameSelector.cxx
 * @brief  A class providing a selection list: implementation file
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 24, 2015
 * @see    NameSelector.h
 */

// our library
#include "larcorealg/TestUtils/NameSelector.h"

// framework libraries
#include "canvas/Utilities/Exception.h"

#include <ostream>

namespace testing {

  //----------------------------------------------------------------------------
  //--- testing::NameSelector
  //----------------------------------------------------------------------------
  NameSelector::Name_t const NameSelector::DefaultName = "*";
  NameSelector::Name_t const NameSelector::ClearAllName = "!";

  //----------------------------------------------------------------------------
  void NameSelector::ParseName(Name_t name) {
    ProcessItem(known_names, name);
  } // NameSelector::ParseName()


  //----------------------------------------------------------------------------
  NameSelector::Response_t NameSelector::Query(Name_t name) const {
    ++(query_registry[name]);
    return LookupResponse(name);
  } // NameSelector::Query()


  //----------------------------------------------------------------------------
  bool NameSelector::Accepted(Name_t name) const {
    Response_t response = Query(name);
    if (response == rsThrow) {
      throw art::Exception(art::errors::Configuration)
        << "NameSelector: name '" << name << "' not configured.";
    }
    return response == rsAccepted;
  } // NameSelector::Accepted()


  //----------------------------------------------------------------------------
  void NameSelector::PrintConfiguration(std::ostream& out) const {

    // resort the known elements
    std::map<Response_t, Names_t> elements;
    for (KnownNames_t::value_type const& element: known_names)
      if (element.first != DefaultName)
        elements[element.second.response].insert(element.first);

    size_t nKnownElements = 0;
    if (!elements[rsAccepted].empty()) {
      auto const& selected_elements = elements[rsAccepted];
      Names_t::const_iterator iName = selected_elements.cbegin(),
        nend = selected_elements.cend();
      out << " accept " << selected_elements.size()
        << ": '" << *(iName++) << "'";
      while (iName != nend) out << ", '" << *(iName++) << "'";
      out << ";";
      nKnownElements += selected_elements.size();
    } // if accepting anything
    if (!elements[rsRejected].empty()) {
      auto const& selected_elements = elements[rsRejected];
      Names_t::const_iterator iName = selected_elements.cbegin(),
        nend = selected_elements.cend();
      out << " reject " << selected_elements.size()
        << ": '" << *(iName++) << "'";
      while (iName != nend) out << ", '" << *(iName++) << "'";
      out << ";";
      nKnownElements += selected_elements.size();
    } // if accepting anything
    if (!elements[rsThrow].empty()) {
      auto const& selected_elements = elements[rsThrow];
      Names_t::const_iterator iName = selected_elements.cbegin(),
        nend = selected_elements.cend();
      out << " throw on " << selected_elements.size()
        << ": '" << *(iName++) << "'";
      while (iName != nend) out << ", '" << *(iName++) << "'";
      out << ";";
      nKnownElements += selected_elements.size();
    } // if accepting anything
    if (nKnownElements > 0) out << " otherwise,";
    switch (DefaultResponse()) {
      case rsAccepted: out << " accept everything"; break;
      case rsRejected: out << " reject everything"; break;
      case rsThrow:    out << " throw on anything"; break;
      default:         out << " I don't know";
    } // switch
  } // NameSelector::PrintConfiguration()


  //----------------------------------------------------------------------------
  NameSelector::Response_t NameSelector::LookupResponse(Name_t name) const {
    KnownNames_t::const_iterator iResponse = known_names.find(name);
    return (iResponse == known_names.end())?
      DefaultResponse(): iResponse->second.response;
  } // NameSelector::LookupResponse()


  //----------------------------------------------------------------------------
  template <>
  void NameSelector::AddFirstName<>(KnownNames_t& name_set, Name_t name) {
    ProcessItem(name_set, name);
  } // NameSelector::AddFirstName<>()


  //----------------------------------------------------------------------------
  void NameSelector::InsertItem
    (KnownNames_t& name_set, Name_t item, Response_t response) const
  {
    name_set[item] = { response };
  } // NameSelector::InsertItem()


  //----------------------------------------------------------------------------
  void NameSelector::InsertItem
    (KnownNames_t& name_set, KnownNames_t::value_type item, Response_t response)
    const
  {
    // response is the instruction we have about how to add the item
    // item.second.response is what the default answer for that item is
    Response_t final_response; // uninitialized
    switch (response) {
      case rsAccepted:
        final_response = item.second.response; // respect the response
        break;
      case rsRejected: // flip the response
        switch (item.second.response) {
          case rsAccepted: final_response = rsRejected; break;
          case rsRejected: final_response = rsAccepted; break;
          default:
            throw art::Exception(art::errors::LogicError)
              << __func__ << ": unexpected code flow: invalid added response";
        } // switch item response
        break;
      default:
        throw art::Exception(art::errors::LogicError)
          << __func__ << ": unexpected code flow: invalid response";
    } // switch response
    InsertItem(name_set, item.first, final_response);
  } // NameSelector::InsertItem(KnownNames_t::value_type)


  //----------------------------------------------------------------------------
  void NameSelector::ProcessItem(KnownNames_t& name_set, Name_t item) const {
    // special: if this name is actually a directive to clear all
    if (item == ClearAllName) {
      ClearNameSet(name_set); // clear everything except the default
      return;
    } // if clear all
    Response_t response = ParseMode(item);
    Definitions_t::const_iterator iDefinition = FindDefinition(item);
    if (iDefinition == definitions.end()) InsertItem(name_set, item, response);
    else {
      for (KnownNames_t::value_type const& element: iDefinition->second)
        InsertItem(name_set, element, response);
    }
  } // NameSelector::ProcessItem()


  //----------------------------------------------------------------------------
  NameSelector::Definitions_t::const_iterator NameSelector::FindDefinition
    (Name_t& item) const
  {
    bool bForceDef = false;
    if (!item.empty()) {
      if (item[0] == '@') {
        bForceDef = true;
        item.erase(0, 1);
      }
    }
    Definitions_t::const_iterator iDefinition = definitions.find(item);
    if ((iDefinition == definitions.end()) && bForceDef) {
      throw art::Exception(art::errors::LogicError)
        << "no set named '" << item << "'\n";
    }
    return iDefinition;
  } // NameSelector::FindDefinition()


  //----------------------------------------------------------------------------
  void NameSelector::ClearNameSet(KnownNames_t& name_set) const {
    KnownNames_t::iterator iName = name_set.begin(), nend = name_set.end();
    while (iName != nend) {
      if (iName->first == DefaultName) ++iName;
      else iName = name_set.erase(iName);
    } // while
  } // NameSelector::ClearNameSet()


  //----------------------------------------------------------------------------
  NameSelector::Names_t NameSelector::QueriedWithStatus
    (Response_t answer) const
  {
    NameSelector::Names_t names;
    for (auto const& query_info: query_registry) {
      if (query_info.first == DefaultName) continue;
      if (LookupResponse(query_info.first) != answer) continue;
      names.insert(query_info.first);
    } // for
    return names;
  } // NameSelector::QueriedWithStatus()


  //----------------------------------------------------------------------------
  bool NameSelector::DoCheckQueryRegistry
    (std::ostream* out /* = nullptr */) const
  {
    Names_t missing;
    for (auto const& elem: known_names) {
      if (query_registry[elem.first] > 0U) continue;
      if (elem.first == DefaultName) continue;
      missing.insert(elem.first);
    } // for
    if (out && !missing.empty()) {
      (*out) << missing.size() << " items not queried:";
      for (Name_t const& name: missing) (*out) << " " << name;
      (*out) << std::endl;
    }
    return missing.empty();
  } // NameSelector::CheckQueryRegistry()


  //----------------------------------------------------------------------------
  NameSelector::Response_t NameSelector::ParseMode
    (Name_t& item, Response_t default_mode /* = rsAccepted */)
  {
    if (item[0] == '+') {
      item.erase(0, 1);
      return rsAccepted;
    }
    if (item[0] == '-') {
      item.erase(0, 1);
      return rsRejected;
    }
    return default_mode;
  } // NameSelector::ParseMode()


  //----------------------------------------------------------------------------

} // namespace testing
