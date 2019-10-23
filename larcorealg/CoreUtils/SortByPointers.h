/**
 * @file   larcorealg/CoreUtils/SortByPointers.h
 * @brief  Silly utility to sort vectors indirectly.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   September 28, 2017
 *
 * This library is header-only.
 */

#ifndef LARCOREALG_COREUTILS_SORTBYPOINTER_H
#define LARCOREALG_COREUTILS_SORTBYPOINTER_H

// LArSoft libraries
#include "larcorealg/CoreUtils/makeValueIndex.h"
#include "larcorealg/CoreUtils/MetaUtils.h" // util::is_unique_ptr_v<>

// C/C++ standard libraries
#include <vector>
#include <algorithm> // std::transform(), std::sort()
#include <memory> // std::addressof()
#include <iterator> // std::back_inserter()
#include <utility> // std::move()
#include <type_traits> // std::add_pointer_t


namespace util {

  namespace details {

    template <typename Coll, typename PtrColl>
    struct MoveFromPointersImpl;

  } // namespace details


  //----------------------------------------------------------------------------
  /**
   * @brief Creates a STL vector with pointers to data from another collection.
   * @tparam Coll type of collection of data
   * @param coll data collection
   * @return a STL vector with pointers to `coll` data elements, with same order
   * 
   * The returned collection hosts pointers to elements with the same
   * constantness as the original collection (usually, non-const).
   */
  template <typename Coll>
  auto makePointerVector(Coll& coll);


  //----------------------------------------------------------------------------
  /**
   * @brief Moves the content from a collection of pointers to one of data.
   * @tparam Coll type of collection of data
   * @tparam PtrColl type of collection of pointers to data
   * @param dest collection to be filled
   * @param src collection with the pointers to data to be moved
   *
   * The data pointed from each pointer in `src` is moved into `dest`.
   * The destination collection is cleared first, and `Coll` must support both
   * `clear()` and `push_back()`
   *
   */
  template <typename Coll, typename PtrColl>
  void MoveFromPointers(Coll& dest, PtrColl& src)
    { details::MoveFromPointersImpl<Coll, PtrColl>::move(dest, src); }


  //----------------------------------------------------------------------------
  /**
   * @brief Applies sorting indirectly, minimizing data copy.
   * @tparam Coll type of collection to be sorted
   * @tparam Sorter type of sorter
   * @param coll collection to be sorted
   * @param sorter functor sorting a vector of pointers (`makePointerVector()`)
   *
   * The `sorter` functor can receive a reference to a vector as the one
   * produced by `makePointerVector(coll)` (that is, a C++ STL vector of
   * pointers to the value type of `Coll`), and sort it "in place".
   * The container `Comp` must implement `push_back()` call in a `std::vector`
   * fashion.
   *
   * The algorithm is equivalent to the following:
   * -# create a parallel vector of pointers to the data
   * -# sort the data pointers (delegating to `sorter`)
   * -# move the data, sorted, from the original collection to a new one
   * -# replace the content of cont with the one from the sorted collection
   *
   * Single elements are moved from the original collection to a new one.
   *
   * The data elements of `Coll` must be moveable, as `Coll` itself must be.
   *
   * @note Use this algorithm only as a last resort, as there are usually better
   *       ways to sort collections than this one, which is not even
   *       particularly optimized.
   *
   */
  template <typename Coll, typename Sorter>
  void SortByPointers(Coll& coll, Sorter sorter);


  //----------------------------------------------------------------------------
  /**
   * @brief Sorts a vector of unique pointers using a C pointer sorter.
   * @tparam Coll type of collection to be sorted
   * @tparam Sorter type of sorter function
   * @param coll collection to be sorted
   * @param sorter sorting procedure
   * 
   * This adapter moves the unique pointers around to match a sorted version of
   * source.
   * This is an expensive procedure, implying the creation of a temporary
   * vector and of additional supporting data: avoid it if at all possible.
   */
  template <typename Coll, typename Sorter>
  void SortUniquePointers(Coll& coll, Sorter&& sorter);
  
  
  //----------------------------------------------------------------------------

} // namespace util


//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
namespace util::details {
  
  //----------------------------------------------------------------------------
  template <typename Coll, typename = void>
  struct PointerVectorMaker {
    
    static auto make(Coll& coll) {
      
      using std::begin, std::end;
      
      using pointer_type = decltype(std::addressof(*begin(coll)));
      using ptr_coll_t = std::vector<pointer_type>;

      auto const n = coll.size();

      //
      // create the collection of pointers to data
      //
      ptr_coll_t ptrs;
      ptrs.reserve(n);
      std::transform(begin(coll), end(coll), std::back_inserter(ptrs),
        [](auto& obj){ return std::addressof(obj); });

      return ptrs;

    } // make()
    
  }; // struct PointerVectorMaker<>
  
  
  template <typename Coll>
  struct PointerVectorMaker
    <Coll, std::enable_if_t<util::is_unique_ptr_v<typename Coll::value_type>>>
  {
    
    
    static auto make(Coll& coll) {
      
      
      using coll_t = Coll;
      using unique_ptr_t = typename coll_t::value_type;
      using value_type = typename unique_ptr_t::element_type;
      using pointer_type = std::add_pointer_t<value_type>;
      using ptr_coll_t = std::vector<pointer_type>;
      
      static_assert(util::is_unique_ptr_v<unique_ptr_t>); // kind of silly now
      
      using std::size;
      auto const n = size(coll);

      //
      // create the collection of pointers to data
      //
      ptr_coll_t ptrs;
      ptrs.reserve(n);
      std::transform(coll.begin(), coll.end(), std::back_inserter(ptrs),
        [](auto& obj){ return obj.get(); });

      return ptrs;

    } // make()
    
  }; // struct PointerVectorMaker<unique_ptr>
  
  
  
  //----------------------------------------------------------------------------
  
  
} // namespace util::details



//------------------------------------------------------------------------------
template <typename Coll>
auto util::makePointerVector(Coll& coll)
  { return details::PointerVectorMaker<Coll>::make(coll); }


//------------------------------------------------------------------------------
template <typename Coll, typename Sorter>
void util::SortByPointers(Coll& coll, Sorter sorter) {

  using coll_t = Coll;

  //
  // create the collection of pointers to data
  //
  auto ptrs = makePointerVector(coll);

  //
  // delegate the sorting by pointers
  //
  sorter(ptrs);

  //
  // create a sorted collection moving the content from the original one
  //
  coll_t sorted;
  MoveFromPointers(sorted, ptrs);

  //
  // replace the old container with the new one
  //
  coll = std::move(sorted);

} // util::SortByPointers()


//------------------------------------------------------------------------------
template <typename Coll, typename Sorter>
void util::SortUniquePointers(Coll& coll, Sorter&& sorter) {

  using Collection_t = Coll;
  using UPtr_t = typename Collection_t::value_type;
  
  static_assert(util::is_unique_ptr_v<UPtr_t>);

  //
  // create the collection of pointers to data
  //
  auto ptrs = makePointerVector(coll);
  
  // data pointer -> index
  auto const ptrIndex = util::makeValueIndex(ptrs);
  
  //
  // delegate the sorting by pointers
  //
  sorter(ptrs);

  //
  // create a sorted collection moving the content from the original one
  //
  Collection_t sorted;
  for (auto const& dataPtr: ptrs) {
    std::size_t const originalIndex = ptrIndex.at(dataPtr);
    sorted.emplace_back(std::move(coll[originalIndex]));
  }
  
  //
  // replace the old container with the new one
  //
  coll = std::move(sorted);

} // util::SortUniquePointers()


//------------------------------------------------------------------------------
namespace util {
  namespace details {

    template <typename Coll, typename PtrColl>
    void moveFromPointersImplBase(Coll& dest, PtrColl& src)
      { for (auto&& ptr: src) dest.push_back(std::move(*ptr)); }


    template <typename Coll, typename PtrColl>
    struct MoveFromPointersImpl {
      static void move(Coll& dest, PtrColl& src)
        {
          dest.clear();
          moveFromPointersImplBase(dest, src);
        }
    }; // struct MoveFromPointersImpl


    template <typename Data, typename PtrColl>
    struct MoveFromPointersImpl<std::vector<Data>, PtrColl> {
      static void move(std::vector<Data>& dest, PtrColl& src)
        {
          dest.clear();
          dest.reserve(src.size());
          moveFromPointersImplBase(dest, src);
        }
    }; // struct MoveFromPointersImpl

  } // namespace details
} // namespace util


//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_SORTBYPOINTER_H
