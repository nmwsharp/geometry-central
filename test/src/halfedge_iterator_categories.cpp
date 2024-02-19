#include "geometrycentral/surface/halfedge_element_types.h"

#include "gtest/gtest.h"

#include <type_traits>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// ============================================================
// =============== Helper Type Traits
// ============================================================

// std::void_t is not in the standard library until c++17
template <typename... Ts>
using void_t = void;

// std::is_swappable is not in the standard library until c++17
template <typename T, typename = void>
struct is_swappable : std::false_type {};
template <typename T>
struct is_swappable<T, void_t<decltype(std::swap(std::declval<T>(), std::declval<T>()))>> : std::true_type {};

template <typename T, typename = void>
struct is_dereferencable : std::false_type {};
template <typename T>
struct is_dereferencable<T, void_t<decltype(*std::declval<T&>())>> : std::true_type {};

template <typename T, typename = void>
struct is_pre_incrementable : std::false_type {};
template <typename T>
struct is_pre_incrementable<T, void_t<decltype(++std::declval<T&>())>> : std::true_type {};

template <typename T, typename = void>
struct is_post_incrementable : std::false_type {};
template <typename T>
struct is_post_incrementable<T, void_t<decltype(std::declval<T&>()++)>> : std::true_type {};

// ============================================================
// =============== Templated Test Setup
// ============================================================

template <typename T>
class StlIteratorCategories : public ::testing::Test {};

using GCRangeTypes = ::testing::Types<
    HalfedgeSet, HalfedgeInteriorSet, HalfedgeExteriorSet, CornerSet, VertexSet, EdgeSet, FaceSet, BoundaryLoopSet,
    NavigationSetBase<VertexAdjacentVertexNavigator>, NavigationSetBase<VertexIncomingHalfedgeNavigator>,
    NavigationSetBase<VertexOutgoingHalfedgeNavigator>, NavigationSetBase<VertexAdjacentCornerNavigator>,
    NavigationSetBase<VertexAdjacentEdgeNavigator>, NavigationSetBase<VertexAdjacentFaceNavigator>,
    NavigationSetBase<EdgeAdjacentHalfedgeNavigator>, NavigationSetBase<EdgeAdjacentInteriorHalfedgeNavigator>,
    NavigationSetBase<EdgeAdjacentFaceNavigator>, NavigationSetBase<FaceAdjacentVertexNavigator>,
    NavigationSetBase<FaceAdjacentHalfedgeNavigator>, NavigationSetBase<FaceAdjacentCornerNavigator>,
    NavigationSetBase<FaceAdjacentEdgeNavigator>, NavigationSetBase<FaceAdjacentFaceNavigator>,
    NavigationSetBase<BoundaryLoopAdjacentVertexNavigator>, NavigationSetBase<BoundaryLoopAdjacentHalfedgeNavigator>,
    NavigationSetBase<BoundaryLoopAdjacentEdgeNavigator>>;
TYPED_TEST_SUITE(StlIteratorCategories, GCRangeTypes);


// ============================================================
// =============== Templated Test Cases
// ============================================================

TYPED_TEST(StlIteratorCategories, InputIteratorConcept) {
  using iterator_t = decltype(std::declval<TypeParam>().begin());

  // LegacyIterator requirements
  EXPECT_TRUE(std::is_copy_constructible<iterator_t>::value);
  EXPECT_TRUE(std::is_copy_assignable<iterator_t>::value);
  EXPECT_TRUE(std::is_destructible<iterator_t>::value);
  EXPECT_TRUE(is_swappable<iterator_t>::value);
  // Check for iterator_traits typedefs
  using value_type = typename std::iterator_traits<iterator_t>::value_type;
  using difference_type = typename std::iterator_traits<iterator_t>::difference_type;
  using pointer = typename std::iterator_traits<iterator_t>::pointer;
  using reference = typename std::iterator_traits<iterator_t>::reference;
  using iterator_category = typename std::iterator_traits<iterator_t>::iterator_category;
  EXPECT_TRUE(is_dereferencable<iterator_t>::value);
  EXPECT_TRUE(is_pre_incrementable<iterator_t>::value);

  // LegacyInputIterator
  // - operator-> is not checked for, the RangeIteratorBase and NavigationIteratorBase create their values on
  //   dereference, so returning a pointer to them is asking for a dangling pointer. operator-> is not required to meet
  //   the c++20 std::input_iterator concept so in practice not having this shouldnt affect the ability to use std
  //   algorithms.
  EXPECT_TRUE((std::is_convertible<decltype(std::declval<iterator_t>() != std::declval<iterator_t>()), bool>::value));
  EXPECT_TRUE((std::is_convertible<decltype(*std::declval<iterator_t>()), value_type>::value));
  EXPECT_TRUE(is_post_incrementable<iterator_t>::value);
  EXPECT_TRUE((std::is_convertible<decltype(*std::declval<iterator_t>()++), value_type>::value));
  EXPECT_TRUE((std::is_base_of<std::input_iterator_tag, iterator_category>::value));
}

TYPED_TEST(StlIteratorCategories, ForwardIteratorConcept) {
  using iterator_t = decltype(std::declval<TypeParam>().begin());

  using value_type = std::iterator_traits<iterator_t>::value_type;
  using reference = std::iterator_traits<iterator_t>::reference;
  using iterator_category = std::iterator_traits<iterator_t>::iterator_category;
  // LegacyForwardIterator requirements
  // Note: There are a number of LegacyForwardIterator requirements that we do not enforce here:
  // - RangeIteratorBase<F>::reference is not a reference type as the values are generated on dereference, returning a
  //   reference does not make sense
  // - Because of the above, the multipass guarantee is not strictly met. The value returned from dereferencing two
  //   iterators that compare equal will compare equal, but they will not be references to the same object. In practice
  //   this doesnt matter as far as using std algorithms
  EXPECT_TRUE(std::is_default_constructible<iterator_t>::value);
  EXPECT_TRUE((std::is_same<decltype(std::declval<iterator_t>()++), iterator_t>::value));
  EXPECT_TRUE((std::is_same<decltype(*std::declval<iterator_t>()++), reference>::value));
  EXPECT_TRUE((std::is_same<std::decay<reference>::type, value_type>::value));
  EXPECT_TRUE((std::is_base_of<std::forward_iterator_tag, iterator_category>::value));
}