#ifndef CARPETX_ARITH_ARR_HXX
#define CARPETX_ARITH_ARR_HXX

#include "defs.hxx"

#include <array>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <tuple>

namespace Arith {

template <typename T, std::size_t N> class arr {

  T elts[N];

public:
  // Member types

  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type &;
  using const_reference = const value_type &;
  using pointer = value_type *;
  using const_pointer = const value_type *;

  class iterator {
    friend class arr;

  public:
    using difference_type = arr::difference_type;
    using value_type = arr::value_type;
    using pointer = arr::pointer;
    using reference = arr::reference;
    using iterator_category = std::random_access_iterator_tag;

  private:
    pointer ptr;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator(pointer ptr)
        : ptr(ptr) {}

  public:
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator() = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator(const iterator &) =
        default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator(iterator &&) =
        default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator &
    operator=(const iterator &) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator &
    operator=(iterator &&) = default;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
    swap(iterator &other) noexcept {
      auto tmp = ptr;
      ptr = other.ptr;
      other.ptr = tmp;
    }

    // LegacyIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference operator*() const {
      return *ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator &operator++() {
      ++ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator operator++(int) {
      auto tmp = *this;
      ++*this;
      return tmp;
    }

    // LegacyInputIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST pointer operator->() const {
      return &**this;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator==(const iterator &other) const {
      return ptr == other.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator!=(const iterator &other) const {
      return !(*this == other);
    }

    // LegacyForwardIterator

    // LegacyBidirectionalIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator &operator--() {
      --ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator operator--(int) {
      auto tmp = *this;
      --*this;
      return tmp;
    }

    // LegacyRandomAccessIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator &
    operator+=(difference_type n) {
      ptr += n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator
    operator+(iterator iter, difference_type n) {
      return iter += n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator
    operator+(difference_type n, iterator iter) {
      return iter += n;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator &
    operator-=(difference_type n) {
      ptr -= n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator
    operator-(iterator iter, difference_type n) {
      return iter -= n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST difference_type
    operator-(const iterator &iter, const iterator &other) {
      return iter.ptr - other.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference
    operator[](difference_type n) const {
      return *(*this + n);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<(const iterator &iter, const iterator &other) {
      return (iter - other) < 0;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>(const iterator &iter, const iterator &other) {
      return other < iter;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>=(const iterator &iter, const iterator &other) {
      return !(iter < other);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<=(const iterator &iter, const iterator &other) {
      return !(iter > other);
    }
  };

  class const_iterator {
    friend class arr;

  public:
    using difference_type = arr::difference_type;
    using value_type = arr::value_type;
    using pointer = arr::const_pointer;
    using reference = arr::const_reference;
    using iterator_category = std::random_access_iterator_tag;

  private:
    pointer ptr;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator(pointer ptr)
        : ptr(ptr) {}

  public:
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator() = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
    const_iterator(const const_iterator &) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
    const_iterator(const_iterator &&) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator &
    operator=(const const_iterator &) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator &
    operator=(const_iterator &&) = default;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
    swap(const_iterator &other) noexcept {
      auto tmp = ptr;
      ptr = other.ptr;
      other.ptr = tmp;
    }

    // LegacyIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference operator*() const {
      return *ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator &
    operator++() {
      ++ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
    operator++(int) {
      auto tmp = *this;
      ++*this;
      return tmp;
    }

    // LegacyInputIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST pointer operator->() const {
      return &**this;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator==(const const_iterator &other) const {
      return ptr == other.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator!=(const const_iterator &other) const {
      return !(*this == other);
    }

    // LegacyForwardIterator

    // LegacyBidirectionalIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator &
    operator--() {
      --ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
    operator--(int) {
      auto tmp = *this;
      --*this;
      return tmp;
    }

    // LegacyRandomAccessIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator &
    operator+=(difference_type n) {
      ptr += n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
    operator+(const_iterator iter, difference_type n) {
      return iter += n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
    operator+(difference_type n, const_iterator iter) {
      return iter += n;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator &
    operator-=(difference_type n) {
      ptr -= n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
    operator-(const_iterator iter, difference_type n) {
      return iter -= n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST difference_type
    operator-(const const_iterator &iter, const const_iterator &other) {
      return iter.ptr - other.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference
    operator[](difference_type n) const {
      return *(*this + n);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<(const const_iterator &iter, const const_iterator &other) {
      return (iter - other) < 0;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>(const const_iterator &iter, const const_iterator &other) {
      return other < iter;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>=(const const_iterator &iter, const const_iterator &other) {
      return !(iter < other);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<=(const const_iterator &iter, const const_iterator &other) {
      return !(iter > other);
    }
  };

  class reverse_iterator {
    friend class arr;

  public:
    using difference_type = arr::difference_type;
    using value_type = arr::value_type;
    using pointer = arr::pointer;
    using reference = arr::reference;
    using iterator_category = std::random_access_iterator_tag;

  private:
    pointer ptr;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator(pointer ptr)
        : ptr(ptr) {}

  public:
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator() = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
    reverse_iterator(const reverse_iterator &) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
    reverse_iterator(reverse_iterator &&) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator &
    operator=(const reverse_iterator &) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator &
    operator=(reverse_iterator &&) = default;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
    swap(reverse_iterator &other) noexcept {
      auto tmp = ptr;
      ptr = other.ptr;
      other.ptr = tmp;
    }

    // LegacyIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference operator*() const {
      return *(ptr - 1);
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator &
    operator++() {
      --ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator
    operator++(int) {
      auto tmp = *this;
      ++*this;
      return tmp;
    }

    // LegacyInputIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST pointer operator->() const {
      return &**this;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator==(const reverse_iterator &other) const {
      return ptr == other.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator!=(const reverse_iterator &other) const {
      return !(*this == other);
    }

    // LegacyForwardIterator

    // LegacyBidirectionalIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator &
    operator--() {
      ++ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator
    operator--(int) {
      auto tmp = *this;
      --*this;
      return tmp;
    }

    // LegacyRandomAccessIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator &
    operator+=(difference_type n) {
      ptr -= n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator
    operator+(reverse_iterator iter, difference_type n) {
      return iter += n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator
    operator+(difference_type n, reverse_iterator iter) {
      return iter += n;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator &
    operator-=(difference_type n) {
      ptr += n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator
    operator-(reverse_iterator iter, difference_type n) {
      return iter -= n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST difference_type
    operator-(const reverse_iterator &iter, const reverse_iterator &other) {
      return other.ptr - iter.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference
    operator[](difference_type n) const {
      return *(*this + n);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<(const reverse_iterator &iter, const reverse_iterator &other) {
      return (iter - other) < 0;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>(const reverse_iterator &iter, const reverse_iterator &other) {
      return other < iter;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>=(const reverse_iterator &iter, const reverse_iterator &other) {
      return !(iter < other);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<=(const reverse_iterator &iter, const reverse_iterator &other) {
      return !(iter > other);
    }
  };

  class const_reverse_iterator {
    friend class arr;

  public:
    using difference_type = arr::difference_type;
    using value_type = arr::value_type;
    using pointer = arr::const_pointer;
    using reference = arr::const_reference;
    using iterator_category = std::random_access_iterator_tag;

  private:
    pointer ptr;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
    const_reverse_iterator(pointer ptr)
        : ptr(ptr) {}

  public:
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator() =
        default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
    const_reverse_iterator(const const_reverse_iterator &) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
    const_reverse_iterator(const_reverse_iterator &&) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator &
    operator=(const const_reverse_iterator &) = default;
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator &
    operator=(const_reverse_iterator &&) = default;

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
    swap(const_reverse_iterator &other) noexcept {
      auto tmp = ptr;
      ptr = other.ptr;
      other.ptr = tmp;
    }

    // LegacyIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference operator*() const {
      return *(ptr - 1);
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator &
    operator++() {
      --ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
    operator++(int) {
      auto tmp = *this;
      ++*this;
      return tmp;
    }

    // LegacyInputIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST pointer operator->() const {
      return &**this;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator==(const const_reverse_iterator &other) const {
      return ptr == other.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator!=(const const_reverse_iterator &other) const {
      return !(*this == other);
    }

    // LegacyForwardIterator

    // LegacyBidirectionalIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator &
    operator--() {
      ++ptr;
      return *this;
    }
    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
    operator--(int) {
      auto tmp = *this;
      --*this;
      return tmp;
    }

    // LegacyRandomAccessIterator

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator &
    operator+=(difference_type n) {
      ptr -= n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
    operator+(const_reverse_iterator iter, difference_type n) {
      return iter += n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
    operator+(difference_type n, const_reverse_iterator iter) {
      return iter += n;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator &
    operator-=(difference_type n) {
      ptr += n;
      return *this;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
    operator-(const_reverse_iterator iter, difference_type n) {
      return iter -= n;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST difference_type
    operator-(const const_reverse_iterator &iter,
              const const_reverse_iterator &other) {
      return other.ptr - iter.ptr;
    }

    constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference
    operator[](difference_type n) const {
      return *(*this + n);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<(const const_reverse_iterator &iter,
              const const_reverse_iterator &other) {
      return (iter - other) < 0;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>(const const_reverse_iterator &iter,
              const const_reverse_iterator &other) {
      return other < iter;
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator>=(const const_reverse_iterator &iter,
               const const_reverse_iterator &other) {
      return !(iter < other);
    }

    friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
    operator<=(const const_reverse_iterator &iter,
               const const_reverse_iterator &other) {
      return !(iter > other);
    }
  };

  // Member functions

  // Constructors

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr() : elts{} {}
  // constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  // arr(std::initializer_list<T> x) {
  //   assert(x.size() == N);
  //   for (std::size_t i = 0; i < N; ++i)
  //     elts[i] = x.begin()[i];
  // }
private:
  template <std::size_t... Is>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr(std::initializer_list<T> x,
                                                     std::index_sequence<Is...>)
      : elts{x.begin()[Is]...} {
    static_assert(sizeof...(Is) == N);
  }

public:
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr(std::initializer_list<T> x)
      : arr(x, (assert(x.size() == N), std::make_index_sequence<N>())) {}

  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  arr(const std::array<T, N> &x)
      : elts{} {
    for (std::size_t i = 0; i < N; ++i)
      elts[i] = x.data()[i];
  }
  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  operator std::array<T, N>() const noexcept {
    std::array<T, N> x;
    for (std::size_t i = 0; i < N; ++i)
      x.data()[i] = elts[i];
    return x;
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr(const arr &) = default;
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr(arr &&) = default;

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr &
  operator=(const arr &) = default;
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr &
  operator=(arr &&) = default;

  // Element access

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference at(size_type pos) {
    if (CCTK_BUILTIN_EXPECT(pos >= N, false))
#if defined __CUDACC__ || defined __HIPCC__
      // asm("trap;");
      assert(0);
#else
      throw std::out_of_range("arr");
#endif
    return elts[pos];
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reference
  at(size_type pos) const {
    if (CCTK_BUILTIN_EXPECT(pos >= N, false))
#if defined __CUDACC__ || defined __HIPCC__
      // asm("trap;");
      assert(0);
#else
      throw std::out_of_range("arr");
#endif
    return elts[pos];
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference
  operator[](size_type pos) {
    return elts[pos];
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reference
  operator[](size_type pos) const {
    return elts[pos];
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference front() {
    return *begin();
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reference front() const {
    return *cbegin();
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reference back() {
    return *rbegin();
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reference back() const {
    return *crbegin();
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST pointer data() noexcept {
    return elts;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_pointer
  data() const noexcept {
    return elts;
  }

  // Iterators

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator begin() noexcept {
    return iterator(&elts[0]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
  begin() const noexcept {
    return const_iterator(&elts[0]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
  cbegin() const noexcept {
    return const_iterator(&elts[0]);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST iterator end() noexcept {
    return iterator(&elts[N]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
  end() const noexcept {
    return const_iterator(&elts[N]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_iterator
  cend() const noexcept {
    return const_iterator(&elts[N]);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator
  rbegin() noexcept {
    return reverse_iterator(&elts[N]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
  rbegin() const noexcept {
    return const_reverse_iterator(&elts[N]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
  crbegin() const noexcept {
    return const_reverse_iterator(&elts[N]);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST reverse_iterator
  rend() noexcept {
    return reverse_iterator(&elts[0]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
  rend() const noexcept {
    return const_reverse_iterator(&elts[0]);
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const_reverse_iterator
  crend() const noexcept {
    return const_reverse_iterator(&elts[0]);
  }

  // Capacity

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool empty() const noexcept {
    return N == 0;
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST size_type
  size() const noexcept {
    return N;
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST size_type
  max_size() const noexcept {
    return N;
  }

  // Operations

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void fill(const T &value) {
    for (std::size_t i = 0; i < N; ++i)
      elts[i] = value;
  }

  // We might not have a device-function std::swap
  /*constexpr*/ ARITH_INLINE /*ARITH_DEVICE*/ ARITH_HOST void
  swap(arr &other) noexcept {
    if (this == &other)
      return;
    for (std::size_t i = 0; i < N; ++i)
      std::swap((*this)[i], other[i]);
  }
};

} // namespace Arith

namespace std {

template <std::size_t I, typename T, std::size_t N>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T &
get(Arith::arr<T, N> &arr) noexcept {
  return arr[I];
}
template <std::size_t I, typename T, std::size_t N>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T &&
get(Arith::arr<T, N> &&arr) noexcept {
  return static_cast<T &&>(arr[I]);
}
template <std::size_t I, typename T, std::size_t N>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const T &
get(const Arith::arr<T, N> &arr) noexcept {
  return arr[I];
}
template <std::size_t I, typename T, std::size_t N>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const T &&
get(const Arith::arr<T, N> &&arr) noexcept {
  return static_cast<const T &&>(arr[I]);
}

} // namespace std

namespace Arith {

namespace detail {
template <typename T, class Tuple, std::size_t... Is>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_from_tuple(const Tuple &t, std::index_sequence<Is...>) {
  return arr<T, sizeof...(Is)>{std::get<Is>(t)...};
}
template <typename T, class Tuple, std::size_t... Is>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_from_tuple(Tuple &&t, std::index_sequence<Is...>) {
  return arr<T, sizeof...(Is)>{std::move(std::get<Is>(t))...};
}
} // namespace detail

template <typename T, std::size_t N, typename Tuple>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr<T, N>
array_from_tuple(const Tuple &t) {
  return detail::array_from_tuple<T>(t, std::make_index_sequence<N>());
}
template <typename T, std::size_t N, typename Tuple>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr<T, N>
array_from_tuple(Tuple &&t) {
  return detail::array_from_tuple<T>(std::move(t),
                                     std::make_index_sequence<N>());
}

namespace detail {
template <typename T, class F, std::size_t... Is>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
construct_array(const F &f, std::index_sequence<Is...>) {
  return arr<T, sizeof...(Is)>{f(Is)...};
}
} // namespace detail

template <typename T, std::size_t N, typename F>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr<T, N>
construct_array(const F &f) {
  return detail::construct_array<T>(f, std::make_index_sequence<N>());
}

namespace detail {
template <typename T, class Tuple, std::size_t... Is, typename U>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_push(const Tuple &t, std::index_sequence<Is...>, const U &x) {
  return arr<T, sizeof...(Is) + 1>{std::get<Is>(t)..., T(x)};
}
template <typename T, class Tuple, std::size_t... Is, typename U>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_push(Tuple &&t, std::index_sequence<Is...>, const U &x) {
  return arr<T, sizeof...(Is) + 1>{std::move(std::get<Is>(t))..., T(x)};
}
} // namespace detail

template <class T, std::size_t N, typename U>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_push(const arr<T, N> &a, const U &e) {
  return detail::array_push<T>(a, std::make_index_sequence<N>(), e);
}
template <class T, std::size_t N, typename U>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto array_push(arr<T, N> &&a,
                                                               const U &e) {
  return detail::array_push<T>(std::move(a), std::make_index_sequence<N>(), e);
}

// template <typename T, std::size_t N, typename F>
// constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr<T, N>
// construct_array(const F &f) {
//   if constexpr (N == 0)
//     return arr<T, N>();
//   if constexpr (N > 0)
//     return array_push<T>(construct_array<T, N - 1>(f), f(N - 1));
// }

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_ARR_HXX
