#ifndef CARPETX_ARITH_DIV_HXX
#define CARPETX_ARITH_DIV_HXX

#include "defs.hxx"

namespace Arith {

// Integer division (and related) functions that round towards negative or
// positive infinity

template <typename T> constexpr ARITH_INLINE T div_floor(T x, T y) {
  // C++ division truncates; we want to round towards -infinity instead
  x = y < 0 ? -x : x;
  y = y < 0 ? -y : y;
  return (x >= 0 ? x : x - y + 1) / y;
}

template <typename T> constexpr ARITH_INLINE T align_floor(T x, T y) {
  return div_floor(x, y) * y;
}

template <typename T> constexpr ARITH_INLINE T mod_floor(T x, T y) {
  return x - align_floor(x, y);
}

template <typename T> constexpr ARITH_INLINE T div_ceil(T x, T y) {
  // C++ division truncates; we want to round towards +infinity instead
  x = y < 0 ? -x : x;
  y = y < 0 ? -y : y;
  return (x >= 0 ? x + y - 1 : x) / y;
}

template <typename T> constexpr ARITH_INLINE T align_ceil(T x, T y) {
  return div_ceil(x, y) * y;
}

template <typename T> constexpr ARITH_INLINE T mod_ceil(T x, T y) {
  return x - align_ceil(x, y);
}

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_DIV_HXX
