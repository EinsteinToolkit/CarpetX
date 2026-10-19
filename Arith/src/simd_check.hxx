#ifndef CARPETX_ARITH_SIMD_CHECK_HXX
#define CARPETX_ARITH_SIMD_CHECK_HXX

// This header implements the `SIMD_CHECK` backend for `Arith::simd`/
// `Arith::simdl` (see `simd.hxx`). It is included only when `SIMD_CHECK` is
// defined, e.g. via `-DSIMD_CHECK=8`.
//
// `SIMD_CHECK=<N>` replaces the hardware vector backend (NSIMD) with a
// software-emulated `N`-wide vector ("pack"), independent of the actual
// target CPU's vector width. This lets a bounds-checking build run
// anywhere, without needing e.g. real AVX-512 hardware to test an 8-wide
// (512-bit) `CCTK_REAL` vector.
//
// The point of doing this in software rather than simply asking NSIMD for
// an N-wide pack is the masked load/store family
// (`maskz_loada`/`maskz_loadu`/`masko_loada`/`masko_loadu`/
// `mask_storea`/`mask_storeu`): on real hardware -- and inside NSIMD's
// portable emulation of masked loads/stores on ISAs that lack a native
// masked-load instruction -- a "masked" access can still *physically*
// touch memory for lanes whose mask bit is false (typically: load the
// full vector, then blend/discard the inactive lanes). That is exactly
// how a `loop_int`/`loop_all` tail vector can read or write past the
// interior box or past the allocated grid function array, even though
// the numerical result only ever depends on the active lanes.
//
// `Loop::PointDesc::mask` (computed via `Arith::mask_for_loop_tail`, see
// `loop.hxx`) is `false` for exactly the lanes that lie beyond the current
// loop's upper bound -- the interior boundary for `loop_int`, or the
// allocated-array boundary for `loop_all`. The masked load/store
// implementations below loop over lanes in software and, whenever a lane's
// mask bit is false, never dereference the corresponding memory location
// at all. That makes it impossible for a masked SIMD access driven by
// `loop_int`/`loop_all` to index outside of the region `mask` describes,
// regardless of what the hardware backend would have done.

#include <cctk.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <type_traits>

namespace Arith {
namespace detail {

// Aborts with a CCTK_ERROR (not a raw out-of-bounds memory access) as soon
// as SIMD_CHECK can prove a masked SIMD access touches an "active" lane
// beyond the bound recorded on the mask itself -- i.e. would index outside
// of loop_int's interior, or outside of loop_all's allocated array.
[[noreturn]] inline void checked_bounds_error(const char *kind,
                                              std::size_t width,
                                              std::size_t lane,
                                              std::size_t extent) {
  char buf[512];
  std::snprintf(
      buf, sizeof buf,
      "SIMD_CHECK: %s SIMD access of width %zu touches lane %zu, but the "
      "mask says only %zu lane(s) are valid there (the current loop_int "
      "interior bound or loop_all allocated-array bound, as recorded by "
      "mask_for_loop_tail()). This looks like a mask/pointer that no "
      "longer matches the current loop position -- use a freshly computed "
      "mask (mask_for_loop_tail(), or PointDesc::mask) with the matching "
      "pointer.",
      kind, width, lane, extent);
  CCTK_ERROR(buf);
}

template <std::size_t N> struct checked_packl;

// A software stand-in for `nsimd::pack<T>`.
template <typename T, std::size_t N> struct checked_pack {
  static_assert(N >= 1);
  using value_type = T;
  static constexpr std::size_t size = N;

  std::array<T, N> v;

  constexpr checked_pack() : v{} {}
  constexpr checked_pack(const checked_pack &) = default;
  constexpr checked_pack(checked_pack &&) = default;
  constexpr checked_pack &operator=(const checked_pack &) = default;
  constexpr checked_pack &operator=(checked_pack &&) = default;

  constexpr checked_pack(const T &a) : v{} {
    for (std::size_t n = 0; n < N; ++n)
      v[n] = a;
  }

  constexpr T &operator[](std::size_t n) { return v[n]; }
  constexpr const T &operator[](std::size_t n) const { return v[n]; }
};

// A software stand-in for `nsimd::packl<T>` (a pack of booleans).
//
// `valid_extent` is how many lanes, counted from lane 0, are actually
// valid -- e.g. the interior bound for loop_int, or the allocated-array
// bound for loop_all. `mask_for_loop_tail()` (below) sets it on the mask
// object it returns; the masked load/store functions check active lanes
// against it. Carrying it on the mask value itself (rather than in a
// thread-local/global) keeps this purely functional: no shared mutable
// state, so it stays correct under OpenMP/GPU threading without any
// synchronization, and a mask default-constructed some other way (e.g. by
// `!`, `&&`, or a comparison operator) is simply "unbounded" -- it can
// never spuriously trip a check meant for a different mask.
template <std::size_t N> struct checked_packl {
  static_assert(N >= 1);
  static constexpr std::size_t size = N;

  std::array<bool, N> v;
  std::size_t valid_extent = std::numeric_limits<std::size_t>::max();

  constexpr checked_packl() : v{} {}
  constexpr checked_packl(const checked_packl &) = default;
  constexpr checked_packl(checked_packl &&) = default;
  constexpr checked_packl &operator=(const checked_packl &) = default;
  constexpr checked_packl &operator=(checked_packl &&) = default;

  constexpr checked_packl(bool a) : v{} {
    for (std::size_t n = 0; n < N; ++n)
      v[n] = a;
  }

  constexpr bool &operator[](std::size_t n) { return v[n]; }
  constexpr bool operator[](std::size_t n) const { return v[n]; }
};

////////////////////////////////////////////////////////////////////////////////
// Bitwise operations on floating-point (or integer) lanes, implemented via
// bit-reinterpretation for floating-point types (matching what real SIMD
// bitwise instructions do to floating-point registers).

template <typename T, typename = void> struct bits_of {
  using type = T;
};
template <> struct bits_of<float> {
  using type = std::uint32_t;
};
template <> struct bits_of<double> {
  using type = std::uint64_t;
};

template <typename T> T checked_bit_not(const T &x) {
  if constexpr (std::is_floating_point_v<T>) {
    using U = typename bits_of<T>::type;
    U u;
    std::memcpy(&u, &x, sizeof(T));
    u = ~u;
    T r;
    std::memcpy(&r, &u, sizeof(T));
    return r;
  } else {
    return ~x;
  }
}

#define ARITH_CHECKED_BITOP(name, op)                                        \
  template <typename T> T name(const T &x, const T &y) {                     \
    if constexpr (std::is_floating_point_v<T>) {                             \
      using U = typename bits_of<T>::type;                                   \
      U ux, uy;                                                              \
      std::memcpy(&ux, &x, sizeof(T));                                       \
      std::memcpy(&uy, &y, sizeof(T));                                       \
      U ur = ux op uy;                                                       \
      T r;                                                                   \
      std::memcpy(&r, &ur, sizeof(T));                                       \
      return r;                                                              \
    } else {                                                                 \
      return x op y;                                                         \
    }                                                                        \
  }
ARITH_CHECKED_BITOP(checked_bit_and, &)
ARITH_CHECKED_BITOP(checked_bit_or, |)
ARITH_CHECKED_BITOP(checked_bit_xor, ^)
#undef ARITH_CHECKED_BITOP

////////////////////////////////////////////////////////////////////////////////
// Operators on checked_pack<T,N>

template <typename T, std::size_t N>
inline checked_pack<T, N> operator+(const checked_pack<T, N> &x) {
  return x;
}
template <typename T, std::size_t N>
inline checked_pack<T, N> operator-(const checked_pack<T, N> &x) {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = -x[n];
  return r;
}
template <typename T, std::size_t N>
inline checked_pack<T, N> operator~(const checked_pack<T, N> &x) {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = checked_bit_not(x[n]);
  return r;
}
template <typename T, std::size_t N>
inline checked_packl<N> operator!(const checked_pack<T, N> &x) {
  checked_packl<N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = !x[n];
  return r;
}

#define ARITH_CHECKED_PACK_BINOP(op)                                         \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> operator op(const checked_pack<T, N> &x,         \
                                        const checked_pack<T, N> &y) {       \
    checked_pack<T, N> r;                                                    \
    for (std::size_t n = 0; n < N; ++n)                                      \
      r[n] = x[n] op y[n];                                                   \
    return r;                                                                \
  }                                                                          \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> operator op(const T &a,                          \
                                        const checked_pack<T, N> &y) {       \
    return checked_pack<T, N>(a) op y;                                      \
  }                                                                          \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> operator op(const checked_pack<T, N> &x,         \
                                        const T &b) {                        \
    return x op checked_pack<T, N>(b);                                      \
  }
ARITH_CHECKED_PACK_BINOP(+)
ARITH_CHECKED_PACK_BINOP(-)
ARITH_CHECKED_PACK_BINOP(*)
ARITH_CHECKED_PACK_BINOP(/)
ARITH_CHECKED_PACK_BINOP(%)
#undef ARITH_CHECKED_PACK_BINOP

#define ARITH_CHECKED_PACK_BITBINOP(op, fn)                                  \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> operator op(const checked_pack<T, N> &x,         \
                                        const checked_pack<T, N> &y) {       \
    checked_pack<T, N> r;                                                    \
    for (std::size_t n = 0; n < N; ++n)                                      \
      r[n] = fn(x[n], y[n]);                                                 \
    return r;                                                                \
  }                                                                          \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> operator op(const T &a,                          \
                                        const checked_pack<T, N> &y) {       \
    return checked_pack<T, N>(a) op y;                                      \
  }                                                                          \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> operator op(const checked_pack<T, N> &x,         \
                                        const T &b) {                        \
    return x op checked_pack<T, N>(b);                                      \
  }
ARITH_CHECKED_PACK_BITBINOP(&, checked_bit_and)
ARITH_CHECKED_PACK_BITBINOP(|, checked_bit_or)
ARITH_CHECKED_PACK_BITBINOP(^, checked_bit_xor)
#undef ARITH_CHECKED_PACK_BITBINOP

#define ARITH_CHECKED_PACK_CMP(op)                                           \
  template <typename T, std::size_t N>                                      \
  inline checked_packl<N> operator op(const checked_pack<T, N> &x,           \
                                      const checked_pack<T, N> &y) {          \
    checked_packl<N> r;                                                      \
    for (std::size_t n = 0; n < N; ++n)                                      \
      r[n] = x[n] op y[n];                                                   \
    return r;                                                                \
  }                                                                          \
  template <typename T, std::size_t N>                                      \
  inline checked_packl<N> operator op(const T &a,                            \
                                      const checked_pack<T, N> &y) {          \
    return checked_pack<T, N>(a) op y;                                      \
  }                                                                          \
  template <typename T, std::size_t N>                                      \
  inline checked_packl<N> operator op(const checked_pack<T, N> &x,           \
                                      const T &b) {                          \
    return x op checked_pack<T, N>(b);                                      \
  }
ARITH_CHECKED_PACK_CMP(==)
ARITH_CHECKED_PACK_CMP(!=)
ARITH_CHECKED_PACK_CMP(<)
ARITH_CHECKED_PACK_CMP(>)
ARITH_CHECKED_PACK_CMP(<=)
ARITH_CHECKED_PACK_CMP(>=)
#undef ARITH_CHECKED_PACK_CMP

////////////////////////////////////////////////////////////////////////////////
// Elementwise math functions, used by simd.hxx's generic (non-hardware)
// code paths, which call e.g. `using std::sqrt; return sqrt(x.elts);`.

#define ARITH_CHECKED_PACK_UNARY_FN(name)                                    \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> name(const checked_pack<T, N> &x) {              \
    using std::name;                                                         \
    checked_pack<T, N> r;                                                    \
    for (std::size_t n = 0; n < N; ++n)                                      \
      r[n] = name(x[n]);                                                     \
    return r;                                                                \
  }
ARITH_CHECKED_PACK_UNARY_FN(abs)
ARITH_CHECKED_PACK_UNARY_FN(sqrt)
ARITH_CHECKED_PACK_UNARY_FN(cbrt)
ARITH_CHECKED_PACK_UNARY_FN(exp)
ARITH_CHECKED_PACK_UNARY_FN(exp2)
ARITH_CHECKED_PACK_UNARY_FN(log)
ARITH_CHECKED_PACK_UNARY_FN(log2)
ARITH_CHECKED_PACK_UNARY_FN(log10)
ARITH_CHECKED_PACK_UNARY_FN(sin)
ARITH_CHECKED_PACK_UNARY_FN(cos)
ARITH_CHECKED_PACK_UNARY_FN(tan)
ARITH_CHECKED_PACK_UNARY_FN(asin)
ARITH_CHECKED_PACK_UNARY_FN(acos)
ARITH_CHECKED_PACK_UNARY_FN(atan)
ARITH_CHECKED_PACK_UNARY_FN(sinh)
ARITH_CHECKED_PACK_UNARY_FN(cosh)
ARITH_CHECKED_PACK_UNARY_FN(tanh)
ARITH_CHECKED_PACK_UNARY_FN(asinh)
ARITH_CHECKED_PACK_UNARY_FN(acosh)
ARITH_CHECKED_PACK_UNARY_FN(atanh)
#undef ARITH_CHECKED_PACK_UNARY_FN

// `signbit` (unlike the functions above) returns a logical pack, not a
// value pack.
template <typename T, std::size_t N>
inline checked_packl<N> signbit(const checked_pack<T, N> &x) {
  using std::signbit;
  checked_packl<N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = signbit(x[n]);
  return r;
}

// `fabs` is a synonym for `abs` on floating-point types, but is looked up
// separately in simd.hxx (`using std::abs;`/`using std::fabs;` at
// different call sites).
template <typename T, std::size_t N>
inline checked_pack<T, N> fabs(const checked_pack<T, N> &x) {
  using std::fabs;
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = fabs(x[n]);
  return r;
}

#define ARITH_CHECKED_PACK_BINARY_FN(name)                                   \
  template <typename T, std::size_t N>                                      \
  inline checked_pack<T, N> name(const checked_pack<T, N> &x,                \
                                 const checked_pack<T, N> &y) {               \
    using std::name;                                                         \
    checked_pack<T, N> r;                                                    \
    for (std::size_t n = 0; n < N; ++n)                                      \
      r[n] = name(x[n], y[n]);                                               \
    return r;                                                                \
  }
ARITH_CHECKED_PACK_BINARY_FN(pow)
ARITH_CHECKED_PACK_BINARY_FN(hypot)
ARITH_CHECKED_PACK_BINARY_FN(copysign)
ARITH_CHECKED_PACK_BINARY_FN(max)
ARITH_CHECKED_PACK_BINARY_FN(min)
#undef ARITH_CHECKED_PACK_BINARY_FN

template <typename T, std::size_t N>
inline checked_packl<N> to_logical(const checked_pack<T, N> &x) {
  checked_packl<N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = bool(x[n]);
  return r;
}

template <typename T, std::size_t N>
inline checked_pack<T, N> fma(const checked_pack<T, N> &x,
                              const checked_pack<T, N> &y,
                              const checked_pack<T, N> &z) {
  using std::fma;
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = fma(x[n], y[n], z[n]);
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// Operators and reductions on checked_packl<N>

template <std::size_t N>
inline checked_packl<N> operator!(const checked_packl<N> &x) {
  checked_packl<N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = !x[n];
  return r;
}

#define ARITH_CHECKED_PACKL_BINOP(op)                                        \
  template <std::size_t N>                                                   \
  inline checked_packl<N> operator op(const checked_packl<N> &x,             \
                                      const checked_packl<N> &y) {            \
    checked_packl<N> r;                                                      \
    for (std::size_t n = 0; n < N; ++n)                                      \
      r[n] = x[n] op y[n];                                                   \
    return r;                                                                \
  }                                                                          \
  template <std::size_t N>                                                   \
  inline checked_packl<N> operator op(bool a, const checked_packl<N> &y) {    \
    return checked_packl<N>(a) op y;                                        \
  }                                                                          \
  template <std::size_t N>                                                   \
  inline checked_packl<N> operator op(const checked_packl<N> &x, bool b) {    \
    return x op checked_packl<N>(b);                                        \
  }
ARITH_CHECKED_PACKL_BINOP(&&)
ARITH_CHECKED_PACKL_BINOP(||)
ARITH_CHECKED_PACKL_BINOP(&)
ARITH_CHECKED_PACKL_BINOP(|)
ARITH_CHECKED_PACKL_BINOP(^)
#undef ARITH_CHECKED_PACKL_BINOP

template <std::size_t N> inline bool all(const checked_packl<N> &x) {
  bool r = true;
  for (std::size_t n = 0; n < N; ++n)
    r = r && x[n];
  return r;
}
template <std::size_t N> inline bool any(const checked_packl<N> &x) {
  bool r = false;
  for (std::size_t n = 0; n < N; ++n)
    r = r || x[n];
  return r;
}

template <std::size_t N>
inline checked_packl<N> andnotl(const checked_packl<N> &x,
                                const checked_packl<N> &y) {
  checked_packl<N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = x[n] && !y[n];
  return r;
}

template <typename T, std::size_t N>
inline checked_pack<T, N> andnotb(const checked_pack<T, N> &x,
                                  const checked_pack<T, N> &y) {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = checked_bit_and(x[n], checked_bit_not(y[n]));
  return r;
}

template <std::size_t N>
inline checked_packl<N> if_else1(const checked_packl<N> &cond,
                                 const checked_packl<N> &x,
                                 const checked_packl<N> &y) {
  checked_packl<N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = cond[n] ? x[n] : y[n];
  return r;
}
template <typename T, std::size_t N>
inline checked_pack<T, N> if_else1(const checked_packl<N> &cond,
                                   const checked_pack<T, N> &x,
                                   const checked_pack<T, N> &y) {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = cond[n] ? x[n] : y[n];
  return r;
}

////////////////////////////////////////////////////////////////////////////////
// The core of `SIMD_CHECK`: bounds-safe loads and stores.
//
// A masked operation below only ever reads or writes `ptr[n]` for lanes `n`
// where `mask[n]` is true. In particular, if `mask` was computed by
// `checked_mask_for_loop_tail(i, n)` (below) for the interior/all box
// bound `n`, then a masked load/store at pointer `ptr` (pointing at grid
// index `i` along the vectorized direction) can never touch `ptr[k]` for
// `i+k >= n`, i.e. never indexes past that bound. This holds unconditionally
// -- it needs no bound check at all, active lanes are simply always safe.
//
// Unmasked load/store (`loadu`/`storea`) have no mask to carry a bound on,
// so unlike the masked family below they cannot be checked here; they trust
// the caller for all N lanes, exactly like a real hardware unmasked SIMD
// load/store would.
//
// `mask.valid_extent` only proves a mask is internally self-consistent: a
// mask computed correctly for one position, then reused unchanged at a
// different (pointer, loop-index) position, still passes that check. Where
// the true bounds of what is being accessed are independently known --
// `GF3D2`/`GF3D5` in `loop.hxx` know their own array's real extent, from
// their `layout`, at every access -- callers should use this instead: it
// takes that ground truth (`i0`, the grid x-index this access starts at,
// and `imax0`, the layout's real upper bound) fresh from the caller on
// every call, so a stale mask/pointer is caught even though the mask
// itself is not stale. This needs no shared state either: `i0`/`imax0` are
// ordinary values passed in by the caller, not anything ambient.
template <std::size_t N>
inline void checked_verify_mask_bounds(const checked_packl<N> &mask, int i0,
                                       int imax0, const char *kind) {
  const std::size_t extent =
      imax0 <= i0 ? std::size_t(0) : std::size_t(imax0 - i0);
  for (std::size_t n = 0; n < N; ++n)
    if (mask[n] && n >= extent)
      checked_bounds_error(kind, N, n, extent);
}

template <typename T, std::size_t N>
inline checked_pack<T, N> checked_loada(const T *ptr) {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = ptr[n];
  return r;
}

template <typename T, std::size_t N>
inline void checked_storea(T *ptr, const checked_pack<T, N> &x) {
  for (std::size_t n = 0; n < N; ++n)
    ptr[n] = x[n];
}

// Masked load/store: never dereferences ptr[n] for an inactive lane (so
// this is unconditionally safe), but still errors out if an *active* lane
// lies beyond `mask.valid_extent` -- which means `mask` is stale, or was
// computed for a different pointer/loop position than the one it is being
// used with here.
template <typename T, std::size_t N>
inline checked_pack<T, N> checked_maskz_loadu(const checked_packl<N> &mask,
                                              const T *ptr) {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n) {
    if (!mask[n]) {
      r[n] = T(0);
      continue;
    }
    if (n >= mask.valid_extent)
      checked_bounds_error("masked load", N, n, mask.valid_extent);
    r[n] = ptr[n];
  }
  return r;
}

template <typename T, std::size_t N>
inline checked_pack<T, N> checked_masko_loadu(const checked_packl<N> &mask,
                                              const T *ptr,
                                              const checked_pack<T, N> &other) {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n) {
    if (!mask[n]) {
      r[n] = other[n];
      continue;
    }
    if (n >= mask.valid_extent)
      checked_bounds_error("masked load", N, n, mask.valid_extent);
    r[n] = ptr[n];
  }
  return r;
}

template <typename T, std::size_t N>
inline void checked_mask_storeu(const checked_packl<N> &mask, T *ptr,
                                const checked_pack<T, N> &x) {
  for (std::size_t n = 0; n < N; ++n) {
    if (!mask[n])
      continue;
    if (n >= mask.valid_extent)
      checked_bounds_error("masked store", N, n, mask.valid_extent);
    ptr[n] = x[n];
  }
}

template <typename T, std::size_t N>
inline checked_pack<T, N> checked_iota() {
  checked_pack<T, N> r;
  for (std::size_t n = 0; n < N; ++n)
    r[n] = T(n);
  return r;
}

// `checked_mask_for_loop_tail(i, n)[lane]` is true iff `i+lane < n`, i.e.
// iff the lane's grid index along the vectorized direction is still below
// the loop's upper bound `n` (the interior bound for `loop_int`, or the
// "all"/allocated-array bound for `loop_all`; see `PointDesc::mask` and
// `GridDescBase::box_int`/`box_all` in `loop.hxx`). The returned mask
// carries that bound (`max(0, n-i)`) as `valid_extent`, for the masked
// load/store functions above to check active lanes against -- purely a
// property of this one returned value, not any shared state.
template <std::size_t N>
inline checked_packl<N> checked_mask_for_loop_tail(int i, int n) {
  checked_packl<N> r;
  r.valid_extent = n <= i ? std::size_t(0) : std::size_t(n - i);
  for (std::size_t lane = 0; lane < N; ++lane)
    r[lane] = (i + int(lane)) < n;
  return r;
}

} // namespace detail
} // namespace Arith

#endif // #ifndef CARPETX_ARITH_SIMD_CHECK_HXX
