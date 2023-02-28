#include "prolongate_3d_rf2.hxx"
#include "timer.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <AMReX_Gpu.H>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
#endif

#include <array>
#include <cassert>
#include <cmath>
#include <mutex>
#include <vector>

namespace CarpetX {
using namespace std;

std::ostream &operator<<(std::ostream &os, const centering_t cent) {
  switch (cent) {
  case centering_t::vertex:
    return os << "vertex";
  case centering_t::cell:
    return os << "cell";
  }
  assert(false);
}

std::ostream &operator<<(std::ostream &os, const interpolation_t intp) {
  switch (intp) {
  case interpolation_t::poly:
    return os << "poly";
  case interpolation_t::cons:
    return os << "cons";
  case interpolation_t::hermite:
    return os << "hermite";
  }
  assert(false);
}

// 1D interpolation coefficients

template <centering_t CENT, interpolation_t INTP, int ORDER, typename T>
struct coeffs1d;

template <typename T> struct coeffs1d<VC, POLY, /*order*/ 1, T> {
  static constexpr array<T, 2> coeffs = {
      +1 / T(2),
      +1 / T(2),
  };
};
template <typename T> struct coeffs1d<VC, POLY, /*order*/ 3, T> {
  static constexpr array<T, 4> coeffs = {
      -1 / T(16),
      +9 / T(16),
      +9 / T(16),
      -1 / T(16),
  };
};
template <typename T> struct coeffs1d<VC, POLY, /*order*/ 5, T> {
  static constexpr array<T, 6> coeffs = {
      +3 / T(256),  -25 / T(256), +75 / T(128),
      +75 / T(128), -25 / T(256), +3 / T(256),
  };
};
template <typename T> struct coeffs1d<VC, POLY, /*order*/ 7, T> {
  static constexpr array<T, 8> coeffs = {
      -5 / T(2048),    +49 / T(2048),  -245 / T(2048), +1225 / T(2048),
      +1225 / T(2048), -245 / T(2048), +49 / T(2048),  -5 / T(2048),
  };
};

template <typename T> struct coeffs1d<CC, POLY, /*order*/ 0, T> {
  static constexpr array<T, 1> coeffs = {
      +1 / T(1),
  };
};
template <typename T> struct coeffs1d<CC, POLY, /*order*/ 1, T> {
  static constexpr array<T, 2> coeffs = {
      +1 / T(4),
      +3 / T(4),
  };
};
template <typename T> struct coeffs1d<CC, POLY, /*order*/ 2, T> {
  static constexpr array<T, 3> coeffs = {
      +5 / T(32),
      +15 / T(16),
      -3 / T(32),
  };
};
template <typename T> struct coeffs1d<CC, POLY, /*order*/ 3, T> {
  static constexpr array<T, 4> coeffs = {
      -5 / T(128),
      +35 / T(128),
      +105 / T(128),
      -7 / T(128),
  };
};
template <typename T> struct coeffs1d<CC, POLY, /*order*/ 4, T> {
  static constexpr array<T, 5> coeffs = {
      -45 / T(2048), +105 / T(512), +945 / T(1024), -63 / T(512), +35 / T(2048),
  };
};

template <typename T> struct coeffs1d<VC, CONS, /*order*/ 0, T> {
  static constexpr array<T, 1> coeffs0 = {
      1 / T(1),
  };
  static constexpr array<T, 0> coeffs1 = {};
};
template <typename T> struct coeffs1d<VC, CONS, /*order*/ 1, T> {
  static constexpr array<T, 1> coeffs0 = {
      1 / T(1),
  };
  static constexpr array<T, 2> coeffs1 = {
      +1 / T(2),
      +1 / T(2),
  };
};
// template <typename T> struct coeffs1d<VC, CONS, /*order*/ 2, T> {
//   static constexpr array<T, 3> coeffs0 = {
//       -1 / T(32),
//       +17 / T(16),
//       -1 / T(32),
//   };
//   static constexpr array<T, 3> coeffs1 = {
//       +13 / T(16),
//       -5 / T(32),
//       +11 / T(32),
//   };
// };
// template <typename T> struct coeffs1d<VC, CONS, /*order*/ 4, T> {
//   static constexpr array<T, 5> coeffs0 = {
//       +7 / T(2048), -23 / T(512), +1109 / T(1024), -23 / T(512), +7 /
//       T(2048),
//   };
//   static constexpr array<T, 5> coeffs1 = {
//       +63 / T(2048), -103 / T(512), +781 / T(1024),
//       +233 / T(512), -97 / T(2048),
//   };
// };

template <typename T> struct coeffs1d<CC, CONS, /*order*/ 0, T> {
  static constexpr array<T, 1> coeffs = {
      +1 / T(1),
  };
};
template <typename T> struct coeffs1d<CC, CONS, /*order*/ 2, T> {
  static constexpr array<T, 3> coeffs = {
      +1 / T(8),
      +1 / T(1),
      -1 / T(8),
  };
};
template <typename T> struct coeffs1d<CC, CONS, /*order*/ 4, T> {
  static constexpr array<T, 5> coeffs = {
      -3 / T(128), +11 / T(64), +1 / T(1), -11 / T(64), +3 / T(128),
  };
};

// 1D interpolation operators

template <centering_t CENT, interpolation_t INTP, int ORDER> struct interp1d;
// static constexpr int required_ghosts;
// template <typename T>
// inline T operator()(const T *restrict const crseptr, const ptrdiff_t di,
//                     const int off) const;

// off=0: on coarse point
// off=1: between coarse points
template <int ORDER> struct interp1d<VC, POLY, ORDER> {
  static_assert(ORDER % 2 == 1);
  static constexpr int required_ghosts = (ORDER + 1) / 2;
  template <typename T>
  CCTK_DEVICE CCTK_HOST inline T operator()(const T *restrict const crseptr,
                                            const ptrdiff_t di,
                                            const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    if (off == 0)
      return crseptr[0];
    constexpr array<T, ORDER + 1> cs = coeffs1d<VC, POLY, ORDER, T>::coeffs;
    const int i0 = (ORDER + 1) / 2 - off;
    constexpr int i0min = (ORDER + 1) / 2 - 1;
    constexpr int i0max = (ORDER + 1) / 2;
    constexpr int imin = 0;
    constexpr int imax = (ORDER + 1) / 2 - 1;
    constexpr int i1min = ORDER - imax;
    constexpr int i1max = ORDER - imin;
    // nvcc doesn't accept the constexpr terms below
#ifndef __CUDACC__
    const auto abs0 = [](auto x) { return x >= 0 ? x : -x; };
    static_assert(abs0(imin - i0min) <= required_ghosts, "");
    static_assert(abs0(imin - i0max) <= required_ghosts, "");
    static_assert(abs0(imax - i0min) <= required_ghosts, "");
    static_assert(abs0(imax - i0max) <= required_ghosts, "");
    static_assert(abs0(i1min - i0min) <= required_ghosts, "");
    static_assert(abs0(i1min - i0max) <= required_ghosts, "");
    static_assert(abs0(i1max - i0min) <= required_ghosts, "");
    static_assert(abs0(i1max - i0max) <= required_ghosts, "");
#endif
    T y = 0;
    // Make use of symmetry in coefficients
    for (int i = 0; i < (ORDER + 1) / 2; ++i) {
      const int i1 = ORDER - i;
#ifdef CCTK_DEBUG
      assert(cs[i1] == cs[i]);
#endif
      y += cs[i] * (crseptr[(i - i0) * di] + crseptr[(i1 - i0) * di]);
    }
#ifdef CCTK_DEBUG
    assert(isfinite(y));
#endif
#ifdef CCTK_DEBUG
    T y1 = 0;
    for (int i = 0; i < ORDER + 1; ++i)
      y1 += cs[i] * crseptr[(i - i0) * di];
    assert(isfinite(y1));
    // Don't check for equality; there can be round-off errors
    // assert(y1 == y);
#endif
    return y;
  }
};

// off=0: left sub-cell
// off=1: right sub-cell
template <int ORDER> struct interp1d<CC, POLY, ORDER> {
  static constexpr int required_ghosts = (ORDER + 1) / 2;
  template <typename T>
  CCTK_DEVICE CCTK_HOST inline T operator()(const T *restrict const crseptr,
                                            const ptrdiff_t di,
                                            const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    constexpr array<T, ORDER + 1> cs = coeffs1d<CC, POLY, ORDER, T>::coeffs;
    constexpr int i0 = (ORDER + 1) / 2;
    constexpr int imin = 0;
    constexpr int imax = ORDER;
    // nvcc doesn't accept the constexpr terms below
#ifndef __CUDACC__
    const auto abs0 = [](auto x) { return x >= 0 ? x : -x; };
    static_assert(abs0(imin - i0) <= required_ghosts, "");
    static_assert(abs0(imax - i0) <= required_ghosts, "");
#endif
    T y = 0;
    if (off == 0)
      for (int i = 0; i < ORDER + 1; ++i)
        y += cs[i] * crseptr[(i - i0) * di];
    else
      for (int i = 0; i < ORDER + 1; ++i)
        // For odd orders, the stencil has an even number of points and is thus
        // offset. This offset moves the stencil right by one point when it is
        // reversed.
        y += cs[i] * crseptr[(ORDER + (ORDER % 2 != 0 ? 1 : 0) - i - i0) * di];
#ifdef CCTK_DEBUG
    assert(isfinite(y));
#endif
    return y;
  }
};

// off=0: on coarse point
// off=1: between coarse points
template <int ORDER> struct interp1d<VC, CONS, ORDER> {
  static constexpr int required_ghosts = 0; // TODO: fix this
  template <typename T>
  CCTK_DEVICE CCTK_HOST inline T operator()(const T *restrict const crseptr,
                                            const ptrdiff_t di,
                                            const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    T y = 0;
    // TODO: use symmetry
    if (off == 0) {
      constexpr int i0 = ORDER / 2;
      constexpr array<T, ORDER / 2 * 2 + 1> cs =
          coeffs1d<VC, CONS, ORDER, T>::coeffs0;
      for (int i = 0; i < ORDER / 2 * 2 + 1; ++i)
        y += cs[i] * crseptr[(i - i0) * di];
    } else {
      constexpr int i0 = (ORDER + 1) / 2;
      constexpr array<T, (ORDER + 1) / 2 * 2> cs =
          coeffs1d<VC, CONS, ORDER, T>::coeffs1;
      for (int i = 0; i < (ORDER + 1) / 2 * 2; ++i)
        y += cs[i] * crseptr[(i - i0) * di];
    }
    return y;
  }
};
// template <int ORDER> struct interp1d<VC, CONS, ORDER> {
//   static_assert(ORDER % 2 == 0);
//   static_assert(ORDER > 0);
//   template <typename T>
//   inline T operator()(const T *restrict const crseptr, const ptrdiff_t di,
//                       const int off) const {
// #ifdef CCTK_DEBUG
//     assert(off == 0 || off == 1);
// #endif
//     constexpr int i0 = (ORDER + 1) / 2;
//     T y = 0;
//     if (off == 0) {
//       constexpr array<T, ORDER + 1> cs = coeffs1d<VC, CONS, ORDER,
//       T>::coeffs0; for (int i = 0; i < ORDER + 1; ++i)
//         y += cs[i] * crseptr[(i - i0) * di];
//     } else {
//       constexpr array<T, ORDER + 1> cs = coeffs1d<VC, CONS, ORDER,
//       T>::coeffs1; for (int i = 0; i < ORDER + 1; ++i)
//         y += cs[i] * crseptr[(i - i0) * di];
//     }
//     return y;
//   }
// };

// off=0: left sub-cell
// off=1: right sub-cell
template <int ORDER> struct interp1d<CC, CONS, ORDER> {
  static_assert(ORDER % 2 == 0, "");
  static constexpr int required_ghosts = (ORDER + 1) / 2;
  template <typename T>
  CCTK_DEVICE CCTK_HOST inline T operator()(const T *restrict const crseptr,
                                            const ptrdiff_t di,
                                            const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    constexpr array<T, ORDER + 1> cs = coeffs1d<CC, CONS, ORDER, T>::coeffs;
    constexpr int i0 = (ORDER + 1) / 2;
    constexpr int imin = 0;
    constexpr int imax = ORDER;
    // nvcc doesn't accept the constexpr terms below
#ifndef __CUDACC__
    const auto abs0 = [](auto x) { return x >= 0 ? x : -x; };
    static_assert(abs0(imin - i0) <= required_ghosts, "");
    static_assert(abs0(imax - i0) <= required_ghosts, "");
#endif
    T y = 0;
    // if (off == 0)
    //   for (int i = 0; i < ORDER + 1; ++i)
    //     y += cs[i] * crseptr[(i - i0) * di];
    // else
    //   for (int i = 0; i < ORDER + 1; ++i)
    //     y += cs[i] * crseptr[-(i - i0) * di];
    if (off == 0) {
      // Make use of symmetry in coefficients
      for (int i = 0; i < ORDER / 2; ++i) {
        const int i1 = ORDER - i;
#ifdef CCTK_DEBUG
        assert(cs[i1] == -cs[i]);
#endif
        y += cs[i] * (crseptr[(i - i0) * di] - crseptr[(i1 - i0) * di]);
      }
      y += cs[ORDER / 2] * crseptr[(ORDER / 2 - i0) * di];
    } else {
      for (int i = 0; i < ORDER / 2; ++i) {
        const int i1 = ORDER - i;
#ifdef CCTK_DEBUG
        assert(cs[i1] == -cs[i]);
#endif
        y += cs[i] * (crseptr[(i1 - i0) * di] - crseptr[(i - i0) * di]);
      }
      y += cs[ORDER / 2] * crseptr[(ORDER / 2 - i0) * di];
    }
    return y;
  }
};

// Test 1d interpolators

template <centering_t CENT, interpolation_t INTP, int ORDER, typename T>
struct test_interp1d;

template <centering_t CENT, int ORDER, typename T>
struct test_interp1d<CENT, POLY, ORDER, T> {
  test_interp1d() {
    for (int order = 0; order <= ORDER; ++order) {
      auto f = [&](T x) { return pown(x, order); };
      constexpr int n = (ORDER + 1) / 2 * 2 + 1;
      array<T, n + 2> ys;
      ys[0] = ys[n + 1] = 0 / T(0);
      constexpr int i0 = n / 2;
      static_assert(interp1d<CENT, POLY, ORDER>::required_ghosts <= i0, "");
      static_assert(interp1d<CENT, POLY, ORDER>::required_ghosts <= n - i0, "");
      for (int i = 0; i < n; ++i) {
        T x = (i - i0) + int(CENT) / T(2);
        T y = f(x);
        ys[i + 1] = y;
      }
      for (int off = 0; off < 2; ++off) {
        T x = int(CENT) / T(4) + off / T(2);
        T y = f(x);
        T y1 = interp1d<CENT, POLY, ORDER>()(&ys[i0 + 1], 1, off);
        // We carefully choose the test problem so that round-off
        // cannot be a problem here
        assert(isfinite(y1));
        assert(y1 == y);
      }
    }
  }
};

template <centering_t CENT, int ORDER, typename T>
struct test_interp1d<CENT, CONS, ORDER, T> {
  test_interp1d() {
    for (int order = 0; order <= ORDER; ++order) {
      // Function f, a polynomial
      // const auto f{[&](T x) { return (order + 1) * pown(x, order); }};
      // Integral of f (antiderivative)
      const auto fint{[&](T x) { return pown(x, order + 1); }};
      constexpr int n = (ORDER + 1) / 2 * 2 + 1;
      if (CENT == CC) {
        array<T, n + 2> xs;
        array<T, n + 2> ys;
        xs[0] = xs[n + 1] = 0 / T(0);
        ys[0] = ys[n + 1] = 0 / T(0);
        constexpr int i0 = n / 2;
        static_assert(interp1d<CENT, CONS, ORDER>::required_ghosts <= i0, "");
        static_assert(interp1d<CENT, CONS, ORDER>::required_ghosts <= n - i0,
                      "");
        for (int i = 0; i < n; ++i) {
          const T x = (i - i0) + int(CENT) / T(2);
          // T y = f(x);
          const T dx = 1;
          const T xlo = x - dx / 2;
          const T xhi = x + dx / 2;
          const T y = fint(xhi) - fint(xlo); // average of f over cell
          xs[i + 1] = x;
          ys[i + 1] = y;
        }
        array<T, 2> x1;
        array<T, 2> y1;
        for (int off = 0; off < 2; ++off) {
          x1[off] = int(CENT) / T(4) + off / T(2);
          y1[off] = interp1d<CENT, CONS, ORDER>()(&ys[i0 + 1], 1, off);
          assert(isfinite(y1[off]));
        }
        assert(y1[0] / 2 + y1[1] / 2 == ys[i0 + 1]);
        const T dx = x1[1] - x1[0];
        const T xlo = x1[0] - dx / 2;
        const T xhi = x1[1] + dx / 2;
        const T yint = fint(xhi) - fint(xlo);
        assert(y1[0] * dx + y1[1] * dx == yint);
      } else {
        // Don't test this, the case (VC,CONS) should not be used
        if (false) {
          array<T, n + 3> xs, ys;
          xs[0] = xs[n + 2] = 0 / T(0);
          ys[0] = ys[n + 2] = 0 / T(0);
          constexpr int i0 = n / 2;
          // TODO
          // static_assert(
          //     interp1d<CENT, CONS, ORDER>::required_ghosts <= i0,
          //     "");
          // static_assert(interp1d<CENT, CONS, ORDER>::required_ghosts
          // <=
          //                   n - i0,
          //               "");
          for (int i = -1; i < n; ++i) {
            const T x = (i - i0) + int(CENT) / T(2);
            // T y = f(x);
            const T dx = 1;
            const T xlo = x - dx / 2;
            const T xhi = x + dx / 2;
            const T y = fint(xhi) - fint(xlo);
            xs[i + 2] = x;
            ys[i + 2] = y;
          }
          array<T, 3> x1;
          array<T, 3> y1;
          for (int off = -1; off < 2; ++off) {
            x1[off + 1] = int(CENT) / T(4) + off / T(2);
            if (off < 0)
              y1[off + 1] =
                  interp1d<CENT, CONS, ORDER>()(&ys[i0 + 1], 1, off + 2);
            else
              y1[off + 1] = interp1d<CENT, CONS, ORDER>()(&ys[i0 + 2], 1, off);
            assert(isfinite(y1[off + 1]));
          }
          const T dx = x1[1] - x1[0];
          const T xlo = x1[0];
          const T xhi = x1[2];
          const T yint = fint(xhi) - fint(xlo);
          if (!(y1[0] / 4 + y1[1] / 2 + y1[2] / 4 == ys[i0 + 2]) ||
              !(y1[0] * dx / 2 + y1[1] * dx + y1[2] * dx / 2 == yint)) {
            cerr << "settings: CENT=" << CENT << " ORDER=" << ORDER
                 << " order=" << order << "\n";
            cerr << "input:\n";
            for (int i = -1; i < n; ++i)
              cerr << "  xs=" << xs[i + 2] << " ys=" << ys[i + 2] << "\n";
            cerr << "output:\n";
            for (int off = -1; off < 2; ++off)
              cerr << "  x1=" << x1[off + 1] << " y1=" << y1[off + 1] << "\n";
            cerr << "xlo=" << xlo << " xhi=" << xhi << " yint=" << yint << "\n";
          }
          assert(y1[0] / 4 + y1[1] / 2 + y1[2] / 4 == ys[i0 + 2]);
          assert(y1[0] * dx / 2 + y1[1] * dx + y1[2] * dx / 2 == yint);
        }
      }
    }
  }
};

////////////////////////////////////////////////////////////////////////////////

template <centering_t CENT, interpolation_t INTP, int ORDER, int D, typename T>
void interp3d(const T *restrict const crseptr,
              const amrex::Box &restrict crsebox, T *restrict const fineptr,
              const amrex::Box &restrict finebox,
              const amrex::Box &restrict targetbox) {
  static test_interp1d<CENT, INTP, ORDER, T> test;

  static_assert(D >= 0 && D < 3, "");

  assert(crseptr);
  assert(crsebox.ok());
  assert(fineptr);
  assert(finebox.ok());
  assert(targetbox.ok());

  const amrex::IntVect first_crseind(finebox.loVect());
  amrex::IntVect next_crseind = first_crseind;
  ++next_crseind.getVect()[D];
  const ptrdiff_t di =
      crsebox.index(next_crseind) - crsebox.index(first_crseind);
  assert(di > 0);

  constexpr int required_ghosts = interp1d<CENT, INTP, ORDER>::required_ghosts;
  {
    const amrex::IntVect fineind(targetbox.loVect());
    amrex::IntVect crseind = fineind;
    crseind.getVect()[D] =
        amrex::coarsen(fineind.getVect()[D], 2) - required_ghosts;
    for (int d = 0; d < 3; ++d)
      assert(crseind.getVect()[d] >= crsebox.loVect()[d]);
    for (int d = 0; d < 3; ++d)
      assert(targetbox.loVect()[d] >= finebox.loVect()[d]);
  }
  {
    const amrex::IntVect fineind(targetbox.hiVect());
    amrex::IntVect crseind = fineind;
    crseind.getVect()[D] =
        amrex::coarsen(fineind.getVect()[D], 2) + required_ghosts;
    for (int d = 0; d < 3; ++d)
      assert(crseind.getVect()[d] <= crsebox.hiVect()[d]);
    for (int d = 0; d < 3; ++d)
      assert(targetbox.hiVect()[d] <= finebox.hiVect()[d]);
  }

  const ptrdiff_t fined0 = finebox.index(amrex::IntVect(0, 0, 0));
  constexpr ptrdiff_t finedi = 1;
  assert(finebox.index(amrex::IntVect(1, 0, 0)) - fined0 == finedi);
  const ptrdiff_t finedj = finebox.index(amrex::IntVect(0, 1, 0)) - fined0;
  const ptrdiff_t finedk = finebox.index(amrex::IntVect(0, 0, 1)) - fined0;

  const ptrdiff_t crsed0 = crsebox.index(amrex::IntVect(0, 0, 0));
  constexpr ptrdiff_t crsedi = 1;
  assert(crsebox.index(amrex::IntVect(1, 0, 0)) - crsed0 == crsedi);
  const ptrdiff_t crsedj = crsebox.index(amrex::IntVect(0, 1, 0)) - crsed0;
  const ptrdiff_t crsedk = crsebox.index(amrex::IntVect(0, 0, 1)) - crsed0;

  const auto kernel = [=] CCTK_DEVICE(
                          const int i, const int j,
                          const int k) CCTK_ATTRIBUTE_ALWAYS_INLINE {
#if 0
        const amrex::IntVect fineind(i, j, k);
        amrex::IntVect crseind = fineind;
        crseind.getVect()[D] = coarsen(fineind.getVect()[D], 2);
        const int off = fineind.getVect()[D] - crseind.getVect()[D] * 2;
        fineptr[finebox.index(fineind)] =
            interp1d<CENT, INTP, ORDER>()(
                &crseptr[crsebox.index(crseind)], di, off);
#ifdef CCTK_DEBUG
        assert(isfinite(fineptr[finebox.index(fineind)]));
#endif
#endif
    // Note: fineind = 2 * coarseind + off
    const int ci = D == 0 ? i >> 1 : i;
    const int cj = D == 1 ? j >> 1 : j;
    const int ck = D == 2 ? k >> 1 : k;
    const int off = (D == 0 ? i : D == 1 ? j : k) & 0x1;
    if (D == 0) {
      // allow vectorization
      const T *restrict const ptr =
          &crseptr[crsed0 + ck * crsedk + cj * crsedj + ci * crsedi];
      const T res0 = interp1d<CENT, INTP, ORDER>()(ptr, di, 0);
      const T res1 = interp1d<CENT, INTP, ORDER>()(ptr, di, 1);
      const T res = off == 0 ? res0 : res1;
      fineptr[fined0 + k * finedk + j * finedj + i * finedi] = res;
    } else {
      fineptr[fined0 + k * finedk + j * finedj + i * finedi] =
          interp1d<CENT, INTP, ORDER>()(
              &crseptr[crsed0 + ck * crsedk + cj * crsedj + ci * crsedi], di,
              off);
    }
#ifdef CCTK_DEBUG
    assert(isfinite(fineptr[fined0 + i * finedi + j * finedj + k * finedk]));
#endif
  };

  // TODO: Use `loop_region` from driver.cxx (to be moved to loop.hxx)
#ifndef __CUDACC__
  // CPU

  for (int k = targetbox.smallEnd()[2]; k <= targetbox.bigEnd()[2]; ++k)
    for (int j = targetbox.smallEnd()[1]; j <= targetbox.bigEnd()[1]; ++j)
      for (int i = targetbox.smallEnd()[0]; i <= targetbox.bigEnd()[0]; ++i)
        kernel(i, j, k);

#else // GPU

  amrex::launch(targetbox, [=] CCTK_DEVICE(const amrex::Box &box)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
#ifdef CCTK_DEBUG
                                 assert(box.bigEnd()[0] == box.smallEnd()[0] &&
                                        box.bigEnd()[1] == box.smallEnd()[1] &&
                                        box.bigEnd()[2] == box.smallEnd()[2]);
#endif
                                 const int i = box.smallEnd()[0];
                                 const int j = box.smallEnd()[1];
                                 const int k = box.smallEnd()[2];

                                 kernel(i, j, k);
                               });

#endif // GPU
}

////////////////////////////////////////////////////////////////////////////////

template <centering_t CENTI, centering_t CENTJ, centering_t CENTK,
          interpolation_t INTPI, interpolation_t INTPJ, interpolation_t INTPK,
          int ORDERI, int ORDERJ, int ORDERK>
prolongate_3d_rf2<CENTI, CENTJ, CENTK, INTPI, INTPJ, INTPK, ORDERI, ORDERJ,
                  ORDERK>::~prolongate_3d_rf2() {}

template <centering_t CENTI, centering_t CENTJ, centering_t CENTK,
          interpolation_t INTPI, interpolation_t INTPJ, interpolation_t INTPK,
          int ORDERI, int ORDERJ, int ORDERK>
amrex::Box prolongate_3d_rf2<CENTI, CENTJ, CENTK, INTPI, INTPJ, INTPK, ORDERI,
                             ORDERJ, ORDERK>::CoarseBox(const amrex::Box &fine,
                                                        int ratio) {
  return CoarseBox(fine, amrex::IntVect(ratio));
}

template <centering_t CENTI, centering_t CENTJ, centering_t CENTK,
          interpolation_t INTPI, interpolation_t INTPJ, interpolation_t INTPK,
          int ORDERI, int ORDERJ, int ORDERK>
amrex::Box
prolongate_3d_rf2<CENTI, CENTJ, CENTK, INTPI, INTPJ, INTPK, ORDERI, ORDERJ,
                  ORDERK>::CoarseBox(const amrex::Box &fine,
                                     const amrex::IntVect &ratio) {
#warning "TODO"
#if 0
  bool error = false;
  for (int d = 0; d < dim; ++d)
    error |= !(fine.type(d) == (indextype()[d] == 0 ? amrex::IndexType::NODE
                                                    : amrex::IndexType::CELL));
  if (error)
    CCTK_VERROR("CENT:%d%d%d CONS:%d%d%d ORDER:%d%d%d fine.type:%d%d%d "
                "indextype:%d%d%d",
                int(CENTI), int(CENTJ), int(CENTK), int(CONSI), int(CONSJ),
                int(CONSK), int(ORDERI), int(ORDERJ), int(ORDERK),
                int(fine.type(0) == amrex::IndexType::CELL),
                int(fine.type(1) == amrex::IndexType::CELL),
                int(fine.type(2) == amrex::IndexType::CELL),
                int(indextype()[0]), int(indextype()[1]), int(indextype()[2]));
  for (int d = 0; d < dim; ++d)
    assert(fine.type(d) == (indextype()[d] == 0 ? amrex::IndexType::NODE
                                                : amrex::IndexType::CELL));
#endif
  for (int d = 0; d < dim; ++d)
    assert(ratio.getVect()[d] == 2);
  amrex::Box crse = amrex::coarsen(fine, 2);
  for (int d = 0; d < dim; ++d)
    crse = crse.grow(d, (order()[d] + 1) / 2);
  return crse;
}

template <centering_t CENTI, centering_t CENTJ, centering_t CENTK,
          interpolation_t INTPI, interpolation_t INTPJ, interpolation_t INTPK,
          int ORDERI, int ORDERJ, int ORDERK>
void prolongate_3d_rf2<CENTI, CENTJ, CENTK, INTPI, INTPJ, INTPK, ORDERI, ORDERJ,
                       ORDERK>::interp(const amrex::FArrayBox &crse,
                                       int crse_comp, amrex::FArrayBox &fine,
                                       int fine_comp, int ncomp,
                                       const amrex::Box &fine_region,
                                       const amrex::IntVect &ratio,
                                       const amrex::Geometry &crse_geom,
                                       const amrex::Geometry &fine_geom,
                                       amrex::Vector<amrex::BCRec> const &bcr,
                                       int actual_comp, int actual_state,
                                       amrex::RunOn gpu_or_cpu) {
  DECLARE_CCTK_PARAMETERS;

  static once_flag have_timers;
  static vector<Timer> timers;

  const int thread_num = omp_get_thread_num();

  call_once(have_timers, [&]() {
    const int num_threads = omp_get_num_threads();
    timers.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i) {
      ostringstream buf;
      buf << "prolongate_3d_rf2<CENT=" << CENTI << CENTJ << CENTK
          << ",INTP=" << INTPI << INTPJ << INTPK << ",ORDER=" << ORDERI
          << ORDERJ << ORDERK << ">[thread=" << i << "]";
      timers.emplace_back(buf.str());
    }
  });

  const Timer &timer = timers.at(thread_num);
  Interval interval(timer);

  for (int d = 0; d < dim; ++d)
    assert(ratio.getVect()[d] == 2);
  // ??? assert(gpu_or_cpu == RunOn::Cpu);

  assert(actual_comp == 0);  // ???
  assert(actual_state == 0); // ???

  // Target box is intersection of fine_region and domain of fine
  const amrex::Box target_region = fine_region & fine.box();
  assert(target_region == fine_region);

  // We prolongate first in the x, then y, then the z direction. Each
  // direction changes the target from coarse-plus-ghosts to fine.
  const amrex::Box source_region = CoarseBox(target_region, 2);
  array<amrex::Box, dim> targets;
  for (int d = 0; d < dim; ++d) {
    targets[d] = d == 0 ? source_region : targets[d - 1];
    targets[d].setRange(d, target_region.loVect()[d], target_region.length(d));
  }
  assert(targets[dim - 1] == target_region);

  // Allocate temporary memory
  // `AsyncArray` frees its memory automatically after the kernels are
  // done using it. However, that mechanism is apparently more
  // expensive than a manual `synchronize`.
  static_assert(dim == 3, "");
  amrex::Gpu::DeviceVector<CCTK_REAL> tmp0(ncomp * targets[0].numPts());
  // amrex::Gpu::AsyncArray<CCTK_REAL> tmp0(ncomp * targets[0].numPts());

  // Check that the input values are finite
#ifdef CCTK_DEBUG
  for (int comp = 0; comp < ncomp; ++comp) {
    const CCTK_REAL *restrict const crseptr = crse.dataPtr(crse_comp + comp);
    for (int k = source_region.loVect()[2]; k <= source_region.hiVect()[2];
         ++k) {
      for (int j = source_region.loVect()[1]; j <= source_region.hiVect()[1];
           ++j) {
#pragma omp simd
        for (int i = source_region.loVect()[0]; i <= source_region.hiVect()[0];
             ++i) {
          const amrex::IntVect ind(i, j, k);
          assert(crse.box().contains(ind));
          assert(isfinite(crseptr[crse.box().index(ind)]));
        }
      }
    }
  }
#endif

  // Initialize the result of the x-prolongation with nan
#ifdef CCTK_DEBUG
  for (int i = 0; i < ncomp * targets[0].numPts(); ++i)
    tmp0.data()[i] = 0.0 / 0.0;
#endif

  // Interpolate in the x-direction
  for (int comp = 0; comp < ncomp; ++comp) {
    const CCTK_REAL *restrict crseptr = crse.dataPtr(crse_comp + comp);
    CCTK_REAL *restrict fineptr = &tmp0.data()[comp * targets[0].numPts()];
    interp3d<CENTI, INTPI, ORDERI, /*D*/ 0>(crseptr, crse.box(), fineptr,
                                            targets[0], targets[0]);
  }

  // Check that the result is finite
#ifdef CCTK_DEBUG
  for (int i = 0; i < ncomp * targets[0].numPts(); ++i)
    assert(isfinite(tmp0.data()[i]));
#endif

  // Allocate temporary memory
  amrex::Gpu::DeviceVector<CCTK_REAL> tmp1(ncomp * targets[1].numPts());
  // amrex::Gpu::AsyncArray<CCTK_REAL> tmp1(ncomp * targets[1].numPts());

  // Initialize the result of the y-prolongation with nan
#ifdef CCTK_DEBUG
  for (int i = 0; i < ncomp * targets[1].numPts(); ++i)
    tmp1.data()[i] = 0.0 / 0.0;
#endif

  // Interpolate in the y-direction
  for (int comp = 0; comp < ncomp; ++comp) {
    const CCTK_REAL *restrict crseptr =
        &tmp0.data()[comp * targets[0].numPts()];
    CCTK_REAL *restrict fineptr = &tmp1.data()[comp * targets[1].numPts()];
    interp3d<CENTJ, INTPJ, ORDERJ, /*D*/ 1>(crseptr, targets[0], fineptr,
                                            targets[1], targets[1]);
  }

  // Check that the result is finite
#ifdef CCTK_DEBUG
  for (int i = 0; i < ncomp * targets[1].numPts(); ++i)
    assert(isfinite(tmp1.data()[i]));
#endif

    // Initialize the result of the z-prolongation with nan
#ifdef CCTK_DEBUG
  for (int comp = 0; comp < ncomp; ++comp) {
    CCTK_REAL *restrict fineptr = fine.dataPtr(fine_comp + comp);
    for (int k = target_region.loVect()[2]; k <= target_region.hiVect()[2];
         ++k) {
      for (int j = target_region.loVect()[1]; j <= target_region.hiVect()[1];
           ++j) {
#pragma omp simd
        for (int i = target_region.loVect()[0]; i <= target_region.hiVect()[0];
             ++i) {
          const amrex::IntVect ind(i, j, k);
          assert(fine.box().contains(ind));
          fineptr[fine.box().index(ind)] = 0.0 / 0.0;
        }
      }
    }
  }
#endif

  // Interpolate in the y-direction
  for (int comp = 0; comp < ncomp; ++comp) {
    const CCTK_REAL *restrict crseptr =
        &tmp1.data()[comp * targets[1].numPts()];
    CCTK_REAL *restrict fineptr = fine.dataPtr(fine_comp + comp);
    interp3d<CENTK, INTPK, ORDERK, /*D*/ 2>(crseptr, targets[1], fineptr,
                                            fine.box(), target_region);
  }

  // Check that the result is finite
#ifdef CCTK_DEBUG
  for (int comp = 0; comp < ncomp; ++comp) {
    CCTK_REAL *restrict fineptr = fine.dataPtr(fine_comp + comp);
    for (int k = target_region.loVect()[2]; k <= target_region.hiVect()[2];
         ++k) {
      for (int j = target_region.loVect()[1]; j <= target_region.hiVect()[1];
           ++j) {
#pragma omp simd
        for (int i = target_region.loVect()[0]; i <= target_region.hiVect()[0];
             ++i) {
          const amrex::IntVect ind(i, j, k);
          assert(fine.box().contains(ind));
          assert(isfinite(fineptr[fine.box().index(ind)]));
        }
      }
    }
  }
#endif

#ifdef __CUDACC__
  amrex::Gpu::synchronize();
  AMREX_GPU_ERROR_CHECK();
#endif
}

template <centering_t CENTI, centering_t CENTJ, centering_t CENTK,
          interpolation_t INTPI, interpolation_t INTPJ, interpolation_t INTPK,
          int ORDERI, int ORDERJ, int ORDERK>
void prolongate_3d_rf2<CENTI, CENTJ, CENTK, INTPI, INTPJ, INTPK, ORDERI, ORDERJ,
                       ORDERK>::interp_face(const amrex::FArrayBox &crse,
                                            int crse_comp,
                                            amrex::FArrayBox &fine,
                                            int fine_comp, int ncomp,
                                            const amrex::Box &fine_region,
                                            const amrex::IntVect &ratio,
                                            const amrex::IArrayBox &solve_mask,
                                            const amrex::Geometry &crse_geom,
                                            const amrex::Geometry &fine_geom,
                                            amrex::Vector<amrex::BCRec> const
                                                &bcr,
                                            int bccomp,
                                            amrex::RunOn gpu_or_cpu) {
  // solve_mask; ???
  assert(bccomp == 0); // ???
  interp(crse, crse_comp, fine, fine_comp, ncomp, fine_region, ratio, crse_geom,
         fine_geom, bcr, 0, 0, gpu_or_cpu);
}

////////////////////////////////////////////////////////////////////////////////

// Self-tests
namespace {
// The constructors for these objects run the self-tests. These definitions will
// make them run at start-up. These tests are very cheap.
static test_interp1d<VC, POLY, 1, CCTK_REAL> test_c0_o1;
static test_interp1d<VC, POLY, 3, CCTK_REAL> test_c0_o3;
static test_interp1d<VC, POLY, 5, CCTK_REAL> test_c0_o5;
static test_interp1d<CC, CONS, 0, CCTK_REAL> test_c1_o0;
static test_interp1d<CC, CONS, 2, CCTK_REAL> test_c1_o2;
static test_interp1d<CC, CONS, 4, CCTK_REAL> test_c1_o4;
} // namespace

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c000_o1;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c001_o1;
prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c010_o1;
prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c011_o1;
prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c100_o1;
prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c101_o1;
prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c110_o1;
prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_3d_rf2_c111_o1;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c000_o3;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c001_o3;
prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c010_o3;
prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c011_o3;
prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c100_o3;
prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c101_o3;
prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c110_o3;
prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_3d_rf2_c111_o3;

prolongate_3d_rf2<VC, VC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c000_o0;
prolongate_3d_rf2<VC, VC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c001_o0;
prolongate_3d_rf2<VC, CC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c010_o0;
prolongate_3d_rf2<VC, CC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c011_o0;
prolongate_3d_rf2<CC, VC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c100_o0;
prolongate_3d_rf2<CC, VC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c101_o0;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c110_o0;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_cons_3d_rf2_c111_o0;

prolongate_3d_rf2<VC, VC, VC, CONS, CONS, CONS, 1, 1, 1>
    prolongate_cons_3d_rf2_c000_o1;
prolongate_3d_rf2<VC, VC, CC, CONS, CONS, CONS, 1, 1, 2>
    prolongate_cons_3d_rf2_c001_o1;
prolongate_3d_rf2<VC, CC, VC, CONS, CONS, CONS, 1, 2, 1>
    prolongate_cons_3d_rf2_c010_o1;
prolongate_3d_rf2<VC, CC, CC, CONS, CONS, CONS, 1, 2, 2>
    prolongate_cons_3d_rf2_c011_o1;
prolongate_3d_rf2<CC, VC, VC, CONS, CONS, CONS, 2, 1, 1>
    prolongate_cons_3d_rf2_c100_o1;
prolongate_3d_rf2<CC, VC, CC, CONS, CONS, CONS, 2, 1, 2>
    prolongate_cons_3d_rf2_c101_o1;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, CONS, 2, 2, 1>
    prolongate_cons_3d_rf2_c110_o1;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 2, 2, 2>
    prolongate_cons_3d_rf2_c111_o1;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1>
    prolongate_ddf_3d_rf2_c000_o1;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 0>
    prolongate_ddf_3d_rf2_c001_o1;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 0, 1>
    prolongate_ddf_3d_rf2_c010_o1;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 0, 0>
    prolongate_ddf_3d_rf2_c011_o1;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 0, 1, 1>
    prolongate_ddf_3d_rf2_c100_o1;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 0, 1, 0>
    prolongate_ddf_3d_rf2_c101_o1;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 0, 0, 1>
    prolongate_ddf_3d_rf2_c110_o1;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0>
    prolongate_ddf_3d_rf2_c111_o1;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3>
    prolongate_ddf_3d_rf2_c000_o3;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 3, 3, 2>
    prolongate_ddf_3d_rf2_c001_o3;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 3, 2, 3>
    prolongate_ddf_3d_rf2_c010_o3;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 3, 2, 2>
    prolongate_ddf_3d_rf2_c011_o3;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 2, 3, 3>
    prolongate_ddf_3d_rf2_c100_o3;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 2, 3, 2>
    prolongate_ddf_3d_rf2_c101_o3;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 2, 2, 3>
    prolongate_ddf_3d_rf2_c110_o3;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 2, 2, 2>
    prolongate_ddf_3d_rf2_c111_o3;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5>
    prolongate_ddf_3d_rf2_c000_o5;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 5, 5, 4>
    prolongate_ddf_3d_rf2_c001_o5;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 5, 4, 5>
    prolongate_ddf_3d_rf2_c010_o5;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 5, 4, 4>
    prolongate_ddf_3d_rf2_c011_o5;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 4, 5, 5>
    prolongate_ddf_3d_rf2_c100_o5;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 4, 5, 4>
    prolongate_ddf_3d_rf2_c101_o5;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 4, 4, 5>
    prolongate_ddf_3d_rf2_c110_o5;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 4, 4, 4>
    prolongate_ddf_3d_rf2_c111_o5;

} // namespace CarpetX
