#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include <defs.hxx>
#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>

namespace Z4c {
using namespace Arith;

#if 0
template <typename T, int D, symm_t symm>
gmat<T, D - 1, symm> calc_submatrix(const gmat<T, D, symm> A, int i0, int j0) {
  return gmat<T, D - 1, symm>([&](int is, int js) {
    const int i = is + (is >= i0);
    const int j = js + (js >= j0);
    return A(i, j);
  });
}

template <typename T, int D, symm_t symm>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const gmat<T, D, symm> &g) {
  if constexpr (D == 0) {
    return one<T>();
  } else {
    T detg = zero<T>();
    for (int i = 0; i < D; ++i)
      detg += bitsign(i) * calc_det(calc_submatrix(g, i, 0));
    return detg;
  }
}

template <typename T, int D, symm_t symm>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr gmat<T, D, symm>
calc_inv(const gmat<T, D, symm> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return gmat<T, D, symm>([&](int i, int j) {
    return bitsign(i) * bitsign(j) * detg1 * calc_det(calc_submatrix(g, j, i));
  });
}
#endif

#if 0

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const smat<T, 0> &g) {
  return one<T>();
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const smat<T, 1> &g) {
  return g(0, 0);
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const smat<T, 2> &g) {
  return g(0, 0) * g(1, 1) - pow2(g(0, 1));
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const smat<T, 3> &g) {
  return 2 * g(0, 1) * g(0, 2) * g(1, 2) - g(2, 2) * pow2(g(0, 1)) -
         g(1, 1) * pow2(g(0, 2)) +
         g(0, 0) * (g(1, 1) * g(2, 2) - pow2(g(1, 2)));
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const smat<T, 4> &g) {
  return 2 * g(0, 0) * g(1, 2) * g(1, 3) * g(2, 3) -
         2 * g(0, 3) *
             (g(0, 2) * (g(1, 2) * g(1, 3) - g(1, 1) * g(2, 3)) +
              g(0, 1) * (-(g(1, 3) * g(2, 2)) + g(1, 2) * g(2, 3))) +
         g(0, 0) * g(1, 1) * g(2, 2) * g(3, 3) +
         2 * g(0, 1) * g(0, 2) * (-(g(1, 3) * g(2, 3)) + g(1, 2) * g(3, 3)) -
         g(2, 2) * g(3, 3) * pow2(g(0, 1)) - g(0, 0) * g(3, 3) * pow2(g(1, 2)) +
         pow2(g(0, 3)) * (-(g(1, 1) * g(2, 2)) + pow2(g(1, 2))) -
         g(0, 0) * g(2, 2) * pow2(g(1, 3)) +
         pow2(g(0, 2)) * (-(g(1, 1) * g(3, 3)) + pow2(g(1, 3))) -
         g(0, 0) * g(1, 1) * pow2(g(2, 3)) + pow2(g(0, 1)) * pow2(g(2, 3));
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const mat<T, 0> &g) {
  return one<T>();
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const mat<T, 1> &g) {
  return g(0, 0);
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const mat<T, 2> &g) {
  return -(g(0, 1) * g(1, 0)) + g(0, 0) * g(1, 1);
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const mat<T, 3> &g) {
  return g(0, 2) * (-(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1)) +
         g(0, 1) * (g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2)) +
         g(0, 0) * (-(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2));
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_det(const mat<T, 4> &g) {
  return g(0, 1) * g(1, 3) * g(2, 2) * g(3, 0) -
         g(0, 1) * g(1, 2) * g(2, 3) * g(3, 0) -
         g(0, 0) * g(1, 3) * g(2, 2) * g(3, 1) +
         g(0, 0) * g(1, 2) * g(2, 3) * g(3, 1) -
         g(0, 1) * g(1, 3) * g(2, 0) * g(3, 2) +
         g(0, 0) * g(1, 3) * g(2, 1) * g(3, 2) +
         g(0, 1) * g(1, 0) * g(2, 3) * g(3, 2) -
         g(0, 0) * g(1, 1) * g(2, 3) * g(3, 2) +
         g(0, 3) * (g(1, 2) * (g(2, 1) * g(3, 0) - g(2, 0) * g(3, 1)) +
                    g(1, 1) * (-(g(2, 2) * g(3, 0)) + g(2, 0) * g(3, 2)) +
                    g(1, 0) * (g(2, 2) * g(3, 1) - g(2, 1) * g(3, 2))) +
         g(0, 1) * g(1, 2) * g(2, 0) * g(3, 3) -
         g(0, 0) * g(1, 2) * g(2, 1) * g(3, 3) -
         g(0, 1) * g(1, 0) * g(2, 2) * g(3, 3) +
         g(0, 0) * g(1, 1) * g(2, 2) * g(3, 3) +
         g(0, 2) * (g(1, 3) * (-(g(2, 1) * g(3, 0)) + g(2, 0) * g(3, 1)) +
                    g(1, 1) * (g(2, 3) * g(3, 0) - g(2, 0) * g(3, 3)) +
                    g(1, 0) * (-(g(2, 3) * g(3, 1)) + g(2, 1) * g(3, 3)));
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 0>
calc_inv(const smat<T, 0> &g, const T &detg) {
  return smat<T, 0>([&](int i, int j) ARITH_INLINE { return one<T>()(); });
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 1>
calc_inv(const smat<T, 1> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return smat<T, 1>{detg1};
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 2>
calc_inv(const smat<T, 2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * smat<T, 2>{
                     g(1, 1),
                     -g(0, 1),

                     g(0, 0),
                 };
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 3>
calc_inv(const smat<T, 3> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * smat<T, 3>{
                     g(1, 1) * g(2, 2) - pow2(g(1, 2)),
                     g(0, 2) * g(1, 2) - g(0, 1) * g(2, 2),
                     -(g(0, 2) * g(1, 1)) + g(0, 1) * g(1, 2),
                     g(0, 0) * g(2, 2) - pow2(g(0, 2)),
                     g(0, 1) * g(0, 2) - g(0, 0) * g(1, 2),
                     g(0, 0) * g(1, 1) - pow2(g(0, 1)),
                 };
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 4>
calc_inv(const smat<T, 4> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * smat<T, 4>{
                     2 * g(1, 2) * g(1, 3) * g(2, 3) - g(3, 3) * pow2(g(1, 2)) -
                         g(2, 2) * pow2(g(1, 3)) +
                         g(1, 1) * (g(2, 2) * g(3, 3) - pow2(g(2, 3))),
                     g(0, 3) * (g(1, 3) * g(2, 2) - g(1, 2) * g(2, 3)) +
                         g(0, 2) * (-(g(1, 3) * g(2, 3)) + g(1, 2) * g(3, 3)) +
                         g(0, 1) * (-(g(2, 2) * g(3, 3)) + pow2(g(2, 3))),
                     g(0, 3) * (-(g(1, 2) * g(1, 3)) + g(1, 1) * g(2, 3)) +
                         g(0, 1) * (-(g(1, 3) * g(2, 3)) + g(1, 2) * g(3, 3)) +
                         g(0, 2) * (-(g(1, 1) * g(3, 3)) + pow2(g(1, 3))),
                     g(0, 2) * (-(g(1, 2) * g(1, 3)) + g(1, 1) * g(2, 3)) +
                         g(0, 1) * (g(1, 3) * g(2, 2) - g(1, 2) * g(2, 3)) +
                         g(0, 3) * (-(g(1, 1) * g(2, 2)) + pow2(g(1, 2))),

                     2 * g(0, 2) * g(0, 3) * g(2, 3) - g(3, 3) * pow2(g(0, 2)) -
                         g(2, 2) * pow2(g(0, 3)) +
                         g(0, 0) * (g(2, 2) * g(3, 3) - pow2(g(2, 3))),
                     -(g(0, 3) * (g(0, 2) * g(1, 3) + g(0, 1) * g(2, 3))) +
                         g(0, 1) * g(0, 2) * g(3, 3) +
                         g(0, 0) * (g(1, 3) * g(2, 3) - g(1, 2) * g(3, 3)) +
                         g(1, 2) * pow2(g(0, 3)),
                     g(0, 1) * g(0, 3) * g(2, 2) -
                         g(0, 2) * (g(0, 3) * g(1, 2) + g(0, 1) * g(2, 3)) +
                         g(0, 0) * (-(g(1, 3) * g(2, 2)) + g(1, 2) * g(2, 3)) +
                         g(1, 3) * pow2(g(0, 2)),

                     2 * g(0, 1) * g(0, 3) * g(1, 3) - g(3, 3) * pow2(g(0, 1)) -
                         g(1, 1) * pow2(g(0, 3)) +
                         g(0, 0) * (g(1, 1) * g(3, 3) - pow2(g(1, 3))),
                     -(g(0, 1) * g(0, 3) * g(1, 2)) +
                         g(0, 2) * (g(0, 3) * g(1, 1) - g(0, 1) * g(1, 3)) +
                         g(0, 0) * (g(1, 2) * g(1, 3) - g(1, 1) * g(2, 3)) +
                         g(2, 3) * pow2(g(0, 1)),

                     2 * g(0, 1) * g(0, 2) * g(1, 2) - g(2, 2) * pow2(g(0, 1)) -
                         g(1, 1) * pow2(g(0, 2)) +
                         g(0, 0) * (g(1, 1) * g(2, 2) - pow2(g(1, 2))),
                 };
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 0>
calc_inv(const mat<T, 0> &g, const T &detg) {
  return mat<T, 0>([&](int i, int j) ARITH_INLINE { return one<T>()(); });
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 1>
calc_inv(const mat<T, 1> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return mat<T, 1>{detg1};
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 2>
calc_inv(const mat<T, 2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * mat<T, 2>{
                     g(1, 1),
                     -g(0, 1),

                     -g(1, 0),
                     g(0, 0),
                 };
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 3>
calc_inv(const mat<T, 3> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * mat<T, 3>{
                     -(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2),
                     g(0, 2) * g(2, 1) - g(0, 1) * g(2, 2),
                     -(g(0, 2) * g(1, 1)) + g(0, 1) * g(1, 2),

                     g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2),
                     -(g(0, 2) * g(2, 0)) + g(0, 0) * g(2, 2),
                     g(0, 2) * g(1, 0) - g(0, 0) * g(1, 2),

                     -(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1),
                     g(0, 1) * g(2, 0) - g(0, 0) * g(2, 1),
                     -(g(0, 1) * g(1, 0)) + g(0, 0) * g(1, 1),
                 };
}

template <typename T>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 4>
calc_inv(const mat<T, 4> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * mat<T, 4>{
                     g(1, 3) * (-(g(2, 2) * g(3, 1)) + g(2, 1) * g(3, 2)) +
                         g(1, 2) * (g(2, 3) * g(3, 1) - g(2, 1) * g(3, 3)) +
                         g(1, 1) * (-(g(2, 3) * g(3, 2)) + g(2, 2) * g(3, 3)),
                     g(0, 3) * (g(2, 2) * g(3, 1) - g(2, 1) * g(3, 2)) +
                         g(0, 2) * (-(g(2, 3) * g(3, 1)) + g(2, 1) * g(3, 3)) +
                         g(0, 1) * (g(2, 3) * g(3, 2) - g(2, 2) * g(3, 3)),
                     g(0, 3) * (-(g(1, 2) * g(3, 1)) + g(1, 1) * g(3, 2)) +
                         g(0, 2) * (g(1, 3) * g(3, 1) - g(1, 1) * g(3, 3)) +
                         g(0, 1) * (-(g(1, 3) * g(3, 2)) + g(1, 2) * g(3, 3)),
                     g(0, 3) * (g(1, 2) * g(2, 1) - g(1, 1) * g(2, 2)) +
                         g(0, 2) * (-(g(1, 3) * g(2, 1)) + g(1, 1) * g(2, 3)) +
                         g(0, 1) * (g(1, 3) * g(2, 2) - g(1, 2) * g(2, 3)),

                     g(1, 3) * (g(2, 2) * g(3, 0) - g(2, 0) * g(3, 2)) +
                         g(1, 2) * (-(g(2, 3) * g(3, 0)) + g(2, 0) * g(3, 3)) +
                         g(1, 0) * (g(2, 3) * g(3, 2) - g(2, 2) * g(3, 3)),
                     g(0, 3) * (-(g(2, 2) * g(3, 0)) + g(2, 0) * g(3, 2)) +
                         g(0, 2) * (g(2, 3) * g(3, 0) - g(2, 0) * g(3, 3)) +
                         g(0, 0) * (-(g(2, 3) * g(3, 2)) + g(2, 2) * g(3, 3)),
                     g(0, 3) * (g(1, 2) * g(3, 0) - g(1, 0) * g(3, 2)) +
                         g(0, 2) * (-(g(1, 3) * g(3, 0)) + g(1, 0) * g(3, 3)) +
                         g(0, 0) * (g(1, 3) * g(3, 2) - g(1, 2) * g(3, 3)),
                     g(0, 3) * (-(g(1, 2) * g(2, 0)) + g(1, 0) * g(2, 2)) +
                         g(0, 2) * (g(1, 3) * g(2, 0) - g(1, 0) * g(2, 3)) +
                         g(0, 0) * (-(g(1, 3) * g(2, 2)) + g(1, 2) * g(2, 3)),

                     g(1, 3) * (-(g(2, 1) * g(3, 0)) + g(2, 0) * g(3, 1)) +
                         g(1, 1) * (g(2, 3) * g(3, 0) - g(2, 0) * g(3, 3)) +
                         g(1, 0) * (-(g(2, 3) * g(3, 1)) + g(2, 1) * g(3, 3)),
                     g(0, 3) * (g(2, 1) * g(3, 0) - g(2, 0) * g(3, 1)) +
                         g(0, 1) * (-(g(2, 3) * g(3, 0)) + g(2, 0) * g(3, 3)) +
                         g(0, 0) * (g(2, 3) * g(3, 1) - g(2, 1) * g(3, 3)),
                     g(0, 3) * (-(g(1, 1) * g(3, 0)) + g(1, 0) * g(3, 1)) +
                         g(0, 1) * (g(1, 3) * g(3, 0) - g(1, 0) * g(3, 3)) +
                         g(0, 0) * (-(g(1, 3) * g(3, 1)) + g(1, 1) * g(3, 3)),
                     g(0, 3) * (g(1, 1) * g(2, 0) - g(1, 0) * g(2, 1)) +
                         g(0, 1) * (-(g(1, 3) * g(2, 0)) + g(1, 0) * g(2, 3)) +
                         g(0, 0) * (g(1, 3) * g(2, 1) - g(1, 1) * g(2, 3)),

                     g(1, 2) * (g(2, 1) * g(3, 0) - g(2, 0) * g(3, 1)) +
                         g(1, 1) * (-(g(2, 2) * g(3, 0)) + g(2, 0) * g(3, 2)) +
                         g(1, 0) * (g(2, 2) * g(3, 1) - g(2, 1) * g(3, 2)),
                     g(0, 2) * (-(g(2, 1) * g(3, 0)) + g(2, 0) * g(3, 1)) +
                         g(0, 1) * (g(2, 2) * g(3, 0) - g(2, 0) * g(3, 2)) +
                         g(0, 0) * (-(g(2, 2) * g(3, 1)) + g(2, 1) * g(3, 2)),
                     g(0, 2) * (g(1, 1) * g(3, 0) - g(1, 0) * g(3, 1)) +
                         g(0, 1) * (-(g(1, 2) * g(3, 0)) + g(1, 0) * g(3, 2)) +
                         g(0, 0) * (g(1, 2) * g(3, 1) - g(1, 1) * g(3, 2)),
                     g(0, 2) * (-(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1)) +
                         g(0, 1) * (g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2)) +
                         g(0, 0) * (-(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2)),
                 };
}

#endif

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_trace(const smat<T, D> &A, const smat<T, D> &gu) {
  return sum_symm<D>([&](int x, int y)
                         ARITH_INLINE { return gu(x, y) * A(x, y); });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<vec<T, D>, D>
calc_dgu(const smat<T, D> &gu, const smat<vec<T, D>, D> &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  return smat<vec<T, D>, D>([&](int a, int b) ARITH_INLINE {
    return vec<T, D>([&](int c) ARITH_INLINE {
      // return sum2sym([&](int x, int y) ARITH_INLINE {
      //   return -gu(a, x) * gu(b, y) * dg(x, y)(c);
      // });
      return sum<D>([&](int x) ARITH_INLINE {
        return -gu(a, x) * sum<D>([&](int y) ARITH_INLINE {
          return gu(b, y) * dg(x, y)(c);
        });
      });
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<vec<T, D>, D>
calc_dAu(const smat<T, D> &gu, const smat<vec<T, D>, D> &dgu,
         const smat<T, D> &A, const smat<vec<T, D>, D> &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  return smat<vec<T, D>, D>([&](int a, int b) ARITH_INLINE {
    return vec<T, D>([&](int c) ARITH_INLINE {
      // return sum2sym([&](int x, int y) ARITH_INLINE {
      //   return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
      //          + gu(a, x) * dgu(b, y)(c) * A(x, y) //
      //          + gu(a, x) * gu(b, y) * dA(x, y)(c);
      // });
      return sum<D>([&](int x) ARITH_INLINE {
        return gu(b, x) * sum<D>([&](int y) ARITH_INLINE {
                 return dgu(a, y)(c) * A(x, y)   //
                        + dgu(b, y)(c) * A(x, y) //
                        + gu(b, y) * dA(x, y)(c);
               });
      });
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<smat<T, D>, D>
calc_gammal(const smat<vec<T, D>, D> &dg) {
  // Gammal_abc
  return vec<smat<T, D>, D>([&](int a) ARITH_INLINE {
    return smat<T, D>([&](int b, int c) ARITH_INLINE {
      return (dg(a, b)(c) + dg(a, c)(b) - dg(b, c)(a)) / 2;
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<smat<T, D>, D>
calc_gamma(const smat<T, D> &gu, const vec<smat<T, D>, D> &Gammal) {
  // Gamma^a_bc
  return vec<smat<T, D>, D>([&](int a) ARITH_INLINE {
    return smat<T, D>([&](int b, int c) ARITH_INLINE {
      return sum<D>([&](int x)
                        ARITH_INLINE { return gu(a, x) * Gammal(x)(b, c); });
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<vec<vec<T, D>, D>, D>
calc_gammalu(const smat<T, D> &gu, const vec<smat<T, D>, D> &Gammal) {
  // Gamma_ab^c
  return vec<vec<vec<T, D>, D>, D>([&](int a) ARITH_INLINE {
    return vec<vec<T, D>, D>([&](int b) ARITH_INLINE {
      return vec<T, D>([&](int c) ARITH_INLINE {
        return sum<D>([&](int x)
                          ARITH_INLINE { return Gammal(a)(b, x) * gu(x, c); });
      });
    });
  });
}

} // namespace Z4c

#endif // #ifndef PHYSICS_HXX
