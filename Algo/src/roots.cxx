#include "roots.hxx"

#include <mat.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <cassert>
#include <cmath>
#include <limits>
#include <tuple>

namespace Algo {

using Arith::mat, Arith::vec;

namespace {
template <typename T> CCTK_DEVICE vec<T, 2> gn(vec<T, 2> x) {
  return vec<T, 2>{x(0) * x(0) - 2, x(0) * x(1) - 2};
}
template <typename T> CCTK_DEVICE
std::pair<vec<T, 2>, mat<T, 2> > gnd(vec<T, 2> x) {
  return {vec<T, 2>{x(0) * x(0) - 2, x(0) * x(1) - 2},
          mat<T, 2>{2 * x(0), x(1), 0, x(0)}};
}
} // namespace

extern "C" void Test_roots(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  auto fn = []CCTK_HOST CCTK_DEVICE (CCTK_REAL x) { return x * x - 2; };
  auto fnd = []CCTK_HOST CCTK_DEVICE (CCTK_REAL x) {
    return std::make_tuple(x * x - 2, 2 * x); };
  auto fnd2 = []CCTK_HOST CCTK_DEVICE (CCTK_REAL x) {
    return std::make_tuple(x * x - 2, 2 * x, 2); };

  {
    const int minbits = std::numeric_limits<double>::digits - 4;
    const int maxiters = 100;
    int iters;
    auto [lo, hi] = bisect(fn, 1.0, 2.0, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    assert(hi >= lo && hi - lo <= std::scalbn(2.0, -minbits));
    assert(fn(lo) <= 0 && fn(hi) >= 0);
    CCTK_VINFO("Test_bisect succeeded in %d iterations", iters);
  }

  {
    const int minbits = std::numeric_limits<double>::digits - 4;
    const int maxiters = 100;
    int iters;
    auto [lo, hi] =
        bracket_and_solve_root(fn, 1.0, 2.0, true, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    assert(hi >= lo && hi - lo <= std::scalbn(2.0, -minbits));
    assert(fn(lo) <= 0 && fn(hi) >= 0);
    CCTK_VINFO("Test_bracket_and_solve_root succeeded in %d iterations", iters);
  }

  {
    const int minbits =
        static_cast<int>(0.6 * std::numeric_limits<double>::digits);
    const int maxiters = 100;
    int iters;
    auto x = newton_raphson(fnd, 1.0, 0.0, 10.0, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    double delta = std::scalbn(1.0, -minbits);
    assert(fn(x - delta) * fn(x + delta) < 0);
    CCTK_VINFO("Test_newton_raphson succeeded in %d iterations", iters);
  }

  {
    const int minbits =
        static_cast<int>(0.6 * std::numeric_limits<double>::digits);
    const int maxiters = 100;
    int iters;
    auto x = halley(fnd2, 1.0, 0.0, 10.0, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    double delta = std::scalbn(1.0, -minbits);
    assert(fn(x - delta) * fn(x + delta) < 0);
    CCTK_VINFO("Test_halley succeeded in %d iterations", iters);
  }

  {
    const int minbits =
        static_cast<int>(0.6 * std::numeric_limits<double>::digits);
    const int maxiters = 100;
    int iters;
    auto x = schroder(fnd2, 1.0, 0.0, 10.0, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    double delta = std::scalbn(1.0, -minbits);
    assert(fn(x - delta) * fn(x + delta) < 0);
    CCTK_VINFO("Test_schroder succeeded in %d iterations", iters);
  }

  {
    CCTK_INT * iters_print =
      (CCTK_INT *) amrex::The_Arena()->alloc(sizeof(CCTK_INT));
    const amrex::Box box(amrex::IntVect(0), amrex::IntVect(0, 0, 0));

    amrex::launch(box, [=] CCTK_DEVICE(amrex::Box const& tbx)
        CCTK_ATTRIBUTE_ALWAYS_INLINE {
      const int minbits = std::numeric_limits<double>::digits - 4;
      const int maxiters = 100;
      int iters;
      auto [lo, hi] = brent(fn, 1.0, 2.0, minbits, maxiters, iters);
      // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
      // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
      assert(iters < maxiters);
      assert(hi >= lo && hi - lo <= std::scalbn(2.0, -minbits));
      assert(fn(lo) <= 0 && fn(hi) >= 0);
      iters_print[0] = iters;
    });
    amrex::Gpu::Device::streamSynchronize();
    CCTK_VINFO("Test_brent succeeded in %d iterations", iters_print[0]);

    amrex::launch(box, [=] CCTK_DEVICE(amrex::Box const& tbx)
        CCTK_ATTRIBUTE_ALWAYS_INLINE {
      const int minbits =
          static_cast<int>(0.6 * std::numeric_limits<double>::digits);
      const int maxiters = 100;
      int iters;
      bool failed;
      auto x = newton_raphson_nd(gnd<CCTK_REAL>, vec<CCTK_REAL, 2>{1.0, 1.0},
                                 vec<CCTK_REAL, 2>{0.0, 0.0},
                                 vec<CCTK_REAL, 2>{10.0, 10.0},
                                 minbits, maxiters, iters, failed);
      // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
      // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
      assert(iters < maxiters);
      assert(sumabs(gn(x)) < 1.0e-9);
      iters_print[0] = iters;
    });
    amrex::Gpu::Device::streamSynchronize();
    CCTK_VINFO("Test_newton_raphson_nd succeeded in %d iterations",
        iters_print[0]);

    amrex::The_Arena()->free(iters_print);
  }

}

} // namespace Algo
