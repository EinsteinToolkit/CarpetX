// TODO: Don't include files from other thorn; create a proper interface
#include "../../CarpetX/src/driver.hxx"

#include <pdesolvers.hxx>

#include <defs.hxx>
#include <dual.hxx>
#include <loop.hxx>
#include <smallvector.hxx>
#include <spvect.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <utility>

namespace Poisson2 {

constexpr int dim = Loop::dim;

extern "C" void Poisson2_Source(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Source;
  DECLARE_CCTK_PARAMETERS;

  const int npoints = 27; // 3 levels

  const std::array<int, dim> indextype = {0, 0, 0};
  const Loop::GF3D2layout layout1(cctkGH, indextype);
  const Loop::GridDescBase grid(cctkGH);

  const Loop::GF3D2<CCTK_REAL> gf_src(layout1, src);

  if (CCTK_EQUALS(source, "constant")) {
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [=](const Loop::PointDesc &p)
                               ARITH_INLINE { gf_src(p.I) = 1; });
  } else if (CCTK_EQUALS(source, "logo")) {
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
          using std::lrint;
          Arith::vect<int, dim> i;
          for (int d = 0; d < dim; ++d)
            i[d] = lrint((p.X[d] / logo_width + 0.5) * npoints - 0.5);
          if (all(i >= Arith::vect<int, dim>::pure(0) &&
                  i < Arith::vect<int, dim>::pure(npoints))) {
            if (all(i / 9 == Arith::vect<int, dim>::pure(1))) {
              gf_src(p.I) = 3; // red box (large)
            } else if (all(i / 3 % 3 == Arith::vect<int, dim>::pure(1))) {
              gf_src(p.I) = 2; // green box (medium)
            } else if (all(i / 1 % 3 == Arith::vect<int, dim>::pure(1))) {
              gf_src(p.I) = 1; // blue box (small)
            } else {
              gf_src(p.I) = 0;
            }
          } else {
            gf_src(p.I) = 0;
          }
        });
  } else {
    CCTK_ERROR("Unknown value for parameter \"source\"");
  }
}

extern "C" void Poisson2_ChooseInitialGuess(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_ChooseInitialGuess;

  // Set up an initial guess only when initializing the coarsest
  // level. Otherwise, keep the data from the previous solve.
  *need_initial_guess = CarpetX::ghext->num_levels() == 1;
}

extern "C" void Poisson2_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Init;

  const std::array<int, dim> indextype = {0, 0, 0};
  const Loop::GF3D2layout layout1(cctkGH, indextype);
  const Loop::GridDescBase grid(cctkGH);

  const Loop::GF3D2<CCTK_REAL> gf_sol(layout1, sol);

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p)
                             ARITH_INLINE { gf_sol(p.I) = 0; });
}

extern "C" void Poisson2_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Sync;
  // do nothing
}

extern "C" void Poisson2_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Boundaries;

  const std::array<int, dim> indextype = {0, 0, 0};
  const Loop::GF3D2layout layout1(cctkGH, indextype);
  const Loop::GridDescBase grid(cctkGH);

  const Loop::GF3D2<CCTK_REAL> gf_sol(layout1, sol);

  grid.loop_bnd<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p)
                             ARITH_INLINE { gf_sol(p.I) = 0; });
}

constexpr std::size_t max_small_size = 16;
template <typename T>
using smallvector1 = Arith::smallvector<T, max_small_size>;
template <typename T>
using spvect1 = Arith::spvect<Arith::vect<int, dim>, T, smallvector1>;
template <typename T> using dual1 = Arith::dual<T, spvect1<T> >;

template <typename T>
dual1<T> residual(const int fdorder, const Loop::GF3D2<const T> &point_type,
                  const Loop::GF3D2<const T> &sol_,
                  const Loop::GF3D2<const T> &src, const Loop::PointDesc &p) {
  const auto sol = [&](const Arith::vect<int, dim> &I) {
    const T val = sol_(I);
    spvect1<T> eps;
    eps.emplace_back(I, 1);
    return dual1<T>(val, std::move(eps));
  };

  assert(point_type(p.I) == 1 /*intr*/ || point_type(p.I) == 5 /*rest*/);
  using Arith::pow2;
  dual1<T> ddsol = 0;
  switch (fdorder) {
  case 2:
    for (int d = 0; d < dim; ++d)
      ddsol += (sol(p.I - p.DI[d]) //
                - 2 * sol(p.I)     //
                + sol(p.I + p.DI[d])) /
               pow2(p.DX[d]);
    break;
  case 4:
    for (int d = 0; d < dim; ++d)
      ddsol += (-1 / 12.0 * sol(p.I - 2 * p.DI[d]) //
                + 4 / 3.0 * sol(p.I - p.DI[d])     //
                - 5 / 2.0 * sol(p.I)               //
                + 4 / 3.0 * sol(p.I + p.DI[d])     //
                - 1 / 12.0 * sol(p.I + 2 * p.DI[d])) /
               pow2(p.DX[d]);
    break;
  default:
    assert(0);
    abort();
  }
  return ddsol - src(p.I);
}

extern "C" void Poisson2_Residual(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Residual;
  DECLARE_CCTK_PARAMETERS;

  const std::array<int, dim> indextype = {0, 0, 0};
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  const Loop::GF3D2layout layout1(cctkGH, indextype);
  const Loop::GridDescBase grid(cctkGH);

  for (int d = 0; d < dim; ++d)
    assert(2 * nghostzones[d] >= fdorder);

  const Loop::GF3D2<const CCTK_REAL> gf_point_type(layout1, point_type);
  const Loop::GF3D2<const CCTK_REAL> gf_sol(layout1, sol);
  const Loop::GF3D2<const CCTK_REAL> gf_src(layout1, src);
  const Loop::GF3D2<CCTK_REAL> gf_res(layout1, res);

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
#if 0
        assert(gf_point_type(p.I) == 1);
        using Arith::pow2;
        CCTK_REAL ddsol = 0;
        switch (fdorder) {
        case 2:
          for (int d = 0; d < dim; ++d)
            ddsol += (gf_sol(p.I - p.DI[d]) //
                      - 2 * gf_sol(p.I)     //
                      + gf_sol(p.I + p.DI[d])) /
                     pow2(p.DX[d]);
          break;
        case 4:
          for (int d = 0; d < dim; ++d)
            ddsol += (-1 / 12.0 * gf_sol(p.I - 2 * p.DI[d]) //
                      + 4 / 3.0 * gf_sol(p.I - p.DI[d])     //
                      - 5 / 2.0 * gf_sol(p.I)               //
                      + 4 / 3.0 * gf_sol(p.I + p.DI[d])     //
                      - 1 / 12.0 * gf_sol(p.I + 2 * p.DI[d])) /
                     pow2(p.DX[d]);
          break;
        default:
          assert(0);
          abort();
        }
        gf_res(p.I) = ddsol - gf_src(p.I);
#else
        gf_res(p.I) = residual(fdorder,gf_point_type,gf_sol,gf_src,p).val;
#endif
                         });
}

extern "C" void Poisson2_ResidualSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_ResidualSync;
  // do nothing
}

extern "C" void Poisson2_ResidualBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_ResidualBoundaries;
  DECLARE_CCTK_PARAMETERS;

  const std::array<int, dim> indextype = {0, 0, 0};
  const Loop::GF3D2layout layout1(cctkGH, indextype);
  const Loop::GridDescBase grid(cctkGH);

  const Loop::GF3D2<const CCTK_REAL> gf_point_type(layout1, point_type);
  const Loop::GF3D2<CCTK_REAL> gf_res(layout1, res);

  grid.loop_bnd<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
                           assert(gf_point_type(p.I) == 2);
                           gf_res(p.I) = 0;
                         });
}

extern "C" void Poisson2_Jacobian(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Poisson2_Jacobian;
  DECLARE_CCTK_PARAMETERS;

  const std::array<int, dim> indextype = {0, 0, 0};
  const Loop::GF3D2layout layout1(cctkGH, indextype);
  const Loop::GridDescBase grid(cctkGH);

  const Loop::GF3D2<const CCTK_REAL> gf_point_type(layout1, point_type);
  const Loop::GF3D2<const CCTK_REAL> gf_idx(layout1, idx);
  const Loop::GF3D2<const CCTK_REAL> gf_sol(layout1, sol);
  const Loop::GF3D2<const CCTK_REAL> gf_src(layout1, src);

  assert(PDESolvers::jacobians.has_value());
  PDESolvers::jacobian_t &J = PDESolvers::jacobians.value().get_local();

  const int fdorder1 = fdorder_jac < 0 ? fdorder : fdorder_jac;

  grid.loop_int<0, 0, 0>(
      grid.nghostzones, [=, &J](const Loop::PointDesc &p) ARITH_INLINE {
        const auto getidx = [&](const Arith::vect<int, dim> &I) {
          return int(gf_idx(I));
        };
        const int idx = getidx(p.I);
#if 0
        assert(gf_point_type(p.I) == 1);
        using Arith::pow2;
        switch (fdorder1) {
        case 2:
          for (int d = 0; d < dim; ++d) {
            assert(getidx(p.I - p.DI[d]) >= 0);
            assert(getidx(p.I + p.DI[d]) >= 0);
            J.add_value(idx, getidx(p.I - p.DI[d]), 1 / pow2(p.DX[d]));
            J.add_value(idx, getidx(p.I), -2 / pow2(p.DX[d]));
            J.add_value(idx, getidx(p.I + p.DI[d]), 1 / pow2(p.DX[d]));
          }
          break;
        case 4:
          for (int d = 0; d < dim; ++d) {
            J.add_value(idx, getidx(p.I - 2 * p.DI[d]),
                        -1 / 12.0 / pow2(p.DX[d]));
            J.add_value(idx, getidx(p.I - p.DI[d]), 4 / 3.0 / pow2(p.DX[d]));
            J.add_value(idx, getidx(p.I), -5 / 2.0 / pow2(p.DX[d]));
            J.add_value(idx, getidx(p.I + p.DI[d]), 4 / 3.0 / pow2(p.DX[d]));
            J.add_value(idx, getidx(p.I + 2 * p.DI[d]),
                        -1 / 12.0 / pow2(p.DX[d]));
          }
          break;
        default:
          assert(0);
          abort();
        }
#else
        const auto row =
            residual(fdorder1, gf_point_type, gf_sol, gf_src, p).eps;
        for (const auto &elt : row)
          J.add_value(idx, getidx(elt.first), elt.second);
#endif
      });
}

} // namespace Poisson2
