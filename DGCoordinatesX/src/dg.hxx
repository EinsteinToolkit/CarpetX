#ifndef DG_HXX
#define DG_HXX

#include <div.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <array>

namespace DG {
using namespace Arith;
using namespace Loop;

template <typename T, int N> struct dg_basis;

// Collocation points are in the range [-1; +1].
// Derivative coefficients include 1 ghost point.

// The derivative operator is non-zero for ghosts only for the boundary points,
// and is zero on the diagonal.

template <typename T> struct dg_basis<T, 2> {
  // Element size
  static constexpr int N = 2;
  // Collocation points
  static constexpr std::array<T, N> coords{-1, 1};
  // Collocation point weights for integration; the weights sum to 2
  static constexpr std::array<T, N> vols{1, 1};
  // Derivative stencils for each collocation point
  using ATN2 = std::array<T, N + 2>;
  static constexpr std::array<ATN2, N> deriv{ATN2{-0.5, 0, 0.5, 0},
                                             ATN2{0, -0.5, 0, 0.5}};
};

template <typename T> struct dg_basis<T, 3> {
  static constexpr int N = 3;
  static constexpr std::array<T, N> coords{-1, 0, 1};
  static constexpr std::array<T, N> vols{0.5, 1, 0.5};
  using ATN2 = std::array<T, N + 2>;
  static constexpr std::array<ATN2, N> deriv{ATN2{-1.5, 0, 2, -0.5, 0}, //
                                             ATN2{0, -0.5, 0, 0.5, 0},  //
                                             ATN2{0, 0.5, -2, 0, 1.5}};
};

template <typename T> struct dg_basis<T, 4> {
  static constexpr int N = 4;
  static constexpr std::array<T, N> coords{
      -1,
      -0.4472135954999579392818347337462552470881,
      0.4472135954999579392818347337462552470881,
      1,
  };
  static constexpr std::array<T, N> vols{
      0.166666666666666666666666666667,
      0.83333333333333333333333333333,
      0.83333333333333333333333333333,
      0.166666666666666666666666666667,
  };
  using ATN2 = std::array<T, N + 2>;
  static constexpr std::array<ATN2, N> deriv{
      ATN2{
          -3,
          0,
          4.0450849718747371205114670859,
          -1.54508497187473712051146708591,
          0.500000000000000000000000000000,
          0,
      },
      ATN2{
          0,
          -0.80901699437494742410229341718,
          0,
          1.1180339887498948482045868344,
          -0.309016994374947424102293417183,
          0,
      },
      ATN2{
          0,
          0.309016994374947424102293417183,
          -1.1180339887498948482045868344,
          0,
          0.80901699437494742410229341718,
          0,
      },
      ATN2{
          0,
          -0.500000000000000000000000000000,
          1.54508497187473712051146708591,
          -4.0450849718747371205114670859,
          0,
          3,
      }};
};

template <typename T> struct dg_basis<T, 8> {
  static constexpr int N = 8;
  static constexpr std::array<T, N> coords{
      -1,
      -0.8717401485096066153374457612206634381038,
      -0.5917001814331423021445107313979531899457,
      -0.2092992179024788687686572603453512552955,
      0.2092992179024788687686572603453512552955,
      0.5917001814331423021445107313979531899457,
      0.8717401485096066153374457612206634381038,
      1,
  };
  static constexpr std::array<T, N> vols{
      0.035714285714285714285714285714, 0.210704227143506039382992065776,
      0.34112269248350436476424067711,  0.41245879465870388156705297140,
      0.41245879465870388156705297140,  0.34112269248350436476424067711,
      0.210704227143506039382992065776, 0.035714285714285714285714285714,
  };
  using ATN2 = std::array<T, N + 2>;
  static constexpr std::array<ATN2, N> deriv{
      ATN2{
          -14,
          0,
          18.9375986071173705131546424941,
          -7.5692898193484870087981787234,
          4.2979081642651752142490066598,
          -2.810188989257949038537763687,
          1.94165942554412230232935244532,
          -1.29768738832023198239705918927,
          0.50000000000000000000000000000,
          0,
      },
      ATN2{
          0,
          -3.209915703002990336815351572,
          0,
          4.5435850645665642173963830608,
          -2.11206121431454222839296010401,
          1.2942320509135015076932056372,
          -0.8694480983314929341612144455,
          0.5735654149402641373099571270,
          -0.21995751477130436303001970391,
          0,
      },
      ATN2{
          0,
          0.7924766813205145268013597361,
          -2.806475794736433436470922316,
          0,
          2.875517405972505217473636209,
          -1.37278583180602848091037564048,
          0.84502255650651048982986561020,
          -0.5370395861576610609376193431,
          0.203284568900592744214055743889,
          0,
      },
      ATN2{
          0,
          -0.3721504357285948658470918627,
          1.0789446887904527025084635997,
          -2.378187233515505824569788251,
          0,
          2.388924359158239214879541510,
          -1.1353580168811114404647212743,
          0.6611573509003112232014379102,
          -0.2433307127237910097078416316,
          0,
      },
      ATN2{
          0,
          0.2433307127237910097078416316,
          -0.6611573509003112232014379102,
          1.1353580168811114404647212743,
          -2.388924359158239214879541510,
          0,
          2.378187233515505824569788251,
          -1.0789446887904527025084635997,
          0.3721504357285948658470918627,
          0,
      },
      ATN2{
          0,
          -0.203284568900592744214055743889,
          0.5370395861576610609376193431,
          -0.84502255650651048982986561020,
          1.37278583180602848091037564048,
          -2.875517405972505217473636209,
          0,
          2.806475794736433436470922316,
          -0.7924766813205145268013597361,
          0,
      },
      ATN2{
          0,
          0.21995751477130436303001970391,
          -0.5735654149402641373099571270,
          0.8694480983314929341612144455,
          -1.2942320509135015076932056372,
          2.11206121431454222839296010401,
          -4.5435850645665642173963830608,
          0,
          3.209915703002990336815351572,
          0,
      },
      ATN2{
          0,
          -0.50000000000000000000000000000,
          1.29768738832023198239705918927,
          -1.94165942554412230232935244532,
          2.810188989257949038537763687,
          -4.2979081642651752142490066598,
          7.5692898193484870087981787234,
          -18.9375986071173705131546424941,
          0,
          14,
      }};
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, int N> struct get_coords {
  vect<T, dim> coord;
  T vol;
  get_coords(const Loop::GridDescBase &grid, const Loop::PointDesc &p) {
    vect<int, dim> i0, i1;
    vect<T, dim> x0, dx;
    for (int d = 0; d < dim; ++d) {
      i1[d] = mod_floor(p.I[d] - grid.nghostzones[d], N);
      i0[d] = p.I[d] - i1[d];
      x0[d] = grid.x0[d] + (grid.lbnd[d] + i0[d] - T(0.5)) * grid.dx[d];
      dx[d] = N * grid.dx[d];
      coord[d] = x0[d] + dx[d] / 2 + dx[d] / 2 * dg_basis<T, N>().coords[i1[d]];
    }
    vol = 1;
    for (int d = 0; d < dim; ++d)
      vol *= dx[d] / (2 * N) * dg_basis<T, N>().coords[i1[d]];
  }
};

template <int N, typename T>
vect<T, dim> get_derivs(const Loop::GridDescBase &grid,
                        const Loop::GF3D2<const T> &u,
                        const Loop::PointDesc &p) {
  vect<T, dim> deriv;
  vect<int, dim> i0, i1;
  vect<T, dim> dx;
  for (int d = 0; d < dim; ++d) {
    i1[d] = mod_floor(p.I[d] - grid.nghostzones[d], N);
    i0[d] = p.I[d] - i1[d];
    dx[d] = N * grid.dx[d];
    deriv[d] = 0;
    for (int j = 0; j < N + 2; ++j) {
      vect<int, dim> J = p.I;
      J[d] = i0[d] + j - 1;
      deriv[d] += dg_basis<T, N>().deriv[i1[d]][j] * u(J);
    }
    deriv[d] /= dx[d] / 2;
  }
  return deriv;
}

} // namespace DG

#endif // #ifndef DG_HXX
