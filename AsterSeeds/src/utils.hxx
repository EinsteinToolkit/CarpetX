#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace AsterSeeds {
using namespace std;
using namespace Loop;

// Second-order average of cell-centered grid functions to edge center
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T calc_avg_c2e(const GF3D2<const T> &gf,
                                     const PointDesc &p, const int dir) {
  constexpr auto DI = PointDesc::DI;
  T gf_avg = 0.0;

  for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
    for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
      for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
        gf_avg += gf(p.I - DI[0] * di - DI[1] * dj - DI[2] * dk);
      }
    }
  }
  return gf_avg / 4.0;
}

} //namespace AsterSeeds
