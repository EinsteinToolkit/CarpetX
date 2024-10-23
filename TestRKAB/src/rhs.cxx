
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

namespace TestRKAB {
using namespace Arith;

// Finite difference helpers
enum class fd_dir : std::size_t { x = 0, y = 1, z = 2 };

template <fd_dir direction, typename T>
static inline auto fd_c_1_4(const Loop::PointDesc &p,
                            const Loop::GF3D2<T> &gf) noexcept -> T {
  constexpr auto d{static_cast<size_t>(direction)};
  const auto num{gf(p.I - 2 * p.DI[d]) - 8.0 * gf(p.I - 1 * p.DI[d]) +
                 8.0 * gf(p.I + 1 * p.DI[d]) - 1.0 * gf(p.I + 2 * p.DI[d])};
  const auto den{1.0 / (12.0 * p.DX[d])};
  return num * den;
}

extern "C" void TestRKAB_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestRKAB_RHS;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        phi_rhs(p.I) = Pi(p.I);
        Pi_rhs(p.I) = fd_c_1_4<fd_dir::x>(p, Dx) + fd_c_1_4<fd_dir::y>(p, Dy) +
                      fd_c_1_4<fd_dir::z>(p, Dz);
        Dx_rhs(p.I) = fd_c_1_4<fd_dir::x>(p, Pi);
        Dy_rhs(p.I) = fd_c_1_4<fd_dir::y>(p, Pi);
        Dz_rhs(p.I) = fd_c_1_4<fd_dir::z>(p, Pi);
      });
}

extern "C" void TestRKAB_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

} // namespace TestRKAB
