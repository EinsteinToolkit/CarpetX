#include <loop_device.hxx>

#include <vect.hxx>
#include <sbp.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace SBPWaveToy {

namespace standing {

// Exact standing wave solution:
// u(t,x,y,z) = A cos(2\pi \omega t) cos(2\pi kx x) cos(2\pi ky y) cos(2\pi kz
// z) \omega = sqrt(kx^2 + ky^2 + kz^2)
template <typename T>
static inline auto uvw(T amplitude, T kx, T ky, T kz, T t, T x, T y, T z) {
  using std::sqrt, std::cos, std::acos;

  const auto pi = acos(-1.0);
  const auto omega = sqrt(kx * kx + ky * ky + kz * kz);

  const auto u = amplitude * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
                 cos(2 * pi * ky * y) * cos(2 * pi * kz * z);

  const auto v = -amplitude * (2 * pi * omega) * sin(2 * pi * omega * t) *
                 cos(2 * pi * kx * x) * cos(2 * pi * ky * y) *
                 cos(2 * pi * kz * z);

  const auto wx = -amplitude * (2 * pi * kx) * cos(2 * pi * omega * t) *
                  sin(2 * pi * kx * x) * cos(2 * pi * ky * y) *
                  cos(2 * pi * kz * z);

  const auto wy = -amplitude * (2 * pi * ky) * cos(2 * pi * omega * t) *
                  cos(2 * pi * kx * x) * sin(2 * pi * ky * y) *
                  cos(2 * pi * kz * z);

  const auto wz = -amplitude * (2 * pi * kz) * cos(2 * pi * omega * t) *
                  cos(2 * pi * kx * x) * cos(2 * pi * ky * y) *
                  sin(2 * pi * kz * z);

  return std::make_tuple(u, v, wx, wy, wz);
}

} // namespace standing

extern "C" void SBPWaveToy_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SBPWaveToy_Initial;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto [lu, lv, lwx, lwy, lwz] =
            standing::uvw(amplitude, standing_wave_kx, standing_wave_ky,
                          standing_wave_kz, cctk_time, p.x, p.y, p.z);

        u(p.I) = lu;
        v(p.I) = lv;
        wx(p.I) = lwx;
        wy(p.I) = lwy;
        wz(p.I) = lwz;
      });
}

extern "C" void SBPWaveToy_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SBPWaveToy_RHS;
  DECLARE_CCTK_PARAMETERS;

  // ---------------------------------------------------------------------------
  // How the SAT boundary penalties are obtained (energy method)
  // ---------------------------------------------------------------------------
  //
  // We evolve the first-order reduction of the scalar wave equation
  //
  //     du/dt   = v
  //     dv/dt   = dx(wx) + dy(wy) + dz(wz)   (= div w = Laplacian u)
  //     dw_i/dt = d_i v
  //
  // The variable u carries no spatial derivative, so it needs no boundary
  // treatment for stability; the hyperbolic part is the (v, w) subsystem.
  // Its continuous energy is
  //
  //     E = 1/2 \int ( v^2 + w_x^2 + w_y^2 + w_z^2 ) dV ,
  //     dE/dt = \oint v ( w . n ) dS                         (n = outward
  //     normal)
  //
  // so energy can only enter/leave through the faces. We reproduce this
  // estimate discretely and then add Simultaneous-Approximation-Term (SAT)
  // penalties at the boundary points to force dE/dt <= 0.
  //
  // The SBP operator factors as  D = H^{-1} Q  with a diagonal norm H and
  // Q + Q^T = B = diag(-1, 0, ..., 0, +1). With E = 1/2 (v^T H v + w^T H w) the
  // SBP property collapses the interior terms and leaves exactly the boundary
  // flux, per direction d:
  //
  //     dE_d/dt = s * v_b * w_{d,b} ,        s = +1 upper face, -1 lower face.
  //
  // Here s is precisely p.BI[d] (it is +1 at the last interior point, -1 at the
  // first, and 0 everywhere else), and the subscript b denotes the boundary
  // point. We cancel this flux with two penalty terms applied only at boundary
  // points. Writing the SAT-scaled coefficient
  //
  //     sigma = 1 / H_00 = 1 / (h_0 * dx) = inv_h0 / dx ,   h_0 = 17/48,
  //
  // (H_00 = h_0 * dx is the norm weight of the outermost boundary row, the only
  // diagonal entry that differs from the interior weight 1) the penalties are
  //
  //     v_rhs   -= sigma * (v - v*)
  //     w_d_rhs -= s * sigma * (v - v*)
  //
  // where v* is the target boundary value of v. Substituting these into the
  // energy balance, the factor 1/H_00 makes H * (H^{-1} e_b) = e_b, so the two
  // penalties contribute exactly  -v_b (v_b - v*)  and  -s v_b w_{d,b}  to
  // dE_d/dt. The cross term  -s v_b w_{d,b}  cancels the SBP flux  +s v_b
  // w_{d,b} identically, leaving
  //
  //     dE_d/dt = -v_b (v_b - v*) .
  //
  // For a Dirichlet wall the boundary value of u is held fixed, so v* =
  // du/dt|_b = 0 and dE_d/dt = -v_b^2 <= 0: the scheme is provably energy
  // dissipative, each face acting independently.
  //
  // "reflecting" and "dirichlet" share this exact (v, w) penalty (both freeze
  // the boundary in time, v* = 0). They differ only in how the *value* of u at
  // the boundary is set:
  //   * "reflecting" leaves u at whatever value it was initialized to (v = 0
  //     simply freezes it there).
  //   * "dirichlet" additionally relaxes u toward a chosen constant g =
  //     dirichlet_value with the matching penalty
  //
  //         u_rhs -= sigma * (u - g)        (boundary points only)
  //
  // This is the SAT for the trivial transport equation du/dt = v. Because u
  // does not appear in E, it cannot affect the (v, w) energy estimate above;
  // it only adds a dissipative pull (du/dt = v - sigma (u - g)) whose
  // equilibrium v = 0, u = g is exactly the fixed-value Dirichlet state.
  // ---------------------------------------------------------------------------

  // SBP operator
  const auto op = SBPOperators::ddst2007::get_op_42();

  // Boundary values
  const bool use_dirichlet = CCTK_EQUALS(boundary_condition, "dirichlet");
  const CCTK_REAL g = dirichlet_value;

  // 1 / h_0
  const CCTK_REAL inv_h0 =
      CCTK_REAL(1) / static_cast<CCTK_REAL>(op.boundary_h(0));

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // SBP derivative terms
        u_rhs(p.I) = v(p.I);

        v_rhs(p.I) =
            op.apply<0>(p, wx) + op.apply<1>(p, wy) + op.apply<2>(p, wz);

        wx_rhs(p.I) = op.apply<0>(p, v);
        wy_rhs(p.I) = op.apply<1>(p, v);
        wz_rhs(p.I) = op.apply<2>(p, v);

        // SAT penalties on physical boundary faces (p.BI[d] != 0). See the
        // derivation above.
        if (use_dirichlet) {
          for (int d = 0; d < Loop::dim; d++) {
            if (p.BI[d] != 0) {
              const auto sigma = inv_h0 / p.DX[d];
              const auto pen = v(p.I) * sigma;

              v_rhs(p.I) -= pen;

              if (d == 0) {
                wx_rhs(p.I) -= p.BI[d] * pen;
              } else if (d == 1) {
                wy_rhs(p.I) -= p.BI[d] * pen;
              } else {
                wz_rhs(p.I) -= p.BI[d] * pen;
              }

              // General fixed-value Dirichlet: pull u toward g. Applied once
              // per face the point belongs to; on an edge/corner the relaxation
              // simply acts in each contributing direction, all toward the same
              // g, which is consistent.
              u_rhs(p.I) -= sigma * (u(p.I) - g);
            }
          }
        }
      });
}

extern "C" void SBPWaveToy_calc_energy_density(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SBPWaveToy_calc_energy_density;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto lv = v(p.I);

        const auto lwx = wx(p.I);
        const auto lwy = wy(p.I);
        const auto lwz = wz(p.I);

        en_dens(p.I) = 0.5 * (lv * lv + lwx * lwx + lwy * lwy + lwz * lwz);
      });
}

extern "C" void SBPWaveToy_calc_evo_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SBPWaveToy_calc_evo_error;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::fabs;

        const auto [eu, ev, ewx, ewy, ewz] =
            standing::uvw(amplitude, standing_wave_kx, standing_wave_ky,
                          standing_wave_kz, cctk_time, p.x, p.y, p.z);

        u_error(p.I) = fabs(u(p.I) - eu);
        v_error(p.I) = fabs(v(p.I) - ev);
        wx_error(p.I) = fabs(wx(p.I) - ewx);
        wy_error(p.I) = fabs(wy(p.I) - ewy);
        wz_error(p.I) = fabs(wz(p.I) - ewz);
      });
}

extern "C" void SBPWaveToy_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SBPWaveToy_Sync;
  DECLARE_CCTK_PARAMETERS;
  // Empty body: the SYNC: state directive in schedule.ccl drives the ghost
  // zone exchange; this routine exists only as a schedule hook for it.
}

} // namespace SBPWaveToy
