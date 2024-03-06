#ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_HXX
#define CARPETX_CARPETX_PROLONGATE_3D_RF2_HXX

#include "driver.hxx"

#include <AMReX_Interpolater.H>

#include <cassert>
#include <iostream>

namespace CarpetX {

enum class centering_t : int { vertex = 0, cell = 1 };
std::ostream &operator<<(std::ostream &os, centering_t cent);
constexpr auto VC = centering_t::vertex;
constexpr auto CC = centering_t::cell;

enum class interpolation_t { poly, hermite, cons, eno };
std::ostream &operator<<(std::ostream &os, interpolation_t cent);
constexpr auto POLY = interpolation_t::poly;
constexpr auto HERMITE = interpolation_t::hermite;
constexpr auto CONS = interpolation_t::cons;
constexpr auto ENO = interpolation_t::eno;

enum class fallback_t { none, linear };
std::ostream &operator<<(std::ostream &os, fallback_t fb);
constexpr auto FB_NONE = fallback_t::none;
constexpr auto FB_LINEAR = fallback_t::linear;

template <centering_t CENTI, centering_t CENTJ, centering_t CENTK,
          interpolation_t INTPI, interpolation_t INTPJ, interpolation_t INTPK,
          int ORDERI, int ORDERJ, int ORDERK, fallback_t FB>
class prolongate_3d_rf2 final : public amrex::Interpolater {

  // Centering must be vertex (0) or cell (1)
  static_assert(CENTI == VC || CENTI == CC);
  static_assert(CENTJ == VC || CENTJ == CC);
  static_assert(CENTK == VC || CENTK == CC);

  // Conservative must be one of the possible choices
  static_assert(INTPI == POLY || INTPI == HERMITE || INTPI == CONS ||
                INTPI == ENO);
  static_assert(INTPJ == POLY || INTPJ == HERMITE || INTPJ == CONS ||
                INTPJ == ENO);
  static_assert(INTPK == POLY || INTPK == HERMITE || INTPK == CONS ||
                INTPK == ENO);

  // Order must be nonnegative
  static_assert(ORDERI >= 0);
  static_assert(ORDERJ >= 0);
  static_assert(ORDERK >= 0);

  // Fallback must be one of the possible choices
  static_assert(FB == FB_NONE || FB == FB_LINEAR);

  static constexpr std::array<centering_t, dim> indextype() {
    return {CENTI, CENTJ, CENTK};
  }
  static constexpr std::array<interpolation_t, dim> interpolation() {
    return {INTPI, INTPJ, INTPK};
  }
  static constexpr std::array<int, dim> order() {
    return {ORDERI, ORDERJ, ORDERK};
  }
  static constexpr fallback_t fallback() { return FB; }

public:
  virtual ~prolongate_3d_rf2() override;

  virtual amrex::Box CoarseBox(const amrex::Box &fine, int ratio) override;
  virtual amrex::Box CoarseBox(const amrex::Box &fine,
                               const amrex::IntVect &ratio) override;

#ifndef AMREX_USE_GPU
private:
#endif
  void interp_per_var(const amrex::FArrayBox &crse, int crse_comp,
                      amrex::FArrayBox &fine, int fine_comp, int ncomp,
                      const amrex::Box &fine_region,
                      const amrex::IntVect &ratio,
                      const amrex::Geometry &crse_geom,
                      const amrex::Geometry &fine_geom,
                      amrex::Vector<amrex::BCRec> const &bcr, int actual_comp,
                      int actual_state, amrex::RunOn gpu_or_cpu);
  void interp_per_group(const amrex::FArrayBox &crse, int crse_comp,
                        amrex::FArrayBox &fine, int fine_comp, int ncomp,
                        const amrex::Box &fine_region,
                        const amrex::IntVect &ratio,
                        const amrex::Geometry &crse_geom,
                        const amrex::Geometry &fine_geom,
                        amrex::Vector<amrex::BCRec> const &bcr, int actual_comp,
                        int actual_state, amrex::RunOn gpu_or_cpu);

public:
  virtual void interp(const amrex::FArrayBox &crse, int crse_comp,
                      amrex::FArrayBox &fine, int fine_comp, int ncomp,
                      const amrex::Box &fine_region,
                      const amrex::IntVect &ratio,
                      const amrex::Geometry &crse_geom,
                      const amrex::Geometry &fine_geom,
                      amrex::Vector<amrex::BCRec> const &bcr, int actual_comp,
                      int actual_state, amrex::RunOn gpu_or_cpu) override;

  virtual void interp_face(const amrex::FArrayBox &crse, int crse_comp,
                           amrex::FArrayBox &fine, int fine_comp, int ncomp,
                           const amrex::Box &fine_region,
                           const amrex::IntVect &ratio,
                           const amrex::IArrayBox &solve_mask,
                           const amrex::Geometry &crse_geom,
                           const amrex::Geometry &fine_geom,
                           amrex::Vector<amrex::BCRec> const &bcr,
                           int actual_bccomp, amrex::RunOn gpu_or_cpu) override;
};

////////////////////////////////////////////////////////////////////////////////

// Polynomial (Lagrange) interpolation

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_poly_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_poly_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c000_o5;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c001_o5;
extern prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c010_o5;
extern prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c011_o5;
extern prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c100_o5;
extern prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c101_o5;
extern prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c110_o5;
extern prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_poly_3d_rf2_c111_o5;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c000_o7;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c001_o7;
extern prolongate_3d_rf2<VC, CC, VC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c010_o7;
extern prolongate_3d_rf2<VC, CC, CC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c011_o7;
extern prolongate_3d_rf2<CC, VC, VC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c100_o7;
extern prolongate_3d_rf2<CC, VC, CC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c101_o7;
extern prolongate_3d_rf2<CC, CC, VC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c110_o7;
extern prolongate_3d_rf2<CC, CC, CC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_poly_3d_rf2_c111_o7;

// Conservative interpolation

extern prolongate_3d_rf2<VC, VC, VC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c000_o0;
extern prolongate_3d_rf2<VC, VC, CC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c001_o0;
extern prolongate_3d_rf2<VC, CC, VC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c010_o0;
extern prolongate_3d_rf2<VC, CC, CC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c011_o0;
extern prolongate_3d_rf2<CC, VC, VC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c100_o0;
extern prolongate_3d_rf2<CC, VC, CC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c101_o0;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c110_o0;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_cons_3d_rf2_c111_o0;

extern prolongate_3d_rf2<VC, VC, VC, CONS, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_cons_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, CONS, CONS, CONS, 1, 1, 2, FB_NONE>
    prolongate_cons_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, CONS, CONS, CONS, 1, 2, 1, FB_NONE>
    prolongate_cons_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, CONS, CONS, CONS, 1, 2, 2, FB_NONE>
    prolongate_cons_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, CONS, CONS, CONS, 2, 1, 1, FB_NONE>
    prolongate_cons_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, CONS, CONS, CONS, 2, 1, 2, FB_NONE>
    prolongate_cons_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, CONS, 2, 2, 1, FB_NONE>
    prolongate_cons_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 2, 2, 2, FB_NONE>
    prolongate_cons_3d_rf2_c111_o1;

// DDF interpolation

// Prolongation operators for discrete differential forms:
// interpolating (non-conservative) for vertex centred directions,
// conservative (with one order lower) for cell centred directions.
extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 0, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 0, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 0, 1, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 0, 1, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 0, 0, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 3, 3, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 3, 2, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 3, 2, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 2, 3, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 2, 3, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 2, 2, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 2, 2, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o5;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 5, 5, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o5;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 5, 4, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o5;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 5, 4, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o5;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 4, 5, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o5;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 4, 5, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o5;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 4, 4, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o5;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 4, 4, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o5;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o7;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 7, 7, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o7;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 7, 6, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o7;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 7, 6, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o7;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 6, 7, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o7;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 6, 7, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o7;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 6, 6, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o7;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 6, 6, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o7;

// Natural interpolation

// Interpolate (non-conservatively) for vertex centred directions,
// conservative for cell centred directions.
extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_natural_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_natural_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c000_o5;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c001_o5;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c010_o5;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c011_o5;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c100_o5;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c101_o5;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c110_o5;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_natural_3d_rf2_c111_o5;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c000_o7;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c001_o7;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c010_o7;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c011_o7;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c100_o7;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c101_o7;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c110_o7;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 7, 7, 7, FB_NONE>
    prolongate_natural_3d_rf2_c111_o7;

// ENO (tensor product) interpolation

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_eno_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 1, 1, 0, FB_NONE>
    prolongate_eno_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 1, 0, 1, FB_NONE>
    prolongate_eno_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 1, 0, 0, FB_NONE>
    prolongate_eno_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 0, 1, 1, FB_NONE>
    prolongate_eno_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 0, 1, 0, FB_NONE>
    prolongate_eno_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 0, 0, 1, FB_NONE>
    prolongate_eno_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 0, 0, 0, FB_NONE>
    prolongate_eno_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_eno_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 3, 3, 2, FB_NONE>
    prolongate_eno_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 3, 2, 3, FB_NONE>
    prolongate_eno_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 3, 2, 2, FB_NONE>
    prolongate_eno_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 2, 3, 3, FB_NONE>
    prolongate_eno_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 2, 3, 2, FB_NONE>
    prolongate_eno_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 2, 2, 3, FB_NONE>
    prolongate_eno_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 2, 2, 2, FB_NONE>
    prolongate_eno_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_eno_3d_rf2_c000_o5;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 5, 5, 4, FB_NONE>
    prolongate_eno_3d_rf2_c001_o5;
extern prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 5, 4, 5, FB_NONE>
    prolongate_eno_3d_rf2_c010_o5;
extern prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 5, 4, 4, FB_NONE>
    prolongate_eno_3d_rf2_c011_o5;
extern prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 4, 5, 5, FB_NONE>
    prolongate_eno_3d_rf2_c100_o5;
extern prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 4, 5, 4, FB_NONE>
    prolongate_eno_3d_rf2_c101_o5;
extern prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 4, 4, 5, FB_NONE>
    prolongate_eno_3d_rf2_c110_o5;
extern prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 4, 4, 4, FB_NONE>
    prolongate_eno_3d_rf2_c111_o5;

// Hermite interpolation

extern prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 1, 1, 1,
                         FB_NONE>
    prolongate_hermite_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 3, 3, 3,
                         FB_NONE>
    prolongate_hermite_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 5, 5, 5,
                         FB_NONE>
    prolongate_hermite_3d_rf2_c000_o5;
extern prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o5;
extern prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o5;
extern prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o5;
extern prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o5;
extern prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o5;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o5;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o5;

// Interpolate polynomially in vertex centred directions and conserve
// with 3rd order accuracy and a linear fallback in cell centred
// directions

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c000_o1;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c001_o1;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c010_o1;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c011_o1;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c100_o1;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c101_o1;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c110_o1;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 1, 1, 1, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c111_o1;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c000_o3;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c001_o3;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c010_o3;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c011_o3;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c100_o3;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c101_o3;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c110_o3;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c111_o3;

extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c000_o5;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 5, 5, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c001_o5;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 5, 3, 5, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c010_o5;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 5, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c011_o5;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 3, 5, 5, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c100_o5;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 3, 5, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c101_o5;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 3, 3, 5, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c110_o5;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c111_o5;

#if 0
extern prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 7, 7, 7, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c000_o7;
extern prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 7, 7, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c001_o7;
extern prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 7, 3, 7, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c010_o7;
extern prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 7, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c011_o7;
extern prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 3, 7, 7, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c100_o7;
extern prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 3, 7, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c101_o7;
extern prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 3, 3, 7, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c110_o7;
extern prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 3, 3, 3, FB_LINEAR>
    prolongate_poly_cons3lfb_3d_rf2_c111_o7;
#endif

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_HXX
