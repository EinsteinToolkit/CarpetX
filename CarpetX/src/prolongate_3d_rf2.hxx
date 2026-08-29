#ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_HXX
#define CARPETX_CARPETX_PROLONGATE_3D_RF2_HXX

#include "driver.hxx"

#include <AMReX_Interpolater.H>

#include <array>
#include <cassert>
#include <iostream>
#include <map>

namespace CarpetX {

enum class centering_t : int { vertex = 0, cell = 1 };
std::ostream &operator<<(std::ostream &os, centering_t cent);
constexpr auto VC = centering_t::vertex;
constexpr auto CC = centering_t::cell;

enum class interpolation_t { poly, hermite, cons, eno, minmod };
std::ostream &operator<<(std::ostream &os, interpolation_t cent);
constexpr auto POLY = interpolation_t::poly;
constexpr auto HERMITE = interpolation_t::hermite;
constexpr auto CONS = interpolation_t::cons;
constexpr auto ENO = interpolation_t::eno;
constexpr auto MINMOD = interpolation_t::minmod;

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
                INTPI == ENO || INTPI == MINMOD);
  static_assert(INTPJ == POLY || INTPJ == HERMITE || INTPJ == CONS ||
                INTPJ == ENO || INTPJ == MINMOD);
  static_assert(INTPK == POLY || INTPK == HERMITE || INTPK == CONS ||
                INTPK == ENO || INTPK == MINMOD);

  // Order must be nonnegative
  static_assert(ORDERI >= 0);
  static_assert(ORDERJ >= 0);
  static_assert(ORDERK >= 0);

  // Minmod is always linear
  static_assert(INTPI == MINMOD ? ORDERI == 1 : true);
  static_assert(INTPJ == MINMOD ? ORDERJ == 1 : true);
  static_assert(INTPK == MINMOD ? ORDERK == 1 : true);

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

extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_poly_3d_rf2;

// Conservative interpolation

extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_cons_3d_rf2;

// DDF interpolation

// Prolongation operators for discrete differential forms:
// interpolating (non-conservative) for vertex centred directions,
// conservative (with one order lower) for cell centred directions.
extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_ddf_3d_rf2;

// Natural interpolation

// Interpolate (non-conservatively) for vertex centred directions,
// conservative for cell centred directions.
extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_natural_3d_rf2;

// ENO (tensor product) interpolation

extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_eno_3d_rf2;

// Minmod (tensor product) interpolation

extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_minmod_3d_rf2;

// Hermite interpolation

extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_hermite_3d_rf2;

// Interpolate polynomially in vertex centred directions and conserve
// with 3rd order accuracy and a linear fallback in cell centred
// directions

extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_poly_cons3lfb_3d_rf2;

// Interpolate polynomially in vertex centred directions and use ENO
// interpolation with 3rd order accuracy and a linear fallback in cell
// centred directions

extern const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_poly_eno3lfb_3d_rf2;

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_PROLONGATE_3D_RF2_HXX
