#ifndef CARPETX_CARPETX_BOUNDARIES_HXX
#define CARPETX_CARPETX_BOUNDARIES_HXX

#include "driver.hxx"
#include "loop_device.hxx"

namespace CarpetX {

struct BoundaryCondition {

  const GHExt::PatchData::LevelData::GroupData &groupdata;
  const GHExt::PatchData &patchdata;
  const amrex::Geometry &geom;
  amrex::FArrayBox &dest;

  // Interior of the domain
  Arith::vect<int, dim> imin, imax;
  Arith::vect<CCTK_REAL, dim> xmin, xmax, dx;

  // Region to fill
  Arith::vect<int, dim> amin, amax;
  // Destination region
  Arith::vect<int, dim> dmin, dmax;

  Loop::GF3D2layout layout;
  CCTK_REAL *restrict destptr;

  BoundaryCondition(const GHExt::PatchData::LevelData::GroupData &groupdata,
                    const amrex::Box &box, amrex::FArrayBox &dest);

  BoundaryCondition(const BoundaryCondition &) = delete;
  BoundaryCondition(BoundaryCondition &&) = delete;
  BoundaryCondition &operator=(const BoundaryCondition &) = delete;
  BoundaryCondition &operator=(BoundaryCondition &&) = delete;

  void apply() const;

  template <int NI, int NJ, int NK> void apply_on_face() const;

  template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI>
  void apply_on_face_symbcx(const Arith::vect<int, dim> &bmin,
                            const Arith::vect<int, dim> &bmax) const;

  template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
            symmetry_t SCJ, boundary_t BCJ>
  void apply_on_face_symbcxy(const Arith::vect<int, dim> &bmin,
                             const Arith::vect<int, dim> &bmax) const;

  template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
            symmetry_t SCJ, boundary_t BCJ, symmetry_t SCK, boundary_t BCK>
  void apply_on_face_symbcxyz(const Arith::vect<int, dim> &bmin,
                              const Arith::vect<int, dim> &bmax) const;
};

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_BOUNDARIES_HXX
