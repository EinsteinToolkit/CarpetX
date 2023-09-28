#ifndef CARPETX_CARPETX_BOUNDARIES_HXX
#define CARPETX_CARPETX_BOUNDARIES_HXX

#include "driver.hxx"
#include "loop_device.hxx"

namespace CarpetX {

namespace boundaries_detail {

// TODO: Move these functions to loop.hxx

template <typename F>
CCTK_ATTRIBUTE_NOINLINE __attribute__((__flatten__, __hot__)) void
loop_region(const F &f, const Arith::vect<int, dim> &imin,
            const Arith::vect<int, dim> &imax) {
  assert(all(imin < imax));

  const amrex::Box box(amrex::IntVect(imin[0], imin[1], imin[2]),
                       amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1));
  amrex::ParallelFor(box, [=] CCTK_DEVICE(const int i, const int j, const int k)
                              CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                const Arith::vect<int, dim> p{i, j, k};
                                f(p);
                              });
}

} // namespace boundaries_detail
using namespace boundaries_detail;

////////////////////////////////////////////////////////////////////////////////

struct BoundaryCondition {
  const GHExt::PatchData::LevelData::GroupData &groupdata;
  const GHExt::PatchData &patchdata;
  const amrex::Geometry &geom;
  amrex::FArrayBox &dest;

  // Interior of the domain: Do not set any points in this region
  Arith::vect<int, dim> imin, imax;
  Arith::vect<CCTK_REAL, dim> xmin, xmax, dx;

  // Destination region
  Arith::vect<int, dim> dmin, dmax;

  Loop::GF3D2layout layout;
  CCTK_REAL *restrict destptr;

  BoundaryCondition(const GHExt::PatchData::LevelData::GroupData &groupdata,
                    amrex::FArrayBox &dest);

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
