#ifndef CARPETX_CARPETX_BOUNDARIES_HXX
#define CARPETX_CARPETX_BOUNDARIES_HXX

#include "driver.hxx"
#include "loop_device.hxx"

namespace CarpetX {

namespace {

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

template <typename F> struct task_t {
  F kernel;

  // Region for this task
  Arith::vect<int, dim> imin, imax;
  int cmin, cmax;

  // Convert the region to an AMReX box
  CCTK_DEVICE CCTK_HOST inline CCTK_ATTRIBUTE_ALWAYS_INLINE amrex::Box
  box() const noexcept {
    assert(all(imin < imax));
    assert(cmin < cmax);
    return amrex::Box(amrex::IntVect(imin[0], imin[1], imin[2]),
                      amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1));
  }
};

template <typename F>
task_t<F> make_task(const F &kernel, const Arith::vect<int, dim> &imin,
                    const Arith::vect<int, dim> &imax, const int cmin,
                    const int cmax) {
  return {kernel, imin, imax, cmin, cmax};
}

template <typename F> void run_task(const task_t<F> &task) {
  amrex::ParallelFor(
      task.box(),
      [kernel = task.kernel, cmin = task.cmin, cmax = task.cmax] CCTK_DEVICE(
          const int i, const int j, const int k) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const Arith::vect<int, dim> p{i, j, k};
        kernel(p, cmin, cmax);
      });
}

template <typename F> using tasks_t = amrex::Vector<task_t<F> >;

template <typename F> void run_tasks(const tasks_t<F> &tasks) {
#ifdef AMREX_USE_GPU
  amrex::ParallelFor(tasks, [=] CCTK_DEVICE(const int i, const int j,
                                            const int k, const task_t<F> &task)
                                CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                  const Arith::vect<int, dim> p{i, j, k};
                                  task.kernel(p, task.cmin, task.cmax);
                                });
#else
  for (const auto &task : tasks) {
    for (int comp = task.cmin; comp < task.cmax; ++comp) {
      amrex::ParallelFor(task.box(),
                         [kernel = task.kernel, comp] CCTK_DEVICE(
                             const int i, const int j, const int k)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               const Arith::vect<int, dim> p{i, j, k};
                               kernel(p, comp, comp + 1);
                             });
    }
  }
#endif
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

struct BoundaryKernel {
#ifdef CCTK_DEBUG
  // Region to fill
  Arith::vect<int, dim> amin, amax;
  // Destination region
  Arith::vect<int, dim> dmin, dmax;
#endif

  // Interior of the domain
  Arith::vect<CCTK_REAL, dim> xmin, dx;

  Loop::GF3D2layout layout;
  CCTK_REAL *restrict destptr;

  Arith::vect<int, dim> inormal;
  Arith::vect<symmetry_t, dim> symmetries;
  Arith::vect<boundary_t, dim> boundaries;

  static constexpr int maxncomps = 10;

  Arith::vect<CCTK_REAL, maxncomps> dirichlet_values;

  Arith::vect<int, dim> neumann_source;

  Arith::vect<int, dim> linear_extrapolation_source;

  Arith::vect<int, dim> robin_source;
  Arith::vect<CCTK_REAL, maxncomps> robin_values;

  Arith::vect<int, dim> reflection_offset;
  Arith::vect<CCTK_REAL, maxncomps> reflection_parities;

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE void
  operator()(const Arith::vect<int, dim> &dst, int cmin, int cmax) const;
};

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

  BoundaryKernel boundary_kernel(
      const Arith::vect<int, dim> inormal,
      const Arith::vect<symmetry_t, dim> symmetries,
      const Arith::vect<boundary_t, dim> boundaries,
      //
      const Arith::vect<CCTK_REAL, BoundaryKernel::maxncomps> &dirichlet_values,
      const Arith::vect<int, dim> &neumann_source,
      const Arith::vect<int, dim> &linear_extrapolation_source,
      const Arith::vect<int, dim> &robin_source,
      const Arith::vect<CCTK_REAL, BoundaryKernel::maxncomps> &robin_values,
      const Arith::vect<int, dim> &reflection_offset,
      const Arith::vect<CCTK_REAL, BoundaryKernel::maxncomps>
          &reflection_parities) const {
    return BoundaryKernel{
#ifdef CCTK_DEBUG
        amin, amax, dmin, dmax,
#endif
        xmin, dx, layout, destptr,
        //
        inormal, symmetries, boundaries,
        //
        dirichlet_values, neumann_source, linear_extrapolation_source,
        robin_source, robin_values, reflection_offset, reflection_parities};
  }

  mutable tasks_t<BoundaryKernel> tasks;

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
