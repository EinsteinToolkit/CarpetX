#ifndef LOOP_DEVICE_HXX
#define LOOP_DEVICE_HXX

#include "loop.hxx"

#include <cctk.h>

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#ifdef AMREX_USE_GPU
#include <AMReX_Gpu.H>
#endif

#include <cassert>

#ifndef HAVE_CAPABILITY_Loop
#error                                                                         \
    "Using #include <loop_device.hxx> requires the capability 'Loop' in the calling thorn. Add this to the thorn's 'configuration.ccl' file."
#endif

namespace Loop {

// The struct GridDescBaseDevice is very similar to GridDescBase, except that it
// requires AMReX include files, and thus "REQUIRES CarpetX" in the using
// thorns' configuration.ccl file.
struct GridDescBaseDevice : GridDescBase {

protected:
  GridDescBaseDevice() : GridDescBase() {}

public:
  GridDescBaseDevice(const GridDescBaseDevice &) = default;
  GridDescBaseDevice(GridDescBaseDevice &&) = default;
  GridDescBaseDevice &operator=(const GridDescBaseDevice &) = default;
  GridDescBaseDevice &operator=(GridDescBaseDevice &&) = default;

  GridDescBaseDevice(const cGH *cctkGH) : GridDescBase(cctkGH) {}

  // Loop over a given box
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  void loop_box_device(const F &f, const vect<int, dim> &restrict bnd_min,
                       const vect<int, dim> &restrict bnd_max,
                       const vect<int, dim> &restrict imin,
                       const vect<int, dim> &restrict imax) const {
#ifndef AMREX_USE_GPU

    return this->template loop_box<CI, CJ, CK, VS>(f, bnd_min, bnd_max, imin,
                                                   imax);

#else

    static_assert(CI == 0 || CI == 1);
    static_assert(CJ == 0 || CJ == 1);
    static_assert(CK == 0 || CK == 1);
    static_assert(VS > 0);

    if (any(imin >= imax))
      return;

    // Run on GPU
    static_assert(VS == 1, "Only vector size 1 is supported on GPUs");

    // For some reason, the arguments imin and imax cannot be captured
    // correctly in CUDA, but copies of them can
    const auto bnd_min1 = bnd_min;
    const auto bnd_max1 = bnd_max;
    const auto imin1 = imin;
    const auto imax1 = imax;

    // Convert to AMReX box
    const amrex::Box box(
        amrex::IntVect(imin[0], imin[1], imin[2]),
        amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1),
        amrex::IntVect(CI ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CJ ? amrex::IndexType::CELL : amrex::IndexType::NODE,
                       CK ? amrex::IndexType::CELL : amrex::IndexType::NODE));

    amrex::ParallelFor(
        box, [=, *this] CCTK_DEVICE(const int i, const int j,
                                    const int k) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vect<int, dim> I = {i, j, k};
          const vect<int, dim> NI =
              vect<int, dim>(I > bnd_max1 - 1) - vect<int, dim>(I < bnd_min1);
          const vect<int, dim> I0 =
              if_else(NI == 0, 0, if_else(NI < 0, bnd_min1, bnd_max1 - 1));
          const vect<int, dim> BI =
              vect<int, dim>(I == bnd_max1 - 1) - vect<int, dim>(I == bnd_min1);
          const PointDesc p =
              point_desc({CI, CJ, CK}, I, NI, I0, BI, imin1[0], imax1[0]);
          f(p);
        });

#endif

#ifdef AMREX_USE_GPU
    static bool have_gpu_sync_after_every_kernel = false;
    static bool gpu_sync_after_every_kernel;
    if (!have_gpu_sync_after_every_kernel) {
      int type;
      const void *const gpu_sync_after_every_kernel_ptr =
          CCTK_ParameterGet("gpu_sync_after_every_kernel", "CarpetX", &type);
      assert(gpu_sync_after_every_kernel_ptr);
      assert(type == PARAMETER_BOOLEAN);
      gpu_sync_after_every_kernel =
          *static_cast<const CCTK_INT *>(gpu_sync_after_every_kernel_ptr);
    }
    if (gpu_sync_after_every_kernel) {
      amrex::Gpu::synchronize();
      AMREX_GPU_ERROR_CHECK();
    }
#endif
  }

  // Loop over all points
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_all_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    vect<int, dim> bnd_min, bnd_max;
    boundary_box<CI, CJ, CK>(group_nghostzones, bnd_min, bnd_max);
    vect<int, dim> imin, imax;
    box_all<CI, CJ, CK>(group_nghostzones, imin, imax);
    loop_box_device<CI, CJ, CK, VS>(f, bnd_min, bnd_max, imin, imax);
  }

  // Loop over all interior points
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_int_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    vect<int, dim> bnd_min, bnd_max;
    boundary_box<CI, CJ, CK>(group_nghostzones, bnd_min, bnd_max);
    vect<int, dim> imin, imax;
    box_int<CI, CJ, CK>(group_nghostzones, imin, imax);
    loop_box_device<CI, CJ, CK, VS>(f, bnd_min, bnd_max, imin, imax);
  }

  // Loop over a part of the domain. Loop over the interior first,
  // then faces, then edges, then corners.
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_there_device(const vect<int, dim> &group_nghostzones,
                    const vect<vect<vect<bool, dim>, dim>, dim> &there,
                    const F &f) const {
    constexpr vect<int, dim> offset{CI, CJ, CK};

    vect<int, dim> bnd_min, bnd_max;
    boundary_box<CI, CJ, CK>(group_nghostzones, bnd_min, bnd_max);
    vect<int, dim> all_min, all_max, int_min, int_max;
    domain_boxes<CI, CJ, CK>(group_nghostzones, all_min, all_max, int_min,
                             int_max);

    for (int rank = dim; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if (there[ni + 1][nj + 1][nk + 1]) {

                const vect<int, dim> inormal{ni, nj, nk};

                vect<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = all_min[d];
                    imax[d] = int_min[d];
                    break;
                  case 0: // interior
                    imin[d] = int_min[d];
                    imax[d] = int_max[d];
                    break;
                  case +1: // upper boundary
                    imin[d] = int_max[d];
                    imax[d] = all_max[d];
                    break;
                  default:
                    assert(0);
                  }

                  using std::max, std::min;
                  imin[d] = max(tmin[d], imin[d]);
                  imax[d] = min(tmax[d], imax[d]);
                }

                loop_box_device<CI, CJ, CK, VS>(f, bnd_min, bnd_max, imin,
                                                imax);
              }
            } // if rank
          }
        }
      }

    } // for rank
  }

  // Loop over all outer boundary points. This excludes ghost faces, but
  // includes ghost edges/corners on non-ghost faces. Loop over faces first,
  // then edges, then corners.
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_bnd_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    vect<int, dim> bnd_min, bnd_max;
    boundary_box<CI, CJ, CK>(group_nghostzones, bnd_min, bnd_max);
    vect<int, dim> all_min, all_max, int_min, int_max;
    domain_boxes<CI, CJ, CK>(group_nghostzones, all_min, all_max, int_min,
                             int_max);

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni != 0 && bbox[ni == -1 ? 0 : 1][0]) ||
                  (nj != 0 && bbox[nj == -1 ? 0 : 1][1]) ||
                  (nk != 0 && bbox[nk == -1 ? 0 : 1][2])) {

                const vect<int, dim> inormal{ni, nj, nk};

                vect<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = all_min[d];
                    imax[d] = int_min[d];
                    break;
                  case 0: // interior
                    imin[d] = int_min[d];
                    imax[d] = int_max[d];
                    break;
                  case +1: // upper boundary
                    imin[d] = int_max[d];
                    imax[d] = all_max[d];
                    break;
                  default:
                    assert(0);
                  }

                  using std::min, std::max;
                  imin[d] = max(tmin[d], imin[d]);
                  imax[d] = min(tmax[d], imax[d]);
                }

                loop_box_device<CI, CJ, CK, VS>(f, bnd_min, bnd_max, imin,
                                                imax);
              }
            } // if rank
          }
        }
      }

    } // for rank
  }

#if 0
  // Loop over all outer ghost points. This includes ghost edges/corners on
  // non-ghost faces. Loop over faces first, then edges, then corners.
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_ghosts_inclusive_device(const vect<int, dim> &group_nghostzones,
                               const F &f) const {
    constexpr vect<int, dim> offset{CI, CJ, CK};

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni != 0 && !bbox[ni == -1 ? 0 : 1][0]) ||
                  (nj != 0 && !bbox[nj == -1 ? 0 : 1][1]) ||
                  (nk != 0 && !bbox[nk == -1 ? 0 : 1][2])) {

                const vect<int, dim> inormal{ni, nj, nk};

                vect<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  const int ghost_offset =
                      nghostzones[d] - group_nghostzones[d];
                  const int begin_bnd = ghost_offset;
                  const int begin_int = nghostzones[d];
                  const int end_int = lsh[d] - offset[d] - nghostzones[d];
                  const int end_bnd = lsh[d] - offset[d] - ghost_offset;
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = begin_bnd;
                    imax[d] = begin_int;
                    break;
                  case 0: // interior
                    imin[d] = begin_int;
                    imax[d] = end_int;
                    break;
                  case +1: // upper boundary
                    imin[d] = end_int;
                    imax[d] = end_bnd;
                    break;
                  default:
                    assert(0);
                  }

                  imin[d] = std::max(tmin[d], imin[d]);
                  imax[d] = std::min(tmax[d], imax[d]);
                }

                loop_box_boundary_device<CI, CJ, CK, VS>(f, imin, imax,
                                                         inormal);
              }
            } // if rank
          }
        }
      }
    } // for rank
  }
#endif

  // Loop over all outer ghost points. This excludes ghost edges/corners on
  // non-ghost faces. Loop over faces first, then edges, then corners.
  template <int CI, int CJ, int CK, int VS = 1, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_ghosts_device(const vect<int, dim> &group_nghostzones,
                     const F &f) const {
    vect<int, dim> bnd_min, bnd_max;
    boundary_box<CI, CJ, CK>(group_nghostzones, bnd_min, bnd_max);
    vect<int, dim> all_min, all_max, int_min, int_max;
    domain_boxes<CI, CJ, CK>(group_nghostzones, all_min, all_max, int_min,
                             int_max);

    for (int rank = dim - 1; rank >= 0; --rank) {

      for (int nk = -1; nk <= +1; ++nk) {
        for (int nj = -1; nj <= +1; ++nj) {
          for (int ni = -1; ni <= +1; ++ni) {
            if ((ni == 0) + (nj == 0) + (nk == 0) == rank) {

              if ((ni == 0 || !bbox[ni == -1 ? 0 : 1][0]) &&
                  (nj == 0 || !bbox[nj == -1 ? 0 : 1][1]) &&
                  (nk == 0 || !bbox[nk == -1 ? 0 : 1][2])) {

                const vect<int, dim> inormal{ni, nj, nk};

                vect<int, dim> imin, imax;
                for (int d = 0; d < dim; ++d) {
                  switch (inormal[d]) {
                  case -1: // lower boundary
                    imin[d] = all_min[d];
                    imax[d] = int_min[d];
                    break;
                  case 0: // interior
                    imin[d] = int_min[d];
                    imax[d] = int_max[d];
                    break;
                  case +1: // upper boundary
                    imin[d] = int_max[d];
                    imax[d] = all_max[d];
                    break;
                  default:
                    assert(0);
                  }

                  using std::min, std::max;
                  imin[d] = max(tmin[d], imin[d]);
                  imax[d] = min(tmax[d], imax[d]);
                }

                loop_box_device<CI, CJ, CK, VS>(f, bnd_min, bnd_max, imin,
                                                imax);
              }
            } // if rank
          }
        }
      }
    } // for rank
  }

  template <int CI, int CJ, int CK, where_t where, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
      std::enable_if_t<(where == where_t::everywhere), void>
      loop_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    loop_all_device<CI, CJ, CK>(group_nghostzones, f);
  }
  template <int CI, int CJ, int CK, where_t where, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
      std::enable_if_t<(where == where_t::interior), void>
      loop_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    loop_int_device<CI, CJ, CK>(group_nghostzones, f);
  }
  template <int CI, int CJ, int CK, where_t where, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
      std::enable_if_t<(where == where_t::boundary), void>
      loop_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    loop_bnd_device<CI, CJ, CK>(group_nghostzones, f);
  }
#if 0
  template <int CI, int CJ, int CK, where_t where, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
      std::enable_if_t<(where == where_t::ghosts_inclusive), void>
      loop_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    loop_ghosts_inclusive_device<CI, CJ, CK>(group_nghostzones, f);
  }
#endif
  template <int CI, int CJ, int CK, where_t where, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
      std::enable_if_t<(where == where_t::ghosts), void>
      loop_device(const vect<int, dim> &group_nghostzones, const F &f) const {
    loop_ghosts_device<CI, CJ, CK>(group_nghostzones, f);
  }

  template <typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_there_device_idx(const vect<int, dim> &indextype,
                        const vect<vect<vect<bool, dim>, dim>, dim> &there,
                        const vect<int, dim> &group_nghostzones,
                        const F &f) const {
    switch (indextype[0] + 2 * indextype[1] + 4 * indextype[2]) {
    case 0b000:
      return noinline([&] {
        return loop_there_device<0, 0, 0>(group_nghostzones, there, f);
      });
    case 0b001:
      return noinline([&] {
        return loop_there_device<1, 0, 0>(group_nghostzones, there, f);
      });
    case 0b010:
      return noinline([&] {
        return loop_there_device<0, 1, 0>(group_nghostzones, there, f);
      });
    case 0b011:
      return noinline([&] {
        return loop_there_device<1, 1, 0>(group_nghostzones, there, f);
      });
    case 0b100:
      return noinline([&] {
        return loop_there_device<0, 0, 1>(group_nghostzones, there, f);
      });
    case 0b101:
      return noinline([&] {
        return loop_there_device<1, 0, 1>(group_nghostzones, there, f);
      });
    case 0b110:
      return noinline([&] {
        return loop_there_device<0, 1, 1>(group_nghostzones, there, f);
      });
    case 0b111:
      return noinline([&] {
        return loop_there_device<1, 1, 1>(group_nghostzones, there, f);
      });
    default:
      assert(0);
    }
  }

  template <where_t where, typename F>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
  loop_device_idx(const vect<int, dim> &indextype,
                  const vect<int, dim> &group_nghostzones, const F &f) const {
    switch (indextype[0] + 2 * indextype[1] + 4 * indextype[2]) {
    case 0b000:
      return noinline(
          [&] { return loop_device<0, 0, 0, where>(group_nghostzones, f); });
    case 0b001:
      return noinline(
          [&] { return loop_device<1, 0, 0, where>(group_nghostzones, f); });
    case 0b010:
      return noinline(
          [&] { return loop_device<0, 1, 0, where>(group_nghostzones, f); });
    case 0b011:
      return noinline(
          [&] { return loop_device<1, 1, 0, where>(group_nghostzones, f); });
    case 0b100:
      return noinline(
          [&] { return loop_device<0, 0, 1, where>(group_nghostzones, f); });
    case 0b101:
      return noinline(
          [&] { return loop_device<1, 0, 1, where>(group_nghostzones, f); });
    case 0b110:
      return noinline(
          [&] { return loop_device<0, 1, 1, where>(group_nghostzones, f); });
    case 0b111:
      return noinline(
          [&] { return loop_device<1, 1, 1, where>(group_nghostzones, f); });
    default:
      assert(0);
    }
  }
};

template <where_t where, typename F>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
loop_device_idx(const cGH *cctkGH, const vect<int, dim> &indextype,
                const vect<int, dim> &nghostzones, const F &f) {
  GridDescBaseDevice(cctkGH).loop_device_idx<where>(indextype, nghostzones, f);
}

} // namespace Loop

#undef CCTK_CENTERING_GRID
#define CCTK_CENTERING_GRID const Loop::GridDescBaseDevice cctk_grid(cctkGH)

#endif // #ifndef LOOP_DEVICE_HXX
