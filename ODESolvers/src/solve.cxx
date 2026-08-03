// TODO: Don't include files from other thorns; create a proper interface
#include "../../CarpetX/src/driver.hxx"
#include "../../CarpetX/src/schedule.hxx"
#include "../../CarpetX/src/timer.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>
#include <util_Table.h>

#include <div.hxx>
#include "cactus_explicit_rk_operations.hxx"
#include "explicit_rk.hxx"
#include "explicit_rk_dense_provider.hxx"
#include "subcycling_group_schema_builder.hxx"
#include "subcycling_ode_provider_registry.hxx"
#include "subcycling_transaction_bridge.hxx"
#include <subcycling_dense_output.hxx>
#include <subcycling_scratch_state_transaction.hxx>
#include <subcycling_step_context.hxx>

#include <AMReX_MultiFab.H>

#if defined _OPENMP || defined __HIPCC__
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace ODESolvers {

////////////////////////////////////////////////////////////////////////////////

// Taken from <https://en.cppreference.com/w/cpp/experimental/make_array>
namespace details {
template <class> struct is_ref_wrapper : std::false_type {};
template <class T>
struct is_ref_wrapper<std::reference_wrapper<T> > : std::true_type {};

template <class T>
using not_ref_wrapper = std::negation<is_ref_wrapper<std::decay_t<T> > >;

template <class D, class...> struct return_type_helper {
  using type = D;
};
template <class... Types>
struct return_type_helper<void, Types...> : std::common_type<Types...> {
  static_assert(std::conjunction_v<not_ref_wrapper<Types>...>,
                "Types cannot contain reference_wrappers when D is void");
};

template <class D, class... Types>
using return_type = std::array<typename return_type_helper<D, Types...>::type,
                               sizeof...(Types)>;
} // namespace details

template <class D = void, class... Types>
constexpr details::return_type<D, Types...> make_array(Types &&...t) {
  return {std::forward<Types>(t)...};
}

////////////////////////////////////////////////////////////////////////////////

// A state vector component, with mfabs for each level, group, and variable
struct statecomp_t {
  using scalar_type = CCTK_REAL;

  statecomp_t() = default;

  statecomp_t(statecomp_t &&) = default;
  statecomp_t &operator=(statecomp_t &&) = default;

  // Don't allow copies because we might own stuff
  statecomp_t(const statecomp_t &) = delete;
  statecomp_t &operator=(const statecomp_t &) = delete;

  std::vector<CarpetX::GHExt::PatchData::LevelData::GroupData *> groupdatas;
  std::vector<amrex::MultiFab *> mfabs;
  std::optional<TransactionStateBackend> scratch_backend;

  static statecomp_t from_scratch(TransactionStateBackend backend) {
    statecomp_t result;
    result.scratch_backend.emplace(std::move(backend));
    return result;
  }

  bool is_scratch() const noexcept { return scratch_backend.has_value(); }
  TransactionStateBackend &scratch() {
    if (!scratch_backend)
      throw std::logic_error("state component is not scratch-backed");
    return *scratch_backend;
  }
  const TransactionStateBackend &scratch() const {
    if (!scratch_backend)
      throw std::logic_error("state component is not scratch-backed");
    return *scratch_backend;
  }
  void replace_scratch(TransactionStateBackend backend) {
    if (!is_scratch() || &scratch().transaction() != &backend.transaction())
      throw std::invalid_argument(
          "scratch state replacement changes backend or owner");
    scratch_backend.emplace(std::move(backend));
  }

  static void init_tmp_mfabs();
  static void free_tmp_mfabs();

  void set_valid(const CarpetX::valid_t valid) const;
  template <std::size_t N>
  static void combine_valids(const statecomp_t &dst, const CCTK_REAL scale,
                             const std::array<CCTK_REAL, N> &factors,
                             const std::array<const statecomp_t *, N> &srcs,
                             const CarpetX::valid_t where);
  void check_valid(const CarpetX::valid_t required,
                   const std::function<std::string()> &why) const;
  void check_valid(const CarpetX::valid_t required,
                   const std::string &why) const {
    check_valid(required, [=]() { return why; });
  }

  statecomp_t copy(const CarpetX::valid_t where) const;
  statecomp_t snapshot_state() const {
    return copy(CarpetX::make_valid_all());
  }
  statecomp_t snapshot_rhs() const {
    return copy(CarpetX::make_valid_int());
  }

  template <std::size_t N>
  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const std::array<CCTK_REAL, N> &factors,
                      const std::array<const statecomp_t *, N> &srcs,
                      const CarpetX::valid_t where);
  template <std::size_t N>
  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const std::array<CCTK_REAL, N> &factors,
                      const std::array<statecomp_t *, N> &srcs,
                      const CarpetX::valid_t where) {
    std::array<const statecomp_t *, N> srcs1;
    for (std::size_t n = 0; n < N; ++n)
      srcs1[n] = srcs[n];
    lincomb(dst, scale, factors, srcs1, where);
  }

  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const std::vector<CCTK_REAL> &factors,
                      const std::vector<const statecomp_t *> &srcs,
                      const CarpetX::valid_t where);
  static void lincomb(
      const statecomp_t &dst, CCTK_REAL scale,
      LinearCombinationView<CCTK_REAL, statecomp_t> combination,
      const CarpetX::valid_t where);
};

template <std::size_t N> using reals = std::array<CCTK_REAL, N>;
template <std::size_t N> using states = std::array<const statecomp_t *, N>;

////////////////////////////////////////////////////////////////////////////////

// Initialize the temporary mfab mechanism
void statecomp_t::init_tmp_mfabs() {
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      if (groupdataptr == nullptr)
        continue;
      const auto &groupdata = *groupdataptr;
      groupdata.init_tmp_mfabs();
    }
  });
}

// Free all temporary mfabs that we might have allocated
void statecomp_t::free_tmp_mfabs() {
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      if (groupdataptr == nullptr)
        continue;
      const auto &groupdata = *groupdataptr;
      groupdata.free_tmp_mfabs();
    }
  });
}

// State that the state vector has valid data in the interior
void statecomp_t::set_valid(const CarpetX::valid_t valid) const {
  if (is_scratch()) {
    auto &backend = const_cast<TransactionStateBackend &>(scratch());
    backend.set_valid(CarpetX::ScratchStateRegion::interior,
                      valid.valid_int);
    backend.set_valid(CarpetX::ScratchStateRegion::outer,
                      valid.valid_outer);
    backend.set_valid(CarpetX::ScratchStateRegion::ghosts,
                      valid.valid_ghosts);
    return;
  }
  for (auto groupdata : groupdatas) {
    for (int vi = 0; vi < groupdata->numvars; ++vi) {
      const int tl = 0;
      groupdata->valid.at(tl).at(vi).set_int(valid.valid_int, [=]() {
        std::ostringstream buf;
        buf << "ODESolvers after lincomb: Mark interior as "
            << (valid.valid_int ? "valid" : "invalid");
        return buf.str();
      });
      groupdata->valid.at(tl).at(vi).set_outer(valid.valid_outer, [=]() {
        std::ostringstream buf;
        buf << "ODESolvers after lincomb: Mark outer boundary as "
            << (valid.valid_outer ? "valid" : "invalid");
        return buf.str();
      });
      groupdata->valid.at(tl).at(vi).set_ghosts(valid.valid_ghosts, [=]() {
        std::ostringstream buf;
        buf << "ODESolvers after lincomb: Mark ghosts as "
            << (valid.valid_int ? "valid" : "invalid");
        return buf.str();
      });
      // TODO: Parallelize over patches, levels, group, variables, and
      // timelevels
      const CarpetX::active_levels_t active_levels(
          groupdata->level, groupdata->level + 1, groupdata->patch,
          groupdata->patch + 1);
      CarpetX::poison_invalid_gf(active_levels, groupdata->groupindex, vi, tl);
    }
  }
}

// Combine validity information from several sources
template <std::size_t N>
void statecomp_t::combine_valids(const statecomp_t &dst, const CCTK_REAL scale,
                                 const std::array<CCTK_REAL, N> &factors,
                                 const std::array<const statecomp_t *, N> &srcs,
                                 const CarpetX::valid_t where) {
  const int ngroups = dst.groupdatas.size();
  for (const auto &src : srcs)
    assert(int(src->groupdatas.size()) == ngroups);
  for (int group = 0; group < ngroups; ++group) {
    const auto &dstgroup = dst.groupdatas.at(group);
    const int nvars = dstgroup->numvars;
    for (const auto &src : srcs) {
      const auto &srcgroup = src->groupdatas.at(group);
      assert(srcgroup->numvars == nvars);
    }
  }

  for (int group = 0; group < ngroups; ++group) {
    const auto &dstgroup = dst.groupdatas.at(group);
    const int nvars = dstgroup->numvars;
    const int tl = 0;
    for (int vi = 0; vi < nvars; ++vi) {
      CarpetX::valid_t valid = where;
      bool did_set_valid = false;
      if (scale != 0) {
        valid &= dstgroup->valid.at(tl).at(vi).get();
        did_set_valid = true;
      }
      for (std::size_t m = 0; m < srcs.size(); ++m) {
        if (factors.at(m) != 0) {
          const auto &src = srcs.at(m);
          const auto &srcgroup = src->groupdatas.at(group);
          valid &= srcgroup->valid.at(tl).at(vi).get();
          did_set_valid = true;
        }
      }
      if (!did_set_valid)
        valid = CarpetX::valid_t(false);
      dstgroup->valid.at(tl).at(vi) = CarpetX::why_valid_t(
          valid, []() { return "Set from RHS in ODESolvers"; });
    }
  }
}

// Ensure a state vector has valid data in the interior
void statecomp_t::check_valid(const CarpetX::valid_t required,
                              const std::function<std::string()> &why) const {
  if (is_scratch()) {
    const bool valid =
        (!required.valid_int ||
         scratch().valid(CarpetX::ScratchStateRegion::interior)) &&
        (!required.valid_outer ||
         scratch().valid(CarpetX::ScratchStateRegion::outer)) &&
        (!required.valid_ghosts ||
         scratch().valid(CarpetX::ScratchStateRegion::ghosts));
    if (!valid)
      throw std::runtime_error(why());
    return;
  }
  for (const auto groupdata : groupdatas) {
    for (int vi = 0; vi < groupdata->numvars; ++vi) {
      const int tl = 0;
      CarpetX::error_if_invalid(*groupdata, vi, tl, required, why);
      // TODO: Parallelize over pathces, levels, group, variables, and
      // timelevels
      const CarpetX::active_levels_t active_levels(
          groupdata->level, groupdata->level + 1, groupdata->patch,
          groupdata->patch + 1);
      CarpetX::check_valid_gf(active_levels, groupdata->groupindex, vi, tl,
                              CarpetX::nan_handling_t::forbid_nans, why);
    }
  }
}

// Copy state vector into newly allocated memory
statecomp_t statecomp_t::copy(const CarpetX::valid_t where) const {
  if (is_scratch()) {
    auto result = statecomp_t::from_scratch(scratch().clone());
    if (!where.valid_int)
      result.scratch().set_valid(CarpetX::ScratchStateRegion::interior, false);
    if (!where.valid_outer)
      result.scratch().set_valid(CarpetX::ScratchStateRegion::outer, false);
    if (!where.valid_ghosts)
      result.scratch().set_valid(CarpetX::ScratchStateRegion::ghosts, false);
    return result;
  }
  const std::size_t size = mfabs.size();
  statecomp_t result;
  result.groupdatas.reserve(size);
  result.mfabs.reserve(size);
  for (std::size_t n = 0; n < size; ++n) {
    const auto groupdata = groupdatas.at(n);
    // This global nan-check doesn't work since we don't care about the
    // boundaries
    // #ifdef CCTK_DEBUG
    //     const auto &x = mfabs.at(n);
    //     if (x->contains_nan())
    //       CCTK_VERROR("statecomp_t::copy.x: Group %s contains nans",
    //                   groupdata->groupname.c_str());
    // #endif
    auto y = groupdata->alloc_tmp_mfab();
    result.groupdatas.push_back(groupdata);
    result.mfabs.push_back(y);
  }
  lincomb(result, 0, make_array(CCTK_REAL(1)), make_array(this), where);
  // This global nan-check doesn't work since we don't care about the boundaries
  // #ifdef CCTK_DEBUG
  //   for (std::size_t n = 0; n < size; ++n) {
  //     const auto groupdata = result.groupdatas.at(n);
  //     const auto &y = result.mfabs.at(n);
  //     if (y->contains_nan())
  //       CCTK_VERROR("statecomp_t::copy.y: Group %s contains nans",
  //                   groupdata->groupname.c_str());
  //   }
  // #endif
  return result;
}

template <std::size_t N>
void statecomp_t::lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                          const std::array<CCTK_REAL, N> &factors,
                          const std::array<const statecomp_t *, N> &srcs,
                          const CarpetX::valid_t where) {
  if (dst.is_scratch()) {
    if (!where.valid_int || where.valid_outer || where.valid_ghosts)
      throw std::invalid_argument(
          "scratch RK arithmetic supports interior validity only");
    std::vector<TransactionStateBackend::LinearTerm> terms;
    terms.reserve(N);
    for (std::size_t n = 0; n < N; ++n) {
      if (srcs[n] == nullptr || !srcs[n]->is_scratch())
        throw std::invalid_argument(
            "scratch RK arithmetic cannot mix live state components");
      terms.push_back({factors[n], &srcs[n]->scratch()});
    }
    const_cast<TransactionStateBackend &>(dst.scratch())
        .linear_combination(scale, terms);
    return;
  }
  for (const auto *const source : srcs)
    if (source == nullptr || source->is_scratch())
      throw std::invalid_argument(
          "live RK arithmetic cannot mix scratch state components");
  const std::size_t size = dst.mfabs.size();
  for (std::size_t n = 0; n < N; ++n)
    assert(srcs[n]->mfabs.size() == size);
  for (std::size_t m = 0; m < size; ++m) {
    const auto ncomp = dst.mfabs.at(m)->nComp();
    const auto ngrowvect = dst.mfabs.at(m)->nGrowVect();
    for (std::size_t n = 0; n < N; ++n) {
      assert(srcs[n]->mfabs.at(m)->nComp() == ncomp);
      assert(srcs[n]->mfabs.at(m)->nGrowVect() == ngrowvect);
    }
  }

  using std::isfinite;
  assert(isfinite(scale));
  const bool read_dst = scale != 0;
  for (std::size_t n = 0; n < N; ++n)
    assert(isfinite(factors[n]));

  statecomp_t::combine_valids(dst, scale, factors, srcs, where);

#ifndef AMREX_USE_GPU
  std::vector<std::function<void()> > tasks;
#endif

  // TODO: Poison ghosts/boundaries

  for (std::size_t m = 0; m < size; ++m) {
    const std::ptrdiff_t ncomps = dst.mfabs.at(m)->nComp();
    const auto mfitinfo = amrex::MFItInfo().DisableDeviceSync();
    for (amrex::MFIter mfi(*dst.mfabs.at(m), mfitinfo); mfi.isValid(); ++mfi) {
      const amrex::Array4<CCTK_REAL> dstvar = dst.mfabs.at(m)->array(mfi);
      std::array<amrex::Array4<const CCTK_REAL>, N> srcvars;
      for (std::size_t n = 0; n < N; ++n)
        srcvars[n] = srcs[n]->mfabs.at(m)->const_array(mfi);
      for (std::size_t n = 0; n < N; ++n) {
        assert(srcvars[n].jstride == dstvar.jstride);
        assert(srcvars[n].kstride == dstvar.kstride);
        assert(srcvars[n].nstride == dstvar.nstride);
      }
      const std::ptrdiff_t nstride = dstvar.nstride;
      const std::ptrdiff_t npoints = nstride * ncomps;

      CCTK_REAL *restrict const dstptr = dstvar.dataPtr();
      std::array<const CCTK_REAL *restrict, N> srcptrs;
      for (std::size_t n = 0; n < N; ++n)
        srcptrs[n] = srcvars[n].dataPtr();

#ifndef AMREX_USE_GPU
      // CPU

      const std::ptrdiff_t ntiles = omp_get_max_threads();
      const std::ptrdiff_t tile_size = Arith::align_ceil(
          Arith::div_ceil(npoints, ntiles), std::ptrdiff_t(64));

      for (std::ptrdiff_t imin = 0; imin < npoints; imin += tile_size) {
        using std::min;
        const std::ptrdiff_t imax = min(npoints, imin + tile_size);

        if (!read_dst && N == 1 && factors[0] == 1) {
          // Copy

          auto task = [=]() {
            std::memcpy(&dstptr[imin], &srcptrs[0][imin],
                        (imax - imin) * sizeof *dstptr);
          };
          tasks.push_back(std::move(task));

        } else if (!read_dst && N >= 1 && factors[0] == 1) {
          // Write without scaling

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = srcptrs[0][i];
              for (std::size_t n = 1; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));

        } else if (!read_dst) {
          // Write

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = 0;
              for (std::size_t n = 0; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));

        } else if (scale == 1) {
          // Update without scaling

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = dstptr[i];
              for (std::size_t n = 0; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));

        } else {
          // Update

          auto task = [=]() {
#pragma omp simd
            for (std::ptrdiff_t i = imin; i < imax; ++i) {
              CCTK_REAL accum = scale * dstptr[i];
              for (std::size_t n = 0; n < N; ++n)
                accum += factors[n] * srcptrs[n][i];
              dstptr[i] = accum;
            }
          };
          tasks.push_back(std::move(task));
        }
      } // for imin

#else
      // GPU

      const CCTK_REAL scale1 = scale;
      assert(npoints < INT_MAX);
      const amrex::Box box(
          amrex::IntVect(0, 0, 0), amrex::IntVect(npoints - 1, 0, 0),
          amrex::IntVect(amrex::IndexType::CELL, amrex::IndexType::CELL,
                         amrex::IndexType::CELL));

      if (!read_dst) {

        amrex::launch(
            box, [=] CCTK_DEVICE(const amrex::Box &box) __attribute__((
                     __always_inline__, __flatten__)) {
              const int i = box.smallEnd()[0];
              // const int j = box.smallEnd()[1];
              // const int k = box.smallEnd()[2];
              CCTK_REAL accum = 0;
              // The ROCM 6.2 compiler can't handle
              // `std::array::operator[]`, so we avoid it via pointers:
              // for (std::size_t n = 0; n < N; ++n)
              //   accum += factors[n] * srcptrs[n][i];
              const CCTK_REAL *restrict const factors_ptr = factors.data();
              const CCTK_REAL *restrict const *restrict const srcptrs_ptr =
                  srcptrs.data();
              for (std::size_t n = 0; n < N; ++n)
                accum += factors_ptr[n] * srcptrs_ptr[n][i];
              dstptr[i] = accum;
            });

      } else {

        amrex::launch(
            box, [=] CCTK_DEVICE(const amrex::Box &box) __attribute__((
                     __always_inline__, __flatten__)) {
              const int i = box.smallEnd()[0];
              // const int j = box.smallEnd()[1];
              // const int k = box.smallEnd()[2];
              CCTK_REAL accum = scale1 * dstptr[i];
              // The ROCM 6.2 compiler can't handle
              // `std::array::operator[]`, so we avoid it via pointers:
              // for (std::size_t n = 0; n < N; ++n)
              //   accum += factors[n] * srcptrs[n][i];
              const CCTK_REAL *restrict const factors_ptr = factors.data();
              const CCTK_REAL *restrict const *restrict const srcptrs_ptr =
                  srcptrs.data();
              for (std::size_t n = 0; n < N; ++n)
                accum += factors_ptr[n] * srcptrs_ptr[n][i];
              dstptr[i] = accum;
            });
      }

#endif
    }
  }

#ifndef AMREX_USE_GPU
  // run all tasks
#pragma omp parallel for schedule(dynamic)
  for (std::size_t i = 0; i < tasks.size(); ++i)
    tasks[i]();
#else
  // wait for all tasks
  amrex::Gpu::synchronize();
  AMREX_GPU_ERROR_CHECK();
#endif
}

namespace detail {
template <std::size_t N>
void call_lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                  const std::vector<CCTK_REAL> &factors,
                  const std::vector<const statecomp_t *> &srcs,
                  const std::vector<std::size_t> &indices,
                  const CarpetX::valid_t where) {
  assert(indices.size() == N);
  std::array<CCTK_REAL, N> factors1;
  std::array<const statecomp_t *, N> srcs1;
  for (std::size_t n = 0; n < N; ++n) {
    factors1[n] = factors.at(indices[n]);
    srcs1[n] = srcs.at(indices[n]);
  }
  statecomp_t::lincomb(dst, scale, factors1, srcs1, where);
}

template <std::size_t N>
void call_lincomb_view(
    const statecomp_t &dst, const CCTK_REAL scale,
    const LinearCombinationView<CCTK_REAL, statecomp_t> combination,
    const CarpetX::valid_t where) {
  assert(combination.size == N);
  std::array<CCTK_REAL, N> factors;
  std::array<const statecomp_t *, N> sources;
  for (std::size_t n = 0; n < N; ++n) {
    factors[n] = combination.factors[n];
    sources[n] = combination.sources[n];
  }
  statecomp_t::lincomb(dst, scale, factors, sources, where);
}
} // namespace detail

void statecomp_t::lincomb(
    const statecomp_t &dst, const CCTK_REAL scale,
    const LinearCombinationView<CCTK_REAL, statecomp_t> combination,
    const CarpetX::valid_t where) {
  switch (combination.size) {
  case 0:
    return detail::call_lincomb_view<0>(dst, scale, combination, where);
  case 1:
    return detail::call_lincomb_view<1>(dst, scale, combination, where);
  case 2:
    return detail::call_lincomb_view<2>(dst, scale, combination, where);
  case 3:
    return detail::call_lincomb_view<3>(dst, scale, combination, where);
  case 4:
    return detail::call_lincomb_view<4>(dst, scale, combination, where);
  case 5:
    return detail::call_lincomb_view<5>(dst, scale, combination, where);
  case 6:
    return detail::call_lincomb_view<6>(dst, scale, combination, where);
  case 7:
    return detail::call_lincomb_view<7>(dst, scale, combination, where);
  case 8:
    return detail::call_lincomb_view<8>(dst, scale, combination, where);
  case 9:
    return detail::call_lincomb_view<9>(dst, scale, combination, where);
  case 10:
    return detail::call_lincomb_view<10>(dst, scale, combination, where);
  case 11:
    return detail::call_lincomb_view<11>(dst, scale, combination, where);
  case 12:
    return detail::call_lincomb_view<12>(dst, scale, combination, where);
  case 13:
    return detail::call_lincomb_view<13>(dst, scale, combination, where);
  case 14:
    return detail::call_lincomb_view<14>(dst, scale, combination, where);
  default:
    CCTK_VERROR("Unsupported explicit RK vector length: %d",
                static_cast<int>(combination.size));
  }
}

void statecomp_t::lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                          const std::vector<CCTK_REAL> &factors,
                          const std::vector<const statecomp_t *> &srcs,
                          const CarpetX::valid_t where) {
  const std::size_t N = factors.size();
  assert(srcs.size() == N);

  std::size_t NNZ = 0;
  for (std::size_t n = 0; n < N; ++n)
    NNZ += factors[n] != 0;
  std::vector<std::size_t> indices;
  indices.reserve(NNZ);
  for (std::size_t n = 0; n < N; ++n)
    if (factors[n] != 0)
      indices.push_back(n);
  assert(indices.size() == NNZ);

  switch (NNZ) {
  case 0:
    return detail::call_lincomb<0>(dst, scale, factors, srcs, indices, where);
  case 1:
    return detail::call_lincomb<1>(dst, scale, factors, srcs, indices, where);
  case 2:
    return detail::call_lincomb<2>(dst, scale, factors, srcs, indices, where);
  case 3:
    return detail::call_lincomb<3>(dst, scale, factors, srcs, indices, where);
  case 4:
    return detail::call_lincomb<4>(dst, scale, factors, srcs, indices, where);
  case 5:
    return detail::call_lincomb<5>(dst, scale, factors, srcs, indices, where);
  case 6:
    return detail::call_lincomb<6>(dst, scale, factors, srcs, indices, where);
  case 7:
    return detail::call_lincomb<7>(dst, scale, factors, srcs, indices, where);
  case 8:
    return detail::call_lincomb<8>(dst, scale, factors, srcs, indices, where);
  case 9:
    return detail::call_lincomb<9>(dst, scale, factors, srcs, indices, where);
  case 10:
    return detail::call_lincomb<10>(dst, scale, factors, srcs, indices, where);
  case 11:
    return detail::call_lincomb<11>(dst, scale, factors, srcs, indices, where);
  case 12:
    return detail::call_lincomb<12>(dst, scale, factors, srcs, indices, where);
  case 13:
    return detail::call_lincomb<13>(dst, scale, factors, srcs, indices, where);
  case 14:
    return detail::call_lincomb<14>(dst, scale, factors, srcs, indices, where);
  case 15:
    return detail::call_lincomb<15>(dst, scale, factors, srcs, indices, where);
  case 16:
    return detail::call_lincomb<16>(dst, scale, factors, srcs, indices, where);
  default:
    CCTK_VERROR("Unsupported vector length: %d", (int)NNZ);
  }
}

////////////////////////////////////////////////////////////////////////////////

// Mark groups as invalid
void mark_invalid(const std::vector<int> &groups) {
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const int gi : groups) {
      auto &groupdata = *leveldata.groupdata.at(gi);
      // Invalidate all variables of the current time level
      const int tl = 0;
      for (auto &why_valid : groupdata.valid.at(tl))
        why_valid = CarpetX::why_valid_t(
            [] { return "ODESolvers updated the state vector"; });
    }
  });
}

CarpetX::StageKind
carpetx_stage_kind(const ExplicitRKStageKind stage_kind) {
  switch (stage_kind) {
  case ExplicitRKStageKind::primary:
    return CarpetX::StageKind::primary;
  case ExplicitRKStageKind::fractional:
    return CarpetX::StageKind::fractional;
  case ExplicitRKStageKind::endpoint_probe:
    return CarpetX::StageKind::endpoint_probe;
  }
  throw std::invalid_argument("explicit RK stage kind is invalid");
}

CarpetX::StagePoint
carpetx_stage_point(const ExplicitRKStagePoint &stage_point,
                    const double stage_time) {
  if (stage_point.parent_fraction.denominator <= 0)
    throw std::invalid_argument(
        "explicit RK exact stage fraction denominator is not positive");
  return {carpetx_stage_kind(stage_point.kind), stage_point.stage_index,
          stage_point.stage_count,
          CarpetX::step_clock_t(stage_point.parent_fraction.numerator,
                                stage_point.parent_fraction.denominator),
          stage_time};
}

void prepare_subcycling_stage(
    const ExplicitRKStagePoint &stage_point, const double stage_time,
    const CarpetX::SubcyclingODEMethod active_method,
    const char *const method_name) {
  if (!CarpetX::step_context_active())
    return;

  try {
    CarpetX::prepare_stage(carpetx_stage_point(stage_point, stage_time),
                           active_method);
  } catch (const std::exception &error) {
    CCTK_VERROR("Subcycling stage preparation failed at t=%.17g for ODE "
                "method \"%s\": %s",
                stage_time, method_name, error.what());
  } catch (...) {
    CCTK_VERROR("Subcycling stage preparation failed at t=%.17g for ODE "
                "method \"%s\" with an unknown exception",
                stage_time, method_name);
  }
}

///////////////////////////////////////////////////////////////////////////////

extern "C" void ODESolvers_InitConstants(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_InitConstants;
  DECLARE_CCTK_PARAMETERS;

  *do_substeps = 0;
}

extern "C" void ODESolvers_Solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_Solve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL saved_time = cctkGH->cctk_time;
  const CCTK_REAL old_time = saved_time - dt;
  auto active_method = CarpetX::SubcyclingODEMethod::rk4;
  const CarpetX::StepContext *active_context = nullptr;
  const SubcyclingODEProviderCapability *active_provider = nullptr;
  CarpetX::ScratchStateTransaction *active_transaction = nullptr;
  if (CarpetX::step_context_active()) {
    active_context = CarpetX::current_step_context();
    active_transaction = CarpetX::current_scratch_state_transaction();
    try {
      const auto &candidate =
          require_subcycling_ode_provider(std::string_view(method));
      if (active_context->method != candidate.dense.method)
        throw std::invalid_argument(
            "ODE provider method differs from the active StepContext");

      const auto &tableau = explicit_rk_tableau(candidate.method);
      if (candidate.dense.method != subcycling_method(candidate.method) ||
          candidate.dense.tableau_fingerprint !=
              explicit_rk_tableau_fingerprint(candidate.method) ||
          candidate.dense.endpoint_order != tableau.endpoint_order ||
          candidate.dense.dense_uniform_order < tableau.endpoint_order ||
          candidate.dense.stage_count !=
              static_cast<int>(tableau.a.size()) ||
          !candidate.dense.arbitrary_theta || !candidate.dense.verified)
        throw std::invalid_argument(
            "ODE provider dense capability differs from its exact tableau");

      active_provider = &candidate;
      active_method = candidate.dense.method;
    } catch (const std::exception &error) {
      CCTK_VERROR("Active subcycling ODE provider validation failed for "
                  "method \"%s\": %s",
                  method, error.what());
    } catch (...) {
      CCTK_VERROR("Active subcycling ODE provider validation failed for "
                  "method \"%s\" with an unknown exception",
                  method);
    }

    const auto interval_time_scale =
        std::max({1.0, std::abs(active_context->begin_time),
                  std::abs(active_context->end_time)});
    const auto interval_tolerance =
        16.0 * std::numeric_limits<double>::epsilon() * interval_time_scale;
    if (!std::isfinite(old_time) || !std::isfinite(saved_time) ||
        !std::isfinite(dt) ||
        std::abs(old_time - active_context->begin_time) > interval_tolerance ||
        std::abs(saved_time - active_context->end_time) > interval_tolerance ||
        std::abs(dt -
                 (active_context->end_time - active_context->begin_time)) >
            interval_tolerance)
      CCTK_VERROR("ODE solver interval [%g,%g] with dt=%g does not match "
                  "active subcycling StepContext [%g,%g]",
                  double(old_time), double(saved_time), double(dt),
                  active_context->begin_time, active_context->end_time);
  }

  static bool did_output = false;
  if (verbose || !did_output)
    CCTK_VINFO("ODE integrator is %s", method);
  did_output = true;

  static CarpetX::Timer timer("ODESolvers::Solve");
  CarpetX::Interval interval(timer);

  const int tl = 0;

  static CarpetX::Timer timer_setup("ODESolvers::Solve::setup");
  std::optional<CarpetX::Interval> interval_setup(timer_setup);

  GroupSchemaBuild group_schema;
  try {
    group_schema =
        build_cactus_group_schema(carpetx_subcycling_enabled());
  } catch (const std::exception &error) {
    CCTK_VERROR("ODE evolved/RHS group schema construction failed: %s",
                error.what());
  } catch (...) {
    CCTK_ERROR("ODE evolved/RHS group schema construction failed with an "
               "unknown exception");
  }

  statecomp_t var, rhs;
  const auto &ordered_group_pairs =
      group_schema.contract.ordered_group_pairs;
  const auto &var_groups = group_schema.evolved_groups;
  const auto &rhs_groups = group_schema.rhs_groups;
  const auto &dep_groups = group_schema.contract.dependent_groups;
  const int nvars = group_schema.evolved_variable_count;
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &pair : ordered_group_pairs) {
      auto &groupdata = *leveldata.groupdata.at(pair.evolved_group);
      auto &rhs_groupdata = *leveldata.groupdata.at(pair.rhs_group);
      assert(groupdata.groupindex == pair.evolved_group);
      assert(rhs_groupdata.groupindex == pair.rhs_group);
      assert(rhs_groupdata.numvars == groupdata.numvars);
      var.groupdatas.push_back(&groupdata);
      var.mfabs.push_back(groupdata.mfab.at(tl).get());
      rhs.groupdatas.push_back(&rhs_groupdata);
      rhs.mfabs.push_back(rhs_groupdata.mfab.at(tl).get());
    }
  });
  if (verbose)
    CCTK_VINFO("  Integrating %d variables", nvars);
  if (nvars == 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "Integrating %d variables", nvars);

  bool scratch_group_pairs_match = true;
  if (active_transaction != nullptr) {
    const auto &certified_pairs = active_transaction->group_pairs();
    scratch_group_pairs_match =
        certified_pairs.size() == ordered_group_pairs.size();
    for (std::size_t pair = 0;
         scratch_group_pairs_match && pair < ordered_group_pairs.size();
         ++pair)
      scratch_group_pairs_match =
          certified_pairs[pair].evolved_group ==
              ordered_group_pairs[pair].evolved_group &&
          certified_pairs[pair].rhs_group ==
              ordered_group_pairs[pair].rhs_group;
  }

  const bool scratch_dependent_groups_match =
      active_transaction == nullptr ||
      active_transaction->dependent_groups() == dep_groups;

  for (const int gi : var_groups)
    assert(std::find(dep_groups.begin(), dep_groups.end(), gi) ==
           dep_groups.end());
  for (const int gi : rhs_groups)
    assert(std::find(var_groups.begin(), var_groups.end(), gi) ==
           var_groups.end());

  interval_setup.reset();

  {
    static CarpetX::Timer timer_alloc_temps("ODESolvers::Solve::alloc_temps");
    CarpetX::Interval interval_alloc_temps(timer_alloc_temps);
    statecomp_t::init_tmp_mfabs();
  }

  static CarpetX::Timer timer_lincomb("ODESolvers::Solve::lincomb");
  static CarpetX::Timer timer_rhs("ODESolvers::Solve::rhs");
  static CarpetX::Timer timer_poststep("ODESolvers::Solve::poststep");

  const auto copy_state = [](const auto &var, const CarpetX::valid_t where) {
    return var.copy(where);
  };
  const auto evaluate_rhs = [&](const int n) {
    CarpetX::Interval interval_rhs(timer_rhs);
    if (verbose)
      CCTK_VINFO("Calculating RHS #%d at t=%g", n, double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
  };
  const auto validate_rhs = [&](const int) {
    rhs.check_valid(CarpetX::make_valid_int(),
                    "ODESolvers after calling ODESolvers_RHS");
  };
  const auto calcrhs = [&](const int n) {
    evaluate_rhs(n);
    validate_rhs(n);
  };
  const auto calcupdate_at_time =
      [&](const int n, const CCTK_REAL stage_time, const CCTK_REAL a0,
          const auto &as, const auto &vars) {
        {
          CarpetX::Interval interval_lincomb(timer_lincomb);
          statecomp_t::lincomb(var, a0, as, vars,
                               CarpetX::make_valid_int());
          var.check_valid(CarpetX::make_valid_int(),
                          "ODESolvers after defining new state vector");
          mark_invalid(dep_groups);
        }
        {
          CarpetX::Interval interval_poststep(timer_poststep);
          *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = stage_time;
          CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
          if (verbose)
            CCTK_VINFO("Calculated new state #%d at t=%g", n,
                       double(cctkGH->cctk_time));
        }
      };
  // t = t_0 + c
  // var = a_0 * var + \Sum_i a_i * var_i
  const auto calcupdate = [&](const int n, const CCTK_REAL c,
                              const CCTK_REAL a0, const auto &as,
                              const auto &vars) {
    calcupdate_at_time(n, old_time + c, a0, as, vars);
  };

  const bool uses_extracted_explicit_rk =
      active_provider != nullptr ||
      (active_context == nullptr &&
       (CCTK_EQUALS(method, "RK4") || CCTK_EQUALS(method, "RKF78") ||
        CCTK_EQUALS(method, "DP87")));
  if (!uses_extracted_explicit_rk) {
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
  }

  if (CCTK_EQUALS(method, "constant")) {

    // y1 = y0

    // do nothing

  } else if (CCTK_EQUALS(method, "Euler")) {

    // k1 = f(y0)
    // y1 = y0 + h k1

    calcrhs(1);
    calcupdate(1, dt, 1.0, reals<1>{dt}, states<1>{&rhs});

  } else if (CCTK_EQUALS(method, "RK2")) {

    // k1 = f(y0)
    // k2 = f(y0 + h/2 k1)
    // y1 = y0 + h k2

    const auto old = copy_state(var, CarpetX::make_valid_all());

    calcrhs(1);
    calcupdate(1, dt / 2, 1.0, reals<1>{dt / 2}, states<1>{&rhs});

    calcrhs(2);
    calcupdate(2, dt, 0.0, reals<2>{1.0, dt}, states<2>{&old, &rhs});

  } else if (CCTK_EQUALS(method, "RK3")) {

    // k1 = f(y0)
    // k2 = f(y0 + h/2 k1)
    // k3 = f(y0 - h k1 + 2 h k2)
    // y1 = y0 + h/6 k1 + 2/3 h k2 + h/6 k3

    const auto old = copy_state(var, CarpetX::make_valid_all());

    calcrhs(1);
    const auto k1 = copy_state(rhs, CarpetX::make_valid_int());
    calcupdate(1, dt / 2, 1.0, reals<1>{dt / 2}, states<1>{&k1});

    calcrhs(2);
    const auto k2 = copy_state(rhs, CarpetX::make_valid_int());
    calcupdate(2, dt, 0.0, reals<3>{1.0, -dt, 2 * dt},
               states<3>{&old, &k1, &k2});

    calcrhs(3);
    calcupdate(3, dt, 0.0, reals<4>{1.0, dt / 6, 2 * dt / 3, dt / 6},
               states<4>{&old, &k1, &k2, &rhs});

  } else if (CCTK_EQUALS(method, "SSPRK3")) {

    // k1 = f(y0)
    // k2 = f(y0 + h k1)
    // k3 = f(y0 + h/4 k1 + h/4 k2)
    // y1 = y0 + h/6 k1 + h/6 k2 + 2/3 h k3

    const auto old = copy_state(var, CarpetX::make_valid_all());

    calcrhs(1);
    const auto k1 = copy_state(rhs, CarpetX::make_valid_int());
    calcupdate(1, dt, 1.0, reals<1>{dt}, states<1>{&k1});

    calcrhs(2);
    const auto k2 = copy_state(rhs, CarpetX::make_valid_int());
    calcupdate(2, dt / 2, 0.0, reals<3>{1.0, dt / 4, dt / 4},
               states<3>{&old, &k1, &k2});

    calcrhs(3);
    calcupdate(3, dt, 0.0, reals<4>{1.0, dt / 6, dt / 6, 2 * dt / 3},
               states<4>{&old, &k1, &k2, &rhs});

  } else if (uses_extracted_explicit_rk) {

    ExplicitRKMethod explicit_method = ExplicitRKMethod::rk4;
    if (active_provider != nullptr) {
      explicit_method = active_provider->method;
    } else {
      if (CCTK_EQUALS(method, "RK4"))
        explicit_method = ExplicitRKMethod::rk4;
      else if (CCTK_EQUALS(method, "RKF78"))
        explicit_method = ExplicitRKMethod::rkf78;
      else if (CCTK_EQUALS(method, "DP87"))
        explicit_method = ExplicitRKMethod::dp87;
      else
        CCTK_VERROR(
            "Internal explicit RK method dispatch failure for \"%s\"",
            method);
    }

    int explicit_update_index = 0;
    CactusExplicitRKOperations operations{
        var,
        rhs,
        [&](const CCTK_REAL stage_time) {
          *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = stage_time;
        },
        evaluate_rhs,
        validate_rhs,
        [&](const int update_index, const CCTK_REAL stage_time,
            const CCTK_REAL destination_scale,
            const LinearCombinationView<CCTK_REAL, statecomp_t> combination) {
          explicit_update_index = update_index;
          {
            CarpetX::Interval interval_lincomb(timer_lincomb);
            statecomp_t::lincomb(var, destination_scale, combination,
                                 CarpetX::make_valid_int());
            var.check_valid(CarpetX::make_valid_int(),
                            "ODESolvers after defining new state vector");
            mark_invalid(dep_groups);
          }
          *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = stage_time;
        },
        [&](statecomp_t &accumulator, const CCTK_REAL factor,
            const statecomp_t &increment) {
          CarpetX::Interval interval_lincomb(timer_lincomb);
          statecomp_t::lincomb(accumulator, 1.0, reals<1>{factor},
                               states<1>{&increment},
                               CarpetX::make_valid_int());
        }};
    operations.stage_preparation_callback =
        [&](const ExplicitRKStagePoint &stage_point,
            const CCTK_REAL stage_time) {
          *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = stage_time;
          prepare_subcycling_stage(stage_point, stage_time, active_method,
                                   method);
        };
    operations.live_post_step_callback = [&](const CCTK_REAL stage_time) {
      CarpetX::Interval interval_poststep(timer_poststep);
      *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = stage_time;
      CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
      if (verbose)
        CCTK_VINFO("Calculated new state #%d at t=%g",
                   explicit_update_index, double(cctkGH->cctk_time));
    };

    try {
      if (active_transaction != nullptr && !scratch_group_pairs_match)
        throw std::invalid_argument(
            "ODESolvers ordered evolved/RHS pairs differ from the active "
            "scratch transaction");
      if (active_transaction != nullptr && !scratch_dependent_groups_match)
        throw std::invalid_argument(
            "ODESolvers dependent groups differ from the active scratch "
            "transaction");
      if (active_transaction == nullptr) {
        advance_explicit_rk(explicit_method, old_time, dt,
                            InitialRHSMode::calculate, operations);
      } else if (active_context != nullptr &&
                 !active_context->require_dense_output) {
        advance_explicit_rk(explicit_method, old_time, dt,
                            InitialRHSMode::calculate, operations);
      } else {
        // Only a parent level needs the independent dense-extension solves.
        if (active_context == nullptr)
          throw std::logic_error(
              "scratch transaction has no active StepContext");

        TransactionPrimaryObserver primary_observer(*active_transaction);
        advance_explicit_rk(explicit_method, old_time, dt,
                            InitialRHSMode::calculate, operations,
                            primary_observer);
        auto primary = primary_observer.take_complete();

        auto left_state =
            statecomp_t::from_scratch(std::move(primary.left_state));
        auto left_rhs = statecomp_t::from_scratch(std::move(primary.left_rhs));
        auto accepted_endpoint = statecomp_t::from_scratch(
            std::move(primary.accepted_endpoint));
        auto scratch_var = left_state.snapshot_state();
        auto scratch_rhs = left_rhs.snapshot_rhs();

        CactusExplicitRKOperations scratch_operations{
            scratch_var,
            scratch_rhs,
            [](const CCTK_REAL) {},
            [](const int) {},
            [&](const int) {
              scratch_rhs.check_valid(
                  CarpetX::make_valid_int(),
                  "ODESolvers scratch RHS after certified execution");
            },
            [&](const int, const CCTK_REAL,
                const CCTK_REAL destination_scale,
                const LinearCombinationView<CCTK_REAL, statecomp_t>
                    combination) {
              statecomp_t::lincomb(scratch_var, destination_scale, combination,
                                   CarpetX::make_valid_int());
              scratch_var.check_valid(
                  CarpetX::make_valid_int(),
                  "ODESolvers after defining scratch state vector");
            },
            [&](statecomp_t &accumulator, const CCTK_REAL factor,
                const statecomp_t &increment) {
              statecomp_t::lincomb(accumulator, 1.0, reals<1>{factor},
                                   states<1>{&increment},
                                   CarpetX::make_valid_int());
            }};
        scratch_operations.stage_preparation_callback =
            [&](const ExplicitRKStagePoint &stage_point,
                const CCTK_REAL stage_time) {
              *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = stage_time;
              prepare_subcycling_stage(stage_point, stage_time, active_method,
                                       method);
            };
        scratch_operations.stage_materialization_callback =
            [&](statecomp_t &state,
                const ExplicitRKStagePoint &, const CCTK_REAL) {
              if (!state.is_scratch() ||
                  &state.scratch().transaction() != active_transaction)
                throw std::invalid_argument(
                    "scratch stage materialization transaction owner differs");
              active_transaction->restore_state(state.scratch().token());
            };

        using ScratchHooks =
            CactusExplicitRKScratchHooks<CCTK_REAL, statecomp_t>;
        ScratchHooks scratch_hooks;
        scratch_hooks.restore_left =
            [&](statecomp_t &destination_state,
                statecomp_t &destination_rhs,
                const statecomp_t &source_state,
                const statecomp_t &source_rhs, const CCTK_REAL) {
              if (!source_state.is_scratch() || !source_rhs.is_scratch() ||
                  &source_state.scratch().transaction() != active_transaction ||
                  &source_rhs.scratch().transaction() != active_transaction)
                throw std::invalid_argument(
                    "scratch restore_left transaction owner differs");
              active_transaction->restore_left(source_state.scratch().token(),
                                               source_rhs.scratch().token());
              destination_state = source_state.snapshot_state();
              destination_rhs = source_rhs.snapshot_rhs();
            };
        scratch_hooks.restore_state =
            [&](statecomp_t &destination, const statecomp_t &source,
                const CCTK_REAL) {
              if (!source.is_scratch() ||
                  &source.scratch().transaction() != active_transaction)
                throw std::invalid_argument(
                    "scratch restore_state transaction owner differs");
              active_transaction->restore_state(source.scratch().token());
              destination = source.snapshot_state();
            };
        scratch_hooks.post_step_after_update =
            [&](statecomp_t &state,
                const ExplicitRKStagePoint &stage_point,
                const CCTK_REAL stage_time) {
              active_transaction->post_step_after_update(
                  *active_context,
                  carpetx_stage_point(stage_point, stage_time));
              state.replace_scratch(TransactionStateBackend::from_token(
                  *active_transaction,
                  active_transaction->capture_scratch_evolved()));
            };
        scratch_hooks.evaluate_rhs =
            [&](statecomp_t &, statecomp_t &rhs_state, const int,
                const ExplicitRKStagePoint &stage_point,
                const CCTK_REAL stage_time) {
              active_transaction->evaluate_rhs(
                  *active_context,
                  carpetx_stage_point(stage_point, stage_time));
              rhs_state.replace_scratch(TransactionStateBackend::from_token(
                  *active_transaction,
                  active_transaction->capture_scratch_rhs()));
            };
        scratch_hooks.probe_endpoint_rhs =
            [&](statecomp_t &, statecomp_t &rhs_state,
                const ExplicitRKStagePoint &stage_point,
                const CCTK_REAL stage_time) {
              active_transaction->evaluate_rhs(
                  *active_context,
                  carpetx_stage_point(stage_point, stage_time));
              rhs_state.replace_scratch(TransactionStateBackend::from_token(
                  *active_transaction,
                  active_transaction->capture_scratch_rhs()));
              return rhs_state.snapshot_rhs();
            };
        scratch_hooks.rhs_evaluation_count =
            [&] { return active_transaction->rhs_evaluation_count(); };
        scratch_operations.scratch_hooks = &scratch_hooks;

        auto samples = collect_reference_dense_samples(
            explicit_method, old_time, dt, left_state, left_rhs,
            accepted_endpoint, scratch_operations);
        if (active_provider == nullptr)
          throw std::logic_error(
              "scratch dense output has no active ODE provider");
        CarpetX::DenseOutputProvider provider(active_provider->dense);
        std::vector<CarpetX::ScratchDenseSampleRef> sample_refs;
        sample_refs.reserve(samples.size());
        for (const auto &sample : samples.samples()) {
          if (!sample.payload.is_scratch())
            throw std::logic_error(
                "reference dense sample escaped the scratch backend");
          const auto kind =
              sample.kind == CarpetX::DenseSampleKind::value
                  ? CarpetX::ScratchDenseSampleKind::value
                  : CarpetX::ScratchDenseSampleKind::raw_derivative;
          sample_refs.push_back(
              {sample.theta, kind, &sample.payload.scratch().token()});
        }
        const CarpetX::DenseIntervalId interval_id{
            active_context->level,
            active_context->begin_clock,
            active_context->end_clock,
            active_context->begin_time,
            active_context->end_time,
            active_provider->dense.method,
            active_provider->dense.tableau_fingerprint};
        active_transaction->commit_dense(*active_context, provider,
                                         interval_id, sample_refs);
      }
    } catch (const std::exception &error) {
      if (active_transaction != nullptr &&
          !active_transaction->discarded())
        active_transaction->discard();
      CCTK_VERROR("Explicit RK method \"%s\" failed: %s", method,
                  error.what());
    } catch (...) {
      if (active_transaction != nullptr &&
          !active_transaction->discarded())
        active_transaction->discard();
      CCTK_VERROR("Explicit RK method \"%s\" failed with an unknown exception",
                  method);
    }

  } else if (CCTK_EQUALS(method, "IMEX122") ||
             CCTK_EQUALS(method, "Implicit Euler")) {
    // Implicit definition:
    //   y1 = y0 + h/2 f(y0) + h/2 g(y1)
    //   y2 = y0 + h f(y1) + h g(y1)

    // Implicit RHS:
    //   u1 = G(u0, h)   where   u1 = u0 + h g(u1)

    // Explicit definition:
    //   k1 = f(y0)
    //   y1 = G(y0 + h/2 k1, h/2)
    //   k'2 = (y1 - y0 - h/2 k1) / (h/2)
    //   k2 = f(y1)
    //   y2 = y0 + h k2 + h k'2

    const auto y0 = var.copy(CarpetX::make_valid_int /*all*/ ());

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
    if (verbose)
      CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    const auto k1 = rhs.copy(CarpetX::make_valid_int());

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    statecomp_t::lincomb(var, 1, make_array(dt / 2), make_array(&rhs),
                         CarpetX::make_valid_int());
    var.check_valid(CarpetX::make_valid_int(),
                    "ODESolvers after defining new state vector");
    mark_invalid(dep_groups);
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_delta_time) = dt / 2;
    if (verbose)
      CCTK_VINFO("Taking implicit step #1 at t=%g with dt=%g",
                 double(cctkGH->cctk_time), double(cctkGH->cctk_delta_time));
    CallScheduleGroup(cctkGH, "ODESolvers_ImplicitStep");
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_delta_time) = dt;

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
    const auto y1 = var.copy(CarpetX::make_valid_int /*all*/ ());

    statecomp_t kprime2;
    statecomp_t::lincomb(kprime2, 0,
                         make_array(-CCTK_REAL(1), +CCTK_REAL(1), -dt / 2),
                         make_array(&y0, &y1, &k1), CarpetX::make_valid_int());

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt;
    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    const auto k2 = rhs.copy(CarpetX::make_valid_int());

    statecomp_t::lincomb(var, 0, make_array(CCTK_REAL(1), dt, dt),
                         make_array(&y0, &k2, &kprime2),
                         CarpetX::make_valid_int());
    var.check_valid(CarpetX::make_valid_int(),
                    "ODESolvers after defining new state vector");
    mark_invalid(dep_groups);

  } else {
    assert(0);
  }

  {
    static CarpetX::Timer timer_free_temps("ODESolvers::Solve::free_temps");
    CarpetX::Interval interval_free_temps(timer_free_temps);
    statecomp_t::free_tmp_mfabs();
  }

  // Reset current time
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = saved_time;

  // TODO: Update time here, and not during time level cycling in the driver
}

} // namespace ODESolvers
