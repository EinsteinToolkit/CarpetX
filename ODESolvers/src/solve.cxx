// TODO: Don't include files from other thorns; create a proper interface
#include "../../CarpetX/src/driver.hxx"
#include "../../CarpetX/src/schedule.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>
#include <util_Table.h>

#include <AMReX_MultiFab.H>

#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

CCTK_REAL ODESolvers_alpha = 0;

namespace ODESolvers {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// Taken from <https://en.cppreference.com/w/cpp/experimental/make_array>
namespace details {
template <class> struct is_ref_wrapper : std::false_type {};
template <class T>
struct is_ref_wrapper<std::reference_wrapper<T> > : std::true_type {};

template <class T>
using not_ref_wrapper = std::negation<is_ref_wrapper<std::decay_t<T> > >;

template <class D, class...> struct return_type_helper { using type = D; };
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

  statecomp_t() = delete;
  statecomp_t(statecomp_t &&) = default;
  statecomp_t &operator=(statecomp_t &&) = default;

  statecomp_t(const cGH *cctkGH) : cctkGH(cctkGH) {}

  // Don't allow copies because we might own stuff
  statecomp_t(const statecomp_t &) = delete;
  statecomp_t &operator=(const statecomp_t &) = delete;

  const cGH *cctkGH;
  vector<string> groupnames;
  vector<int> groupids;
  vector<amrex::MultiFab *> mfabs;

private:
  // These will be automatically freed when this object is deallocated
  vector<unique_ptr<amrex::MultiFab> > owned_stuff;

public:
  void check_valid(const function<string()> &why) const;
  void check_valid(const string &why) const {
    check_valid([=]() { return why; });
  }

  statecomp_t copy() const;

  template <size_t N>
  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const array<CCTK_REAL, N> &factors,
                      const array<const statecomp_t *, N> &srcs);
  template <size_t N>
  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const array<CCTK_REAL, N> &factors,
                      const array<statecomp_t *, N> &srcs) {
    array<const statecomp_t *, N> srcs1;
    for (size_t n = 0; n < N; ++n)
      srcs1[n] = srcs[n];
    lincomb(dst, scale, factors, srcs1);
  }

  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const vector<CCTK_REAL> &factors,
                      const vector<const statecomp_t *> &srcs);
};

////////////////////////////////////////////////////////////////////////////////

// Ensure a state vector has valid data everywhere
void statecomp_t::check_valid(const function<string()> &why) const {
  for (const int groupid : groupids) {
    if (groupid >= 0) {
      assert(CarpetX::active_levels);
      CarpetX::active_levels->loop([&](const auto &leveldata) {
        const auto &groupdata = *leveldata.groupdata.at(groupid);
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          const int tl = 0;
          CarpetX::check_valid(groupdata, vi, tl, nan_handling_t::forbid_nans,why);
        }
      });
    }
  }
}

// Copy state vector into newly allocated memory
statecomp_t statecomp_t::copy() const {
  check_valid("before copy");
  const size_t size = mfabs.size();
  statecomp_t result(cctkGH);
  result.groupnames = groupnames;
  result.groupids.resize(size, -1);
  result.mfabs.reserve(size);
  result.owned_stuff.reserve(size);
  for (size_t n = 0; n < size; ++n) {
    const auto &x = mfabs.at(n);
#ifdef CCTK_DEBUG
    if (x->contains_nan())
      CCTK_VERROR("statecomp_t::copy.x: Group %s contains nans",
                  groupnames.at(n).c_str());
#endif
    auto y = make_unique<amrex::MultiFab>(x->boxArray(), x->DistributionMap(),
                                          x->nComp(), x->nGrowVect());
    amrex::MultiFab::Copy(*y, *x, 0, 0, y->nComp(), y->nGrowVect());
#ifdef CCTK_DEBUG
    if (y->contains_nan())
      CCTK_VERROR("statecomp_t::copy.y: Group %s contains nans",
                  result.groupnames.at(n).c_str());
#endif
    result.mfabs.push_back(y.get());
    result.owned_stuff.push_back(move(y));
  }
  result.check_valid("after copy");
  return result;
}

template <size_t N>
void statecomp_t::lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                          const array<CCTK_REAL, N> &factors,
                          const array<const statecomp_t *, N> &srcs) {
  const size_t size = dst.mfabs.size();
  for (size_t n = 0; n < N; ++n)
    assert(srcs[n]->mfabs.size() == size);
  for (size_t m = 0; m < size; ++m) {
    const auto ncomp = dst.mfabs.at(m)->nComp();
    const auto ngrowvect = dst.mfabs.at(m)->nGrowVect();
    for (size_t n = 0; n < N; ++n) {
      assert(srcs[n]->mfabs.at(m)->nComp() == ncomp);
      assert(srcs[n]->mfabs.at(m)->nGrowVect() == ngrowvect);
    }
  }

  const bool read_dst = scale != 0;
  if (read_dst)
    dst.check_valid("before lincomb, destination");
  for (size_t n = 0; n < N; ++n)
    srcs[n]->check_valid([=]() {
      ostringstream buf;
      buf << "before lincomb, source #" << n;
      return buf.str();
    });

#ifndef AMREX_USE_GPU
  vector<function<void()> > tasks;
#endif

  for (size_t m = 0; m < size; ++m) {
    const size_t ncomp = dst.mfabs.at(m)->nComp();
    const auto mfitinfo = amrex::MFItInfo().DisableDeviceSync();
    for (amrex::MFIter mfi(*dst.mfabs.at(m), mfitinfo); mfi.isValid(); ++mfi) {
      const amrex::Array4<CCTK_REAL> dstvar = dst.mfabs.at(m)->array(mfi);
      array<amrex::Array4<const CCTK_REAL>, N> srcvars;
      for (size_t n = 0; n < N; ++n)
        srcvars[n] = srcs[n]->mfabs.at(m)->const_array(mfi);
      for (size_t n = 0; n < N; ++n) {
        assert(srcvars[n].jstride == dstvar.jstride);
        assert(srcvars[n].kstride == dstvar.kstride);
        assert(srcvars[n].nstride == dstvar.nstride);
      }
      const ptrdiff_t nstride = dstvar.nstride;

      CCTK_REAL *restrict const dstptr = dstvar.dataPtr();
      array<const CCTK_REAL *restrict, N> srcptrs;
      for (size_t n = 0; n < N; ++n)
        srcptrs[n] = srcvars[n].dataPtr();
      const ptrdiff_t npoints = nstride * ncomp;

#ifndef AMREX_USE_GPU
      // CPU

      if (!read_dst) {

        auto task = [=]() {
#pragma omp simd
          for (ptrdiff_t i = 0; i < npoints; ++i) {
            CCTK_REAL accum = 0;
            for (size_t n = 0; n < N; ++n)
              accum += factors[n] * srcptrs[n][i];
            dstptr[i] = accum;
          }
        };
        tasks.emplace_back(move(task));

      } else {

        auto task = [=]() {
#pragma omp simd
          for (ptrdiff_t i = 0; i < npoints; ++i) {
            CCTK_REAL accum = scale * dstptr[i];
            for (size_t n = 0; n < N; ++n)
              accum += factors[n] * srcptrs[n][i];
            dstptr[i] = accum;
          }
        };
        tasks.emplace_back(move(task));
      }

#else
      // GPU

      const CCTK_REAL scale1 = scale;
      assert(npoints < INT_MAX);
      const amrex::Box box(
          amrex::IntVect(0, 0, 0), amrex::IntVect(npoints - 1, 0, 0),
          amrex::IntVect(amrex::IndexType::CELL, amrex::IndexType::CELL,
                         amrex::IndexType::CELL));

      if (!read_dst) {

        amrex::launch(box, [=] CCTK_DEVICE(const amrex::Box &box)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 const int i = box.smallEnd()[0];
                                 // const int j = box.smallEnd()[1];
                                 // const int k = box.smallEnd()[2];
                                 CCTK_REAL accum = 0;
                                 for (size_t n = 0; n < N; ++n)
                                   accum += factors[n] * srcptrs[n][i];
                                 dstptr[i] = accum;
                               });

      } else {

        amrex::launch(box, [=] CCTK_DEVICE(const amrex::Box &box)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 const int i = box.smallEnd()[0];
                                 // const int j = box.smallEnd()[1];
                                 // const int k = box.smallEnd()[2];
                                 CCTK_REAL accum = scale1 * dstptr[i];
                                 for (size_t n = 0; n < N; ++n)
                                   accum += factors[n] * srcptrs[n][i];
                                 dstptr[i] = accum;
                               });
      }

#endif
    }
  }

#ifndef AMREX_USE_GPU
  // run all tasks
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < tasks.size(); ++i)
    tasks[i]();
#else
  // wait for all tasks
  amrex::Gpu::synchronize();
  AMREX_GPU_ERROR_CHECK();
#endif

  dst.check_valid("after lincomb, destination");
}

namespace detail {
template <size_t N>
void call_lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                  const vector<CCTK_REAL> &factors,
                  const vector<const statecomp_t *> &srcs,
                  const vector<size_t> &indices) {
  assert(indices.size() == N);
  array<CCTK_REAL, N> factors1;
  array<const statecomp_t *, N> srcs1;
  for (size_t n = 0; n < N; ++n) {
    factors1[n] = factors.at(indices[n]);
    srcs1[n] = srcs.at(indices[n]);
  }
  statecomp_t::lincomb(dst, scale, factors1, srcs1);
}
} // namespace detail

void statecomp_t::lincomb(const statecomp_t &dst, const CCTK_REAL scale,
                          const vector<CCTK_REAL> &factors,
                          const vector<const statecomp_t *> &srcs) {
  const size_t N = factors.size();
  assert(srcs.size() == N);

  size_t NNZ = 0;
  for (size_t n = 0; n < N; ++n)
    NNZ += factors[n] != 0;
  vector<size_t> indices;
  indices.reserve(NNZ);
  for (size_t n = 0; n < N; ++n)
    if (factors[n] != 0)
      indices.push_back(n);
  assert(indices.size() == NNZ);

  switch (NNZ) {
  case 0:
    return detail::call_lincomb<0>(dst, scale, factors, srcs, indices);
  case 1:
    return detail::call_lincomb<1>(dst, scale, factors, srcs, indices);
  case 2:
    return detail::call_lincomb<2>(dst, scale, factors, srcs, indices);
  case 3:
    return detail::call_lincomb<3>(dst, scale, factors, srcs, indices);
  case 4:
    return detail::call_lincomb<4>(dst, scale, factors, srcs, indices);
  case 5:
    return detail::call_lincomb<5>(dst, scale, factors, srcs, indices);
  case 6:
    return detail::call_lincomb<6>(dst, scale, factors, srcs, indices);
  case 7:
    return detail::call_lincomb<7>(dst, scale, factors, srcs, indices);
  case 8:
    return detail::call_lincomb<8>(dst, scale, factors, srcs, indices);
  case 9:
    return detail::call_lincomb<9>(dst, scale, factors, srcs, indices);
  default:
    CCTK_ERROR("Unsupported vector length");
  }
}

////////////////////////////////////////////////////////////////////////////////

int get_group_rhs(const int gi) {
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  vector<char> rhs_buf(1000);
  const int iret =
      Util_TableGetString(tags, rhs_buf.size(), rhs_buf.data(), "rhs");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    rhs_buf[0] = '\0'; // default: empty (no RHS)
  } else if (iret >= 0) {
    // do nothing
  } else {
    assert(0);
  }

  const string str(rhs_buf.data());
  if (str.empty())
    return -1; // No RHS specified

  auto str1 = str;
  if (str1.find(':') == string::npos) {
    const char *impl = CCTK_GroupImplementationI(gi);
    str1 = string(impl) + "::" + str1;
  }
  const int gi1 = CCTK_GroupIndex(str1.c_str());
  assert(gi1 >= 0); // Check RHSs are valid groups
  const int flux = gi1;

  assert(flux != gi);

  return flux;
}

///////////////////////////////////////////////////////////////////////////////

extern "C" void ODESolvers_Solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_Solve;
  DECLARE_CCTK_PARAMETERS;

  static bool did_output = false;
  if (verbose || !did_output)
    CCTK_VINFO("Integrator is %s", method);
  did_output = true;

  const CCTK_REAL dt = cctk_delta_time;
  const int tl = 0;

  statecomp_t var(cctkGH), rhs(cctkGH);
  int nvars = 0;
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      // TODO: add support for evolving grid scalars
      if (groupdataptr == nullptr)
        continue;

      const auto &groupdata = *groupdataptr;
      const int rhs_gi = get_group_rhs(groupdata.groupindex);
      if (rhs_gi >= 0) {
        assert(rhs_gi != groupdata.groupindex);
        const auto &rhs_groupdata = *leveldata.groupdata.at(rhs_gi);
        assert(rhs_groupdata.numvars == groupdata.numvars);
        var.groupnames.push_back(CCTK_GroupName(groupdata.groupindex));
        var.groupids.push_back(groupdata.groupindex);
        var.mfabs.push_back(groupdata.mfab.at(tl).get());
        rhs.groupnames.push_back(CCTK_GroupName(rhs_groupdata.groupindex));
        rhs.groupids.push_back(rhs_groupdata.groupindex);
        rhs.mfabs.push_back(rhs_groupdata.mfab.at(tl).get());
        if (leveldata.level == active_levels->min_level)
          nvars += groupdata.numvars;
      }
    }
  });
  if (verbose)
    CCTK_VINFO("  Integrating %d variables", nvars);
  if (nvars == 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "Integrating %d variables", nvars);

  const CCTK_REAL saved_time = cctkGH->cctk_time;
  const CCTK_REAL old_time = cctkGH->cctk_time - dt;

  if (CCTK_EQUALS(method, "constant")) {

    // y1 = y0

    // do nothing

  } else if (CCTK_EQUALS(method, "Euler")) {

    // k1 = f(y0)
    // y1 = y0 + h k1

    // Calculate first RHS
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
    if (verbose)
      CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 1, make_array(dt), make_array(&rhs));

  } else if (CCTK_EQUALS(method, "RK2")) {

    // k1 = f(y0)
    // k2 = f(y0 + h/2 k1)
    // y1 = y0 + h k2

    const auto old = var.copy();
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;

    // Step 1
    if (verbose)
      CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    // rhs = RHS(var)
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Add scaled RHS to state vector
    // lincomb(dest, alpha, beta^i, src^i)
    // dest = alpha * dest + beta^i * src^i
    statecomp_t::lincomb(var, 1, make_array(dt / 2), make_array(&rhs));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    // var = poststep(var)
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    // Step 2
    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    // rhs = RHS(var)
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Calculate new state vector
    statecomp_t::lincomb(var, 0, make_array(1.0, dt), make_array(&old, &rhs));

  } else if (CCTK_EQUALS(method, "RK3")) {

    // k1 = f(y0)
    // k2 = f(y0 + h/2 k1)
    // k3 = f(y0 - h k1 + 2 h k2)
    // y1 = y0 + h/6 k1 + 2/3 h k2 + h/6 k3

    const auto old = var.copy();
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;

    // Step 1
    if (verbose)
      CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    const auto k1 = rhs.copy();

    // Step 2

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 1, make_array(dt / 2), make_array(&k1));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto k2 = rhs.copy();

    // Step 3

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 0, make_array(1.0, -dt, 2 * dt),
                         make_array(&old, &k1, &k2));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #3 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto &k3 = rhs;

    // Calculate new state vector
    statecomp_t::lincomb(var, 0, make_array(1.0, dt / 6, 2 * dt / 3, dt / 6),
                         make_array(&old, &k1, &k2, &k3));

  } else if (CCTK_EQUALS(method, "SSPRK3")) {

    // k1 = f(y0)
    // k2 = f(y0 + h k1)
    // k3 = f(y0 + h/4 k1 + h/4 k2)
    // y1 = y0 + h/6 k1 + h/6 k2 + 2/3 h k3

    const auto old = var.copy();
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;

    // Step 1
    if (verbose)
      CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    const auto k1 = rhs.copy();

    // Step 2

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 0, make_array(dt), make_array(&k1));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto k2 = rhs.copy();

    // Step 3

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 0, make_array(1.0, dt / 4, dt / 4),
                         make_array(&old, &k1, &k2));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #3 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto &k3 = rhs;

    // Calculate new state vector
    statecomp_t::lincomb(var, 0, make_array(1.0, dt / 6, dt / 6, 2 * dt / 3),
                         make_array(&old, &k1, &k2, &k3));

  } else if (CCTK_EQUALS(method, "RK4")) {

    // k1 = f(y0)
    // k2 = f(y0 + h/2 k1)
    // k3 = f(y0 + h/2 k2)
    // k4 = f(y0 + h k3)
    // y1 = y0 + h/6 k1 + h/3 k2 + h/3 k3 + h/6 k4

    const auto old = var.copy();
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;

    // Step 1
    if (verbose)
      CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    const auto k1 = rhs.copy();

    // Step 2

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 1, make_array(dt / 2), make_array(&rhs));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto k2 = rhs.copy();

    // Step 3

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 0, make_array(1.0, dt / 2),
                         make_array(&old, &rhs));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #3 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto k3 = rhs.copy();

    // Step 4

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 0, make_array(1.0, dt), make_array(&old, &rhs));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #4 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto &k4 = rhs;

    // Calculate new state vector
    statecomp_t::lincomb(var, 0,
                         make_array(1.0, dt / 6, dt / 3, dt / 3, dt / 6),
                         make_array(&old, &k1, &k2, &k3, &k4));

  } else if (CCTK_EQUALS(method, "RKF78")) {

    typedef CCTK_REAL T;
    const auto R = [](T x, T y) { return x / y; };
    const tuple<vector<tuple<T, vector<T> > >, vector<T> > tableau{
        {
            {/* 1 */ 0, {}},                                           //
            {/* 2 */ R(2, 27), {R(2, 27)}},                            //
            {/* 3 */ R(1, 9), {R(1, 36), R(3, 36)}},                   //
            {/* 4 */ R(1, 6), {R(1, 24), 0, R(3, 24)}},                //
            {/* 5 */ R(5, 12), {R(20, 48), 0, R(-75, 48), R(75, 48)}}, //
            {/* 6 */ R(1, 2), {R(1, 20), 0, 0, R(5, 20), R(4, 20)}},   //
            {/* 7 */ R(5, 6),
             {R(-25, 108), 0, 0, R(125, 108), R(-260, 108), R(250, 108)}}, //
            {/* 8 */ R(1, 6),
             {R(31, 300), 0, 0, 0, R(61, 225), R(-2, 9), R(13, 900)}}, //
            {/* 9 */ R(2, 3),
             {2, 0, 0, R(-53, 6), R(704, 45), R(-107, 9), R(67, 90), 3}}, //
            {/* 10 */ R(1, 3),
             {R(-91, 108), 0, 0, R(23, 108), R(-976, 135), R(311, 54),
              R(-19, 60), R(17, 6), R(-1, 12)}}, //
            {/* 11 */ 1,
             {R(2383, 4100), 0, 0, R(-341, 164), R(4496, 1025), R(-301, 82),
              R(2133, 4100), R(45, 82), R(45, 164), R(18, 41)}}, //
                                                                 // {/* 12 */ 0,
            //  {R(3, 205), 0, 0, 0, 0, R(-6, 41), R(-3, 205), R(-3, 41), R(3,
            //  41),
            //   R(6, 41)}}, //
            // {/* 13 */ 1,
            //  {R(-1777, 4100), 0, 0, R(-341, 164), R(4496, 1025), R(-289, 82),
            //   R(2193, 4100), R(51, 82), R(33, 164), R(12, 41), 0, 1}}, //
        },
        {
            R(41, 840), 0, 0, 0, 0, R(34, 105), R(9, 35), R(9, 35), R(9, 280),
            R(9, 280), R(41, 840),
            // 0,
            // 0,
        }};

    // Check Butcher tableau
    const size_t nsteps = get<0>(tableau).size();
    {
      for (size_t step = 0; step < nsteps; ++step) {
        // TODO: Could allow <=
        assert(get<1>(get<0>(tableau).at(step)).size() == step);
        const auto &c = get<0>(get<0>(tableau).at(step));
        const auto &as = get<1>(get<0>(tableau).at(step));
        T x = 0;
        for (const auto &a : as)
          x += a;
        assert(fabs(x - c) <= 10 * numeric_limits<T>::epsilon());
      }
      // TODO: Could allow <=
      assert(get<1>(tableau).size() == nsteps);
      const auto &bs = get<1>(tableau);
      T x = 0;
      for (const auto &b : bs)
        x += b;
      assert(fabs(x - 1) <= 10 * numeric_limits<T>::epsilon());
    }

    const auto old = var.copy();

    vector<statecomp_t> ks;
    ks.reserve(nsteps);
    for (size_t step = 0; step < nsteps; ++step) {
      const auto &c = get<0>(get<0>(tableau).at(step));
      const auto &as = get<1>(get<0>(tableau).at(step));
      // Set current time
      *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + c * dt;
      // Add scaled RHS to state vector
      vector<CCTK_REAL> factors;
      vector<const statecomp_t *> srcs;
      factors.reserve(as.size() + 1);
      srcs.reserve(as.size() + 1);
      factors.push_back(1);
      srcs.push_back(&old);
      for (size_t i = 0; i < as.size(); ++i) {
        if (as.at(i) != 0) {
          factors.push_back(as.at(i) * dt);
          srcs.push_back(&ks.at(i));
        }
      }
      statecomp_t::lincomb(var, 0, factors, srcs);
      // TODO: Deallocate ks that are not needed any more
      CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
      if (verbose)
        CCTK_VINFO("Calculating RHS #%d at t=%g", int(step + 1),
                   double(cctkGH->cctk_time));
      CallScheduleGroup(cctkGH, "ODESolvers_RHS");
      ks.push_back(rhs.copy());
    }

    // Calculate new state vector
    const auto &bs = get<1>(tableau);
    vector<CCTK_REAL> factors;
    vector<const statecomp_t *> srcs;
    factors.reserve(bs.size() + 1);
    srcs.reserve(bs.size() + 1);
    factors.push_back(1);
    srcs.push_back(&old);
    for (size_t i = 0; i < bs.size(); ++i) {
      if (bs.at(i) != 0) {
        factors.push_back(bs.at(i) * dt);
        srcs.push_back(&ks.at(i));
      }
    }
    statecomp_t::lincomb(var, 0, factors, srcs);

  } else if (CCTK_EQUALS(method, "DP87")) {

    typedef CCTK_REAL T;
    const auto R = [](T x, T y) { return x / y; };
    // These coefficients are taken from the Einstein Toolkit, thorn
    // CactusNumerical/MoL, file RK87.c, written by Peter Diener,
    // following P. J. Prince and J. R. Dormand, Journal of
    // Computational and Applied Mathematics, volume 7, no 1, 1981
    const tuple<vector<vector<T> >, vector<T> > tableau{
        {
            {/*1*/},                                    //
            {/*2*/ R(1, 18)},                           //
            {/*3*/ R(1, 48), R(1, 16)},                 //
            {/*4*/ R(1, 32), 0, R(3, 32)},              //
            {/*5*/ R(5, 16), 0, -R(75, 64), R(75, 64)}, //
            {/*6*/ R(3, 80), 0, 0, R(3, 16), R(3, 20)}, //
            {/*7*/ R(29443841, 614563906), 0, 0, R(77736538, 692538347),
             -R(28693883, 1125000000), R(23124283, 1800000000)}, //
            {/*8*/ R(16016141, 946692911), 0, 0, R(61564180, 158732637),
             R(22789713, 633445777), R(545815736, 2771057229),
             -R(180193667, 1043307555)}, //
            {/*9*/ R(39632708, 573591083), 0, 0, -R(433636366, 683701615),
             -R(421739975, 2616292301), R(100302831, 723423059),
             R(790204164, 839813087), R(800635310, 3783071287)}, //
            {/*10*/ R(246121993, 1340847787), 0, 0,
             -R(37695042795, 15268766246), -R(309121744, 1061227803),
             -R(12992083, 490766935), R(6005943493, 2108947869),
             R(393006217, 1396673457), R(123872331, 1001029789)}, //
            {/*11*/ -R(1028468189, 846180014), 0, 0, R(8478235783, 508512852),
             R(1311729495, 1432422823), -R(10304129995, 1701304382),
             -R(48777925059, 3047939560), R(15336726248, 1032824649),
             -R(45442868181, 3398467696), R(3065993473, 597172653)}, //
            {/*12*/ R(185892177, 718116043), 0, 0, -R(3185094517, 667107341),
             -R(477755414, 1098053517), -R(703635378, 230739211),
             R(5731566787, 1027545527), R(5232866602, 850066563),
             -R(4093664535, 808688257), R(3962137247, 1805957418),
             R(65686358, 487910083)}, //
            {/*13*/ R(403863854, 491063109), 0, 0, -R(5068492393, 434740067),
             -R(411421997, 543043805), R(652783627, 914296604),
             R(11173962825, 925320556), -R(13158990841, 6184727034),
             R(3936647629, 1978049680), -R(160528059, 685178525),
             R(248638103, 1413531060), 0}, //
        },
        {R(14005451, 335480064), 0, 0, 0, 0, -R(59238493, 1068277825),
         R(181606767, 758867731), R(561292985, 797845732),
         -R(1041891430, 1371343529), R(760417239, 1151165299),
         R(118820643, 751138087), -R(528747749, 2220607170), R(1, 4)}};

    // Check Butcher tableau
    const size_t nsteps = get<0>(tableau).size();
    {
      for (size_t step = 0; step < nsteps; ++step)
        // TODO: Could allow <=
        assert(get<0>(tableau).at(step).size() == step);
      // TODO: Could allow <=
      assert(get<1>(tableau).size() == nsteps);
      const auto &bs = get<1>(tableau);
      T x = 0;
      for (const auto &b : bs)
        x += b;
      assert(fabs(x - 1) <= 10 * numeric_limits<T>::epsilon());
    }

    const auto old = var.copy();

    vector<statecomp_t> ks;
    ks.reserve(nsteps);
    for (size_t step = 0; step < nsteps; ++step) {
      const auto &as = get<0>(tableau).at(step);
      T c = 0;
      for (const auto &a : as)
        c += a;
      // Set current time
      *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + c * dt;
      // Add scaled RHS to state vector
      vector<CCTK_REAL> factors;
      vector<const statecomp_t *> srcs;
      factors.reserve(as.size() + 1);
      srcs.reserve(as.size() + 1);
      factors.push_back(1);
      srcs.push_back(&old);
      for (size_t i = 0; i < as.size(); ++i) {
        if (as.at(i) != 0) {
          factors.push_back(as.at(i) * dt);
          srcs.push_back(&ks.at(i));
        }
      }
      statecomp_t::lincomb(var, 0, factors, srcs);
      // TODO: Deallocate ks that are not needed any more
      CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
      if (verbose)
        CCTK_VINFO("Calculating RHS #%d at t=%g", int(step + 1),
                   double(cctkGH->cctk_time));
      CallScheduleGroup(cctkGH, "ODESolvers_RHS");
      ks.push_back(rhs.copy());
    }

    // Calculate new state vector
    const auto &bs = get<1>(tableau);
    vector<CCTK_REAL> factors;
    vector<const statecomp_t *> srcs;
    factors.reserve(bs.size() + 1);
    srcs.reserve(bs.size() + 1);
    factors.push_back(1);
    srcs.push_back(&old);
    for (size_t i = 0; i < bs.size(); ++i) {
      if (bs.at(i) != 0) {
        factors.push_back(bs.at(i) * dt);
        srcs.push_back(&ks.at(i));
      }
    }
    statecomp_t::lincomb(var, 0, factors, srcs);

  } else if (CCTK_EQUALS(method, "Implicit Euler")) {

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

    const auto y0 = var.copy();

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
    if (verbose)
      CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    const auto k1 = rhs.copy(); 

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    statecomp_t::lincomb(var, 1, make_array(dt / 2), make_array(&rhs));
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
    const auto y1 = var.copy();

    statecomp_t kprime2(cctkGH);
    statecomp_t::lincomb(kprime2, 0, make_array(-1.0, +1.0, -dt / 2),
                         make_array(&y0, &y1, &k1));

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt;
    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    const auto k2 = rhs.copy();

    statecomp_t::lincomb(var, 0, make_array(1.0, dt, dt),
                         make_array(&y0, &k2, &kprime2));

  } else if (CCTK_EQUALS(method, "IMEX_122")) {
    // code structure forked from "RK2"
    
    // function 'f' -> non-stiff part
    // k1_hat = f(y0)
    // Step 1 : 
    // beta = y0 + (dt/2)*k1_hat
    // alpha = dt/2
    // y1 -> get from 'UserSolvedFunction_G(beta)'
    // y1 = (dt/2)*g(y1) - beta

    // Step 2:
    // k2_hat = f(y1)
    // k1 = (y1-beta)/(dt/2)
    // y_np1 = y0 + dt*k2_hat + dt*k1

    // lincomb(dest, a, b^i, src^i)
    // dest = a * dest + b^i * src^i

    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
    const auto old_var = var.copy();
    CallScheduleGroup(cctkGH, "ODESolvers_NonStiffRHS");
    const auto old_non_stiff_rhs = rhs.copy();

    // Step 1
    statecomp_t::lincomb(var,0,make_array(1.0,dt/2),make_array(&old_var, &old_non_stiff_rhs));     // here 'var' = beta
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt/2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    const auto beta = var.copy();
    ODESolvers_alpha = dt/2;
    CallScheduleGroup(cctkGH, "ODESolvers_ImplicitStep");
    // time 
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
    const auto y1_var = var.copy(); //copy y1

    //Step 2
    CallScheduleGroup(cctkGH, "ODESolvers_NonStiffRHS");
    const auto y1_non_stiff_rhs = rhs.copy(); //calculates f(y1)
    // y_np1 = y0 + dt*f(y1) + 2*(y1-beta)
    statecomp_t::lincomb(y1_var,1.0,make_array(-1.0),make_array(&beta)); // compute (y1 - beta) and store it in 'y1_var'
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
    statecomp_t::lincomb(var,0,make_array(1.0,dt,2.0),make_array(&old_var,&y1_non_stiff_rhs,&y1_var));
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt/2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    // reset cctk_time and delta time to original
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_delta_time) = dt;
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;

    //ToDo: Add verbose for each step

  }else if (CCTK_EQUALS(method, "IMEX_ARS343")) {
    /* 
    lincomb(dest, a, b^i, src^i)
     dest = a * dest + b^i * src^i 
    */



    // function 'f' -> non-stiff RHS
    // cctk::delta_time for alpha
    // *const_cast<CCTK_REAL *>(&cctkGH->cctk_delta_time) = dt;
    //--------------Stage 1----------------
    // k1_hat = f(y0)
    // alpha = dt*a(1,1) = dt*Gamma
    // beta_1 = y0 + dt*(a_hat(2,1)*k1_hat)
    // y1 -> get from 'UserSolvedFunction_G(alpha,beta_1)' call
    // compute g(y1) -> k1 = (y1-beta_1)/alpha )
    //--------------Stage 2----------------
    // k2_hat = f(y1)
    // alpha = dt*a(2,2) = dt*Gamma
    // beta_2 = y0 + dt*( a(2,1)*k1 +
    //                    a_hat(3,1)*k1_hat +
    //                    a_hat(3,2)*k2_hat )
    // y2 -> get from 'UserSolvedFunction_G(alpha,beta_2)' call
    // compute g(y2) -> k2 = (y2-beta_2)/alpha )
    //--------------Stage 3----------------
    // k3_hat = f(y2)
    // alpha = dt*a(3,3) = dt*Gamma
    // beta_3 = y0 + dt*( a(3,1)*k1 +
    //                    a(3,2)*k2 +
    //                    a_hat(4,1)*k1_hat +
    //                    a_hat(4,2)*k2_hat +
    //                    a_hat(4,3)*k3_hat )
    // y3 = -> get from 'UserSolvedFunction_G(alpha,beta_3)' call
    // compute g(y3) -> k3 = (y3-beta_3)/alpha )
    //-----Calculate new state vector------
    // k4_hat = f(y3)
    //
    // y_np1 = y0 + dt*( b(1)*k1 +
    //                   b(2)*k2 + 
    //                   b(3)*k3 +
    //                   b_hat(1)*k1_hat +
    //                   b_hat(2)*k2_hat +
    //                   b_hat(3)*k3_hat + 
    //                   b_hat(4)*k4_hat )


    // Define butcher table data:
    const double Gamma = 0.4358665215;
    const double delta = -0.644373171;
    const double eta = 0.3966543747;
    const double mu = 0.5529291479;

    const auto implicit_butcher_table_a = [](int i,int j, const auto Gamma, const auto delta, const auto eta, const auto mu) {
      const double a_imp[4][4] = {{0,0,0,0},
                           {0,Gamma,0,0},
                           {0,(1 - Gamma)/2,Gamma,0},
                           {0,(1-Gamma-delta),delta,Gamma}};
      return a_imp[i][j]; };

    const auto explicit_butcher_table_a_hat = [](int i,int j, const auto Gamma, const auto delta, const auto eta, const auto mu) {
      i = i-1;
      j = j-1;
      const double a_exp[4][4] = {{0,0,0,0},
                           {Gamma,0,0,0},
                           {((1 + Gamma)/2) - eta,eta,0,0},
                           {(1-2*mu),mu,mu,0}};
      return a_exp[i][j]; };

    const auto implicit_b = [](int i, const auto Gamma, const auto delta, const auto eta, const auto mu) {
      const double b[4] = {0,(1-Gamma-delta),delta,Gamma};
      return b[i]; };

    const auto explicit_b_hat = [](int i, const auto Gamma, const auto delta, const auto eta, const auto mu) {
      i = i - 1;
      const double b_hat[4] = {0,(1-Gamma-delta),delta,Gamma};
      return b_hat[i]; };

    const auto implicit_c = [](int i, const auto Gamma, const auto delta, const auto eta, const auto mu) {
      const double c[4] = {0,Gamma,((1+Gamma)/2),1};
      return c[i]; };

    const auto explicit_c_hat = [](int i, const auto Gamma, const auto delta, const auto eta, const auto mu) {
      i = i - 1;
      const double c_hat[4] = {0,Gamma,((1+Gamma)/2),1};
      return c_hat[i]; };


    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
    //-----------------Stage 1---------------------
    const auto y0_var = var.copy();
    CallScheduleGroup(cctkGH, "ODESolvers_NonStiffRHS");
    const auto k1_hat = rhs.copy();

    ODESolvers_alpha = dt*Gamma;
    statecomp_t beta_product(cctkGH);
    beta_product = var.copy(); // copy dummy data into new state vector since lincomb requires new state vector data shape
    statecomp_t::lincomb(beta_product,0,make_array(explicit_butcher_table_a_hat(2,1,Gamma,delta,eta,mu)),make_array(&k1_hat));
    statecomp_t::lincomb(var,0,make_array(1.0,dt),make_array(&y0_var, &beta_product));     // here 'var' = beta
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    auto beta = var.copy();
    
    CallScheduleGroup(cctkGH, "ODESolvers_ImplicitStep");
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
    const auto y1_var = var.copy();
    statecomp_t k1(cctkGH);
    k1 = var.copy(); // copy dummy data into new state vector
    statecomp_t::lincomb(k1,0,make_array(1.0/ ODESolvers_alpha,-1.0/ ODESolvers_alpha),make_array(&y1_var,&beta));
    //-----------------Stage 2---------------------

    CallScheduleGroup(cctkGH, "ODESolvers_NonStiffRHS");
    const auto k2_hat = rhs.copy();

    statecomp_t::lincomb(beta_product,0,
    make_array(
      implicit_butcher_table_a(2,1,Gamma,delta,eta,mu),
      explicit_butcher_table_a_hat(3,1,Gamma,delta,eta,mu),
      explicit_butcher_table_a_hat(3,2,Gamma,delta,eta,mu)),
    make_array(&k1,&k1_hat,&k2_hat));
    statecomp_t::lincomb(var,0,make_array(1.0,dt),make_array(&y0_var, &beta_product));     // here 'var' = beta
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    beta = var.copy();

    CallScheduleGroup(cctkGH, "ODESolvers_ImplicitStep");
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
    //compute g(y2) -> k2:
    const auto y2_var = var.copy(); 
    statecomp_t k2(cctkGH);
    k2 = var.copy(); // copy dummy data into new state vector
    statecomp_t::lincomb(k2,0,make_array(1.0/ ODESolvers_alpha,-1.0/ ODESolvers_alpha),make_array(&y2_var,&beta));
    //-----------------Stage 3---------------------

    CallScheduleGroup(cctkGH, "ODESolvers_NonStiffRHS");
    const auto k3_hat = rhs.copy();

    statecomp_t::lincomb(beta_product,0,
      make_array(
        implicit_butcher_table_a(3,1,Gamma,delta,eta,mu),
        implicit_butcher_table_a(3,2,Gamma,delta,eta,mu),
        explicit_butcher_table_a_hat(4,1,Gamma,delta,eta,mu),
        explicit_butcher_table_a_hat(4,2,Gamma,delta,eta,mu),
        explicit_butcher_table_a_hat(4,3,Gamma,delta,eta,mu)
      ),
    make_array(&k1,&k2,&k1_hat,&k2_hat,&k3_hat));
    statecomp_t::lincomb(var,0,make_array(1.0,dt),make_array(&y0_var, &beta_product));     // here 'var' = beta
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    beta = var.copy();
    CallScheduleGroup(cctkGH, "ODESolvers_ImplicitStep");
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
    const auto y3_var = var.copy();
    statecomp_t k3(cctkGH);
    k3 = var.copy(); // copy dummy data into new state vector
    statecomp_t::lincomb(k3,0,make_array(1.0/ ODESolvers_alpha,-1.0/ ODESolvers_alpha),make_array(&y3_var,&beta));
    //-----------Calculate new state vector-------------

    CallScheduleGroup(cctkGH, "ODESolvers_NonStiffRHS");
    const auto k4_hat = rhs.copy();
    
    statecomp_t ynp1_product(cctkGH);
    ynp1_product = var.copy(); // copy dummy data into new state vector
    statecomp_t::lincomb(ynp1_product,0,
      make_array(
        implicit_b(1,Gamma,delta,eta,mu),
        implicit_b(2,Gamma,delta,eta,mu),
        implicit_b(3,Gamma,delta,eta,mu),
        explicit_b_hat(1,Gamma,delta,eta,mu),
        explicit_b_hat(2,Gamma,delta,eta,mu),
        explicit_b_hat(3,Gamma,delta,eta,mu),
        explicit_b_hat(4,Gamma,delta,eta,mu)
      ),
      make_array(&k1,&k2,&k3,&k1_hat,&k2_hat,&k3_hat,&k4_hat));
    statecomp_t::lincomb(var,0,make_array(1.0,dt),make_array(&y0_var,&ynp1_product));
    

    // reset cctk_time and delta time to original
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_delta_time) = dt;
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
    
    //ToDo: Add verbose for each step

  } else {
    assert(0);
  }

  // Reset current time
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = saved_time;
  // Apply last boundary conditions
  CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
  if (verbose)
    CCTK_VINFO("Calculated new state at t=%g", double(cctkGH->cctk_time));

  // TODO: Update time here, and not during time level cycling in the driver
}

} // namespace ODESolvers
