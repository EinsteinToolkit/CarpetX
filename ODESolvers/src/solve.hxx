#ifndef CARPETX_ODESOLVERS_SOLVE_HXX
#define CARPETX_ODESOLVERS_SOLVE_HXX

// TODO: Don't include files from other thorns; create a proper interface
#include "../../CarpetX/src/driver.hxx"
#include "../../CarpetX/src/schedule.hxx"
#include "../../CarpetX/src/timer.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>
#include <util_Table.h>

#include <div.hxx>

#include <AMReX_MultiFab.H>

#if defined _OPENMP || defined __HIPCC__
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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

  statecomp_t() = default;

  statecomp_t(statecomp_t &&) = default;
  statecomp_t &operator=(statecomp_t &&) = default;

  // Don't allow copies because we might own stuff
  statecomp_t(const statecomp_t &) = delete;
  statecomp_t &operator=(const statecomp_t &) = delete;

  vector<GHExt::PatchData::LevelData::GroupData *> groupdatas;
  vector<amrex::MultiFab *> mfabs;
  // Timelevel of every mfab/groupdata entry. Subcycling aliases the previous
  // step at tl=1; all other paths use tl=0.
  int timelevel = 0;

  static void init_tmp_mfabs();
  static void free_tmp_mfabs();

  void set_valid(const valid_t valid) const;
  template <size_t N>
  static void combine_valids(const statecomp_t &dst, const CCTK_REAL scale,
                             const array<CCTK_REAL, N> &factors,
                             const array<const statecomp_t *, N> &srcs,
                             const valid_t where);
  void check_valid(const valid_t required, const function<string()> &why) const;
  void check_valid(const valid_t required, const string &why) const {
    check_valid(required, [=]() { return why; });
  }

  statecomp_t copy(const valid_t where) const;

  template <size_t N>
  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const array<CCTK_REAL, N> &factors,
                      const array<const statecomp_t *, N> &srcs,
                      const valid_t where);
  template <size_t N>
  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const array<CCTK_REAL, N> &factors,
                      const array<statecomp_t *, N> &srcs,
                      const valid_t where) {
    array<const statecomp_t *, N> srcs1;
    for (size_t n = 0; n < N; ++n)
      srcs1[n] = srcs[n];
    lincomb(dst, scale, factors, srcs1, where);
  }

  static void lincomb(const statecomp_t &dst, CCTK_REAL scale,
                      const vector<CCTK_REAL> &factors,
                      const vector<const statecomp_t *> &srcs,
                      const valid_t where);
};

template <std::size_t N> using reals = std::array<CCTK_REAL, N>;
template <std::size_t N> using states = std::array<const statecomp_t *, N>;

////////////////////////////////////////////////////////////////////////////////

inline int groupindex(const int other_gi, std::string gn) {
  // If the group name does not contain a colon, then prefix the current group's
  // implementation or thorn name
  if (gn.find(':') == std::string::npos) {
    const char *const thorn_or_impl = CCTK_GroupImplementationI(other_gi);
    assert(thorn_or_impl);
    const char *const impl = CCTK_ThornImplementation(thorn_or_impl);
    const char *const thorn = CCTK_ImplementationThorn(thorn_or_impl);
    assert(impl || thorn);
    const char *prefix;
    if (!impl) {
      prefix = thorn;
    } else if (!thorn) {
      prefix = impl;
    } else {
      assert(strcmp(impl, thorn) == 0);
      prefix = impl;
    }
    std::ostringstream buf;
    buf << prefix << "::" + gn;
    gn = buf.str();
  }
  const int gi = CCTK_GroupIndex(gn.c_str());
  return gi;
}

inline int get_group_rhs(const int gi) {
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  std::vector<char> rhs_buf(1000);
  const int iret =
      Util_TableGetString(tags, rhs_buf.size(), rhs_buf.data(), "rhs");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    rhs_buf[0] = '\0'; // default: empty (no RHS)
  } else if (iret >= 0) {
    // do nothing
  } else {
    assert(0);
  }

  const std::string str(rhs_buf.data());
  if (str.empty())
    return -1; // No RHS specified

  const int rhs = groupindex(gi, str);
  if (rhs < 0)
    CCTK_VERROR("Variable group \"%s\" declares a RHS group \"%s\". "
                "That group does not exist.",
                CCTK_FullGroupName(gi), str.c_str());
  assert(rhs != gi);

  return rhs;
}

inline std::vector<int> get_group_dependents(const int gi) {
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  std::vector<char> dependents_buf(1000);
  const int iret = Util_TableGetString(tags, dependents_buf.size(),
                                       dependents_buf.data(), "dependents");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    dependents_buf[0] = '\0'; // default: empty (no DEPENDENTS)
  } else if (iret >= 0) {
    // do nothing
  } else {
    assert(0);
  }

  std::vector<int> dependents;
  const std::string str(dependents_buf.data());
  std::size_t pos = 0;
  for (;;) {
    // Skip white space
    while (pos < str.size() && std::isspace(str[pos]))
      ++pos;
    if (pos == str.size())
      break;
    // Read group name
    const std::size_t group_begin = pos;
    while (pos < str.size() && !std::isspace(str[pos]))
      ++pos;
    const std::size_t group_end = pos;
    const std::string groupname =
        str.substr(group_begin, group_end - group_begin);
    const int dep_gi = groupindex(gi, groupname);
    if (dep_gi < 0)
      CCTK_VERROR("Variable group \"%s\" declares a dependent group \"%s\". "
                  "That group does not exist.",
                  CCTK_FullGroupName(gi), groupname.c_str());
    dependents.push_back(dep_gi);
  }

  return dependents;
}

// Mark groups as invalid
inline void mark_invalid(const std::vector<int> &groups) {
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const int gi : groups) {
      auto &groupdata = *leveldata.groupdata.at(gi);
      // Invalidate all variables of the current time level
      const int tl = 0;
      for (auto &why_valid : groupdata.valid.at(tl))
        why_valid =
            why_valid_t([] { return "ODESolvers updated the state vector"; });
    }
  });
}

} // namespace ODESolvers

#endif // #ifndef CARPETX_ODESOLVERS_SOLVE_HXX
