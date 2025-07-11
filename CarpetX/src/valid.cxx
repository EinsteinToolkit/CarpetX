#include "driver.hxx"
#include "schedule.hxx"
#include "timer.hxx"
#include "valid.hxx"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace CarpetX {

////////////////////////////////////////////////////////////////////////////////

std::string valid_t::explanation() const {
  const auto valstr = [](bool v) { return v ? "valid" : "invalid"; };
  std::ostringstream buf;
  buf << "\n"
      << "  The interior is " << valstr(valid_int) << ".\n"
      << "  The outer boundary is " << valstr(valid_outer) << ".\n"
      << "  The ghost zones are " << valstr(valid_ghosts) << ".\n";
  return buf.str();
}

std::ostream &operator<<(std::ostream &os, const valid_t v) {
  auto str = [](bool v) { return v ? "VAL" : "INV"; };
  return os << "[int:" << str(v.valid_int) << ",outer:" << str(v.valid_outer)
            << ",ghosts:" << str(v.valid_ghosts) << "]";
}

valid_t::operator string() const {
  std::ostringstream buf;
  buf << *this;
  return buf.str();
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const valid_t v) {
  yaml << YAML::LocalTag("valid-1.0.0");
  yaml << YAML::Flow << YAML::BeginMap;
  yaml << YAML::Key << "int" << YAML::Value << v.valid_int;
  yaml << YAML::Key << "outer" << YAML::Value << v.valid_outer;
  yaml << YAML::Key << "ghosts" << YAML::Value << v.valid_ghosts;
  yaml << YAML::EndMap;
  return yaml;
}

std::string why_valid_t::explanation() const {
  const auto valstr = [](bool v) { return v ? "valid" : "invalid"; };
  std::ostringstream buf;
  buf << "\n"
      << "  The interior is " << valstr(valid.valid_int)
      << " because: " << why_int() << ".\n"
      << "  The outer boundary is " << valstr(valid.valid_outer)
      << " because: " << why_outer() << ".\n"
      << "  The ghost zones are " << valstr(valid.valid_ghosts)
      << " because: " << why_ghosts() << ".\n";
  return buf.str();
}

std::ostream &operator<<(std::ostream &os, const why_valid_t &why) {
  return os << why.valid << ","
            << "why{int:" << why.why_int() << ","
            << "outer:" << why.why_outer() << ","
            << "ghosts:" << why.why_ghosts() << "}";
}

why_valid_t::operator string() const {
  std::ostringstream buf;
  buf << *this;
  return buf.str();
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const why_valid_t &why) {
  yaml << YAML::LocalTag("why_valid-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "valid" << YAML::Value << why.valid;
  yaml << YAML::Key << "why" << YAML::Value << YAML::BeginMap;
  yaml << YAML::Key << "int" << YAML::Value << why.why_int();
  yaml << YAML::Key << "outer" << YAML::Value << why.why_outer();
  yaml << YAML::Key << "ghosts" << YAML::Value << why.why_ghosts();
  yaml << YAML::EndMap;
  yaml << YAML::EndMap;
  return yaml;
}

////////////////////////////////////////////////////////////////////////////////

// valid/invalid flags

// Ensure grid functions are valid
void error_if_invalid(const GHExt::PatchData::LevelData::GroupData &groupdata,
                      int vi, int tl, const valid_t &required,
                      const std::function<std::string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VERROR("%s: Grid function \"%s\" is invalid on patch %d, refinement "
                "level %d, time level %d; required: %s, found: %s",
                msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
                groupdata.patch, groupdata.level, tl,
                required.explanation().c_str(),
                groupdata.valid.at(tl).at(vi).explanation().c_str());
}
void warn_if_invalid(const GHExt::PatchData::LevelData::GroupData &groupdata,
                     int vi, int tl, const valid_t &required,
                     const std::function<std::string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VWARN(CCTK_WARN_ALERT,
               "%s: Grid function \"%s\" is invalid on patch %d, refinement "
               "level %d, time level %d; required: %s, found: %s",
               msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
               groupdata.patch, groupdata.level, tl,
               required.explanation().c_str(),
               groupdata.valid.at(tl).at(vi).explanation().c_str());
}

// Ensure arrays are valid
void error_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata,
                      int vi, int tl, const valid_t &required,
                      const std::function<std::string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VERROR(
        "%s: Array \"%s\" is invalid on time level %d; required %s, found %s",
        msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi), tl,
        required.explanation().c_str(),
        groupdata.valid.at(tl).at(vi).explanation().c_str());
}
void warn_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                     int tl, const valid_t &required,
                     const std::function<std::string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VWARN(
        CCTK_WARN_ALERT,
        "%s: Array \"%s\" is invalid on time level %d; required %s, found %s",
        msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi), tl,
        required.explanation().c_str(),
        groupdata.valid.at(tl).at(vi).explanation().c_str());
}

////////////////////////////////////////////////////////////////////////////////

// Poison values to catch uninitialized variables

// Poison grid functions
void poison_invalid_gf(const active_levels_t &active_levels, const int gi,
                       const int vi, const int tl) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  static Timer timer("poison_invalid<GF>");
  Interval interval(timer);

  const poison_value_t<CCTK_REAL> poison_value;
  CCTK_REAL poison;
  poison_value.set_to_poison(poison);

  active_levels.loop_parallel([&](const int patch, const int level,
                                  const int index, const int component,
                                  const cGH *restrict const cctkGH) {
    const auto &patchdata = ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    auto &restrict groupdata = *leveldata.groupdata.at(gi);

    const valid_t &valid = groupdata.valid.at(tl).at(vi).get();
    if (valid.valid_all())
      return;

    const Loop::GridDescBaseDevice grid(cctkGH);
    const Loop::GF3D2layout layout(cctkGH, groupdata.indextype);
    const Loop::GF3D2<CCTK_REAL> gf(
        layout, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, tl, groupdata.firstvarindex + vi)));

    if (!valid.valid_any()) {
      grid.loop_device_idx<where_t::everywhere>(
          groupdata.indextype, groupdata.nghostzones,
          [=] CCTK_DEVICE(const Loop::PointDesc &p)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { gf(p.I) = poison; });
    } else {
      if (!valid.valid_int)
        grid.loop_device_idx<where_t::interior>(
            groupdata.indextype, groupdata.nghostzones,
            [=] CCTK_DEVICE(const Loop::PointDesc &p)
                CCTK_ATTRIBUTE_ALWAYS_INLINE { gf(p.I) = poison; });
      if (!valid.valid_outer)
        grid.loop_device_idx<where_t::boundary>(
            groupdata.indextype, groupdata.nghostzones,
            [=] CCTK_DEVICE(const Loop::PointDesc &p)
                CCTK_ATTRIBUTE_ALWAYS_INLINE { gf(p.I) = poison; });
      if (!valid.valid_ghosts)
        grid.loop_device_idx<where_t::ghosts>(
            groupdata.indextype, groupdata.nghostzones,
            [=] CCTK_DEVICE(const Loop::PointDesc &p)
                CCTK_ATTRIBUTE_ALWAYS_INLINE { gf(p.I) = poison; });
    }
  });

  // Synchronize because we access GPU memory on the CPU
  synchronize();
}

// Poison arrays
void poison_invalid_ga(const int gi, const int vi, const int tl) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  static Timer timer("poison_invalid<GA>");
  Interval interval(timer);

  auto &restrict globaldata = ghext->globaldata;
  auto &restrict arraygroupdata = *globaldata.arraygroupdata.at(gi);
  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);

  const valid_t &valid = arraygroupdata.valid.at(tl).at(vi).get();
  if (valid.valid_all())
    return;

  if (!valid.valid_int) {
    int dimension = arraygroupdata.dimension;
    const int *gsh = arraygroupdata.gsh;
    int n_elems = 1;
    for (int i = 0; i < dimension; i++)
      n_elems *= gsh[i];
    // TODO: use AnyScalarTypeRef for this?
    assert(group.vartype == CCTK_VARIABLE_COMPLEX ||
           group.vartype == CCTK_VARIABLE_REAL ||
           group.vartype == CCTK_VARIABLE_INT);
    switch (group.vartype) {
    case CCTK_VARIABLE_COMPLEX: {
      CCTK_COMPLEX *restrict ptr = static_cast<CCTK_COMPLEX *>(
          arraygroupdata.data.at(tl).data_at(vi * n_elems));
      const poison_value_t<CCTK_COMPLEX> poison_value;
      poison_value.set_to_poison(ptr, n_elems);
    } break;
    case CCTK_VARIABLE_REAL: {
      CCTK_REAL *restrict ptr = static_cast<CCTK_REAL *>(
          arraygroupdata.data.at(tl).data_at(vi * n_elems));
      const poison_value_t<CCTK_REAL> poison_value;
      poison_value.set_to_poison(ptr, n_elems);
    } break;
    case CCTK_VARIABLE_INT: {
      CCTK_INT *restrict ptr = static_cast<CCTK_INT *>(
          arraygroupdata.data.at(tl).data_at(vi * n_elems));
      const poison_value_t<CCTK_INT> poison_value;
      poison_value.set_to_poison(ptr, n_elems);
    } break;
    }
  }
}

// Ensure grid functions are not poisoned
void check_valid_gf(const active_levels_t &active_levels, const int gi,
                    const int vi, const int tl,
                    const nan_handling_t nan_handling1,
                    const std::function<std::string()> &msg) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  static Timer timer("check_valid<GF>");
  Interval interval(timer);

#warning "TODO"
  constexpr nan_handling_t nan_handling = nan_handling_t::forbid_nans;

  const auto is_poison = [] CCTK_DEVICE CCTK_HOST(
                             const CCTK_REAL val) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    poison_value_t<CCTK_REAL> const poison_value;
    if (poison_value.is_poison(val))
      return true;
    using std::isnan;
    if (nan_handling != nan_handling_t::allow_nans && isnan(val))
      return true;
    return false;
  };

  amrex::FArrayBox poison_found(
      amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0)), 1,
      amrex::The_Async_Arena());
#ifdef AMREX_USE_GPU
  constexpr auto run_on = amrex::RunOn::Device;
#else
  constexpr auto run_on = amrex::RunOn::Host;
#endif
  poison_found.operator= <run_on>(0.0);
  CCTK_REAL *restrict const poison_found_ptr = poison_found.dataPtr();

  active_levels.loop_parallel([&](const int patch, const int level,
                                  const int index, const int component,
                                  const cGH *restrict const cctkGH) {
    const auto &patchdata = ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    auto &restrict groupdata = *leveldata.groupdata.at(gi);

    const valid_t &valid = groupdata.valid.at(tl).at(vi).get();
    if (!valid.valid_any())
      return;

    const Loop::GridDescBaseDevice grid(cctkGH);
    const Loop::GF3D2layout layout(cctkGH, groupdata.indextype);
    const Loop::GF3D2<const CCTK_REAL> gf(
        layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, tl, groupdata.firstvarindex + vi)));

    const auto update_poison_found =
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (CCTK_BUILTIN_EXPECT(!is_poison(gf(p.I)), true))
            return;
#pragma omp atomic write
          *poison_found_ptr = 0.0 / 0.0;
        };

    if (valid.valid_all()) {
      grid.loop_device_idx<where_t::everywhere>(
          groupdata.indextype, groupdata.nghostzones, update_poison_found);
    } else {
      if (valid.valid_int)
        grid.loop_device_idx<where_t::interior>(
            groupdata.indextype, groupdata.nghostzones, update_poison_found);
      if (valid.valid_outer)
        grid.loop_device_idx<where_t::boundary>(
            groupdata.indextype, groupdata.nghostzones, update_poison_found);
      if (valid.valid_ghosts)
        grid.loop_device_idx<where_t::ghosts>(
            groupdata.indextype, groupdata.nghostzones, update_poison_found);
    }
  });
  synchronize();

  if (!poison_found.contains_nan<run_on>())
    return;

  std::size_t nan_count{0};
  std::array<int, 3> nan_imin, nan_imax;
  std::array<CCTK_REAL, 3> nan_xmin, nan_xmax;
  for (int d = 0; d < 3; ++d) {
    nan_imin[d] = std::numeric_limits<int>::max();
    nan_imax[d] = std::numeric_limits<int>::min();
    nan_xmin[d] = +1.0 / 0.0;
    nan_xmax[d] = -1.0 / 0.0;
  }

  struct info_t {
    where_t where;
    int patch, level, component;
    vect<int, dim> I;
    vect<CCTK_REAL, dim> X;
    CCTK_REAL val;
  };
  std::vector<info_t> infos;

  active_levels.loop_serially([&](const int patch, const int level,
                                  const int index, const int component,
                                  const cGH *restrict const cctkGH) {
    const auto &patchdata = ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    auto &restrict groupdata = *leveldata.groupdata.at(gi);

    const valid_t &valid = groupdata.valid.at(tl).at(vi).get();
    if (!valid.valid_any())
      return;

    const Loop::GridDescBaseDevice grid(cctkGH);
    const Loop::GF3D2layout layout(cctkGH, groupdata.indextype);
    const Loop::GF3D2<const CCTK_REAL> gf(
        layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, tl, groupdata.firstvarindex + vi)));

    const auto update_nan_count = [&](const Loop::PointDesc &p,
                                      const where_t where) {
      if (CCTK_BUILTIN_EXPECT(!is_poison(gf(p.I)), true))
        return;

      ++nan_count;
      nan_imin[0] = min(nan_imin[0], grid.lbnd[0] + p.i);
      nan_imin[1] = min(nan_imin[1], grid.lbnd[1] + p.j);
      nan_imin[2] = min(nan_imin[2], grid.lbnd[2] + p.k);
      nan_imax[0] = max(nan_imax[0], grid.lbnd[0] + p.i);
      nan_imax[1] = max(nan_imax[1], grid.lbnd[1] + p.j);
      nan_imax[2] = max(nan_imax[2], grid.lbnd[2] + p.k);
      nan_xmin[0] = fmin(nan_xmin[0], p.x);
      nan_xmin[1] = fmin(nan_xmin[1], p.y);
      nan_xmin[2] = fmin(nan_xmin[2], p.z);
      nan_xmax[0] = fmax(nan_xmax[0], p.x);
      nan_xmax[1] = fmax(nan_xmax[1], p.y);
      nan_xmax[2] = fmax(nan_xmax[2], p.z);

      infos.push_back(
          info_t{where, p.patch, p.level, p.component, p.I, p.X, gf(p.I)});
    };

    if (valid.valid_all()) {
      grid.loop_idx(where_t::everywhere, groupdata.indextype,
                    groupdata.nghostzones, [&](const Loop::PointDesc &p) {
                      update_nan_count(p, where_t::everywhere);
                    });
    } else {
      if (valid.valid_int)
        grid.loop_idx(where_t::interior, groupdata.indextype,
                      groupdata.nghostzones, [&](const Loop::PointDesc &p) {
                        update_nan_count(p, where_t::interior);
                      });
      if (valid.valid_outer)
        grid.loop_idx(where_t::boundary, groupdata.indextype,
                      groupdata.nghostzones, [&](const Loop::PointDesc &p) {
                        update_nan_count(p, where_t::boundary);
                      });
      if (valid.valid_ghosts)
        grid.loop_idx(where_t::ghosts, groupdata.indextype,
                      groupdata.nghostzones, [&](const Loop::PointDesc &p) {
                        update_nan_count(p, where_t::ghosts);
                      });
    }
  });

  const auto &patchdata0 = ghext->patchdata.at(0);
  const auto &leveldata0 = patchdata0.leveldata.at(0);
  const auto &groupdata0 = *leveldata0.groupdata.at(gi);
  CCTK_VWARN(CCTK_WARN_ALERT,
             "%s: Grid function \"%s\" contains %td nans, infinities, or "
             "poison in box (%g,%g,%g):(%g,%g,%g); expected valid %s",
             msg().c_str(), CCTK_FullVarName(groupdata0.firstvarindex + vi),
             std::size_t(nan_count), double(nan_xmin[0]), double(nan_xmin[1]),
             double(nan_xmin[2]), double(nan_xmax[0]), double(nan_xmax[1]),
             double(nan_xmax[2]),
             groupdata0.valid.at(tl).at(vi).explanation().c_str());

  std::sort(infos.begin(), infos.end(), [](const info_t &a, const info_t &b) {
    if (a.level < b.level)
      return true;
    if (a.level > b.level)
      return false;
    if (a.patch < b.patch)
      return true;
    if (a.patch > b.patch)
      return false;
    if (a.component < b.component)
      return true;
    if (a.component > b.component)
      return false;
    const std::less<vect<int, dim> > lt;
    return lt(reversed(a.I), reversed(b.I));
  });

  std::ostringstream buf;
  buf << setprecision(std::numeric_limits<CCTK_REAL>::digits10 + 1);
  for (const auto &info : infos)
    buf << "\n"
        << info.where << " level " << info.level << " patch " << info.patch
        << " component " << info.component << " " << info.I << " " << info.X
        << " " << info.val;
  CCTK_WARN(CCTK_WARN_ALERT, buf.str().c_str());

  CCTK_VERROR("%s: Grid function \"%s\" contains nans, infinities, or poison; "
              "expected valid %s",
              msg().c_str(), CCTK_FullVarName(groupdata0.firstvarindex + vi),
              groupdata0.valid.at(tl).at(vi).explanation().c_str());
}

// Ensure arrays are not poisoned
void check_valid_ga(const int gi, const int vi, const int tl,
                    const nan_handling_t nan_handling1,
                    const std::function<std::string()> &msg) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  static Timer timer("check_valid<GA>");
  Interval interval(timer);

  auto &restrict globaldata = ghext->globaldata;
  auto &restrict arraygroupdata = *globaldata.arraygroupdata.at(gi);
  const valid_t &valid = arraygroupdata.valid.at(tl).at(vi).get();

  // arrays have no boundary so we expect them to alway be valid
  assert(valid.valid_outer && valid.valid_ghosts);
  if (!valid.valid_int)
    return;

  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);

#warning "TODO"
  using std::isnan;
  constexpr nan_handling_t nan_handling = nan_handling_t::forbid_nans;

  std::size_t nan_count{0};

  int dimension = arraygroupdata.dimension;
  const int *gsh = arraygroupdata.gsh;
  int n_elems = 1;
  for (int i = 0; i < dimension; i++)
    n_elems *= gsh[i];
  // TODO: use AnyScalarTypeRef for this?
  switch (group.vartype) {
  case CCTK_VARIABLE_COMPLEX: {
    const poison_value_t<CCTK_COMPLEX> poison_value;
    const CCTK_COMPLEX *restrict const ptr =
        static_cast<const CCTK_COMPLEX *const>(
            arraygroupdata.data.at(tl).data_at(vi * n_elems));
    for (int i = 0; i < n_elems; i++) {
      if (CCTK_BUILTIN_EXPECT(
              poison_value.is_poison(ptr[i]) ||
                  (nan_handling != nan_handling_t::allow_nans &&
                   (isnan(ptr[i].real()) || isnan(ptr[i].imag()))),
              false))
        ++nan_count;
    }
  } break;
  case CCTK_VARIABLE_REAL: {
    const poison_value_t<CCTK_REAL> poison_value;
    const CCTK_REAL *restrict const ptr = static_cast<const CCTK_REAL *const>(
        arraygroupdata.data.at(tl).data_at(vi * n_elems));
    for (int i = 0; i < n_elems; i++) {
      if (CCTK_BUILTIN_EXPECT(
              poison_value.is_poison(ptr[i]) ||
                  (nan_handling != nan_handling_t::allow_nans && isnan(ptr[i])),
              false))
        ++nan_count;
    }
  } break;
  case CCTK_VARIABLE_INT: {
    const poison_value_t<CCTK_INT> poison_value;
    const CCTK_INT *restrict const ptr = static_cast<const CCTK_INT *const>(
        arraygroupdata.data.at(tl).data_at(vi * n_elems));
    for (int i = 0; i < n_elems; i++) {
      if (CCTK_BUILTIN_EXPECT(poison_value.is_poison(ptr[i]), false))
        ++nan_count;
    }
  } break;
  }

  if (nan_count == 0)
    return;

  CCTK_VERROR("%s: Grid array \"%s\" has %td nans on time level %d; "
              "expected valid %s",
              msg().c_str(),
              CCTK_FullVarName(arraygroupdata.firstvarindex + vi), nan_count,
              tl, arraygroupdata.valid.at(tl).at(vi).explanation().c_str());
}

////////////////////////////////////////////////////////////////////////////////

// Checksums to catch illegal modifications

checksums_t calculate_checksums(
    const std::vector<std::vector<std::vector<valid_t> > > &will_write) {
  DECLARE_CCTK_PARAMETERS;

  checksums_t checksums;

  if (!poison_undefined_values)
    return checksums;

  static Timer timer("calculate_checksums");
  Interval interval(timer);

  assert(active_levels);
  active_levels->loop_parallel([&](const int patch, const int level,
                                   const int index, const int component,
                                   const cGH *restrict const cctkGH) {
    const auto &patchdata = ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    for (const auto &groupdataptr : leveldata.groupdata) {
      if (!groupdataptr)
        continue;
      const auto &restrict groupdata = *groupdataptr;

      const Loop::GridDescBaseDevice grid(cctkGH);
      const Loop::GF3D2layout layout(cctkGH, groupdata.indextype);

      for (int vi = 0; vi < groupdata.numvars; ++vi) {
        for (int tl = 0; tl < int(groupdata.valid.size()); ++tl) {
          const tiletag_t tiletag{patch, level, component, groupdata.groupindex,
                                  vi,    tl};

          const auto &valid = groupdata.valid.at(tl).at(vi).get();
          // No information given for this timelevel; assume not written
          if (tl >= int(will_write.at(groupdata.groupindex).at(vi).size()))
            continue;
          const auto &wr = will_write.at(groupdata.groupindex).at(vi).at(tl);
          valid_t to_check = valid & ~wr;

          // Check only those variables which are valid, and where
          // some part (but not everything) is written
          if (!(wr.valid_any() && to_check.valid_any()))
            continue;

          const Loop::GF3D2<const CCTK_REAL> gf(
              layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(
                          cctkGH, tl, groupdata.firstvarindex + vi)));

          checksum_t checksum(to_check);
          checksum.add(tiletag);
          const auto add_point = [&](const Loop::PointDesc &p) {
            checksum.add(gf(p.I));
          };

          if (to_check.valid_int)
            grid.loop_idx(where_t::interior, groupdata.indextype,
                          groupdata.nghostzones, add_point);

          if (to_check.valid_outer)
            grid.loop_idx(where_t::boundary, groupdata.indextype,
                          groupdata.nghostzones, add_point);

          if (to_check.valid_ghosts)
            grid.loop_idx(where_t::ghosts, groupdata.indextype,
                          groupdata.nghostzones, add_point);

#pragma omp critical(CarpetX_calculate_checksums)
          checksums[tiletag] = checksum;
        }
      }
    }
  });

  return checksums;
}

void check_checksums(const checksums_t &checksums,
                     const std::function<std::string()> &where) {
  DECLARE_CCTK_PARAMETERS;

  if (!poison_undefined_values)
    return;
  if (checksums.empty())
    return;

  static Timer timer("check_checksums");
  Interval interval(timer);

  assert(active_levels);
  active_levels->loop_parallel([&](const int patch, const int level,
                                   const int index, const int component,
                                   const cGH *restrict const cctkGH) {
    const auto &patchdata = ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    for (const auto &groupdataptr : leveldata.groupdata) {
      if (!groupdataptr)
        continue;
      const auto &restrict groupdata = *groupdataptr;

      const Loop::GridDescBaseDevice grid(cctkGH);
      const Loop::GF3D2layout layout(cctkGH, groupdata.indextype);

      for (int vi = 0; vi < groupdata.numvars; ++vi) {
        for (int tl = 0; tl < int(groupdata.valid.size()); ++tl) {
          const tiletag_t tiletag{patch, level, component, groupdata.groupindex,
                                  vi,    tl};

          if (!checksums.count(tiletag))
            continue;

          const auto &old_checksum = checksums.at(tiletag);
          const auto &did_check = old_checksum.where;
          assert(did_check.valid_any());

          const Loop::GF3D2<const CCTK_REAL> gf(
              layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(
                          cctkGH, tl, groupdata.firstvarindex + vi)));

          checksum_t checksum(did_check);
          checksum.add(tiletag);
          const auto add_point = [&](const Loop::PointDesc &p) {
            checksum.add(gf(p.I));
          };

          if (did_check.valid_int)
            grid.loop_idx(where_t::interior, groupdata.indextype,
                          groupdata.nghostzones, add_point);

          if (did_check.valid_outer)
            grid.loop_idx(where_t::boundary, groupdata.indextype,
                          groupdata.nghostzones, add_point);

          if (did_check.valid_ghosts)
            grid.loop_idx(where_t::ghosts, groupdata.indextype,
                          groupdata.nghostzones, add_point);

          if (checksum != old_checksum)
#pragma omp critical
            CCTK_VERROR(
                "%s: Checksum mismatch: variable %s, tile %s, "
                "int:%d,outer:%d,ghosts:%d, old checksum %s, new checksum %s",
                where().c_str(),
                CCTK_FullVarName(groupdata.firstvarindex + tiletag.vi),
                std::string(tiletag).c_str(), int(did_check.valid_int),
                int(did_check.valid_outer), int(did_check.valid_ghosts),
                std::string(old_checksum).c_str(),
                std::string(checksum).c_str());
        }
      }
    }
  });
}

} // namespace CarpetX
