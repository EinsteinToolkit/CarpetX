#include "driver.hxx"
#include "loop_device.hxx"
#include "schedule.hxx"
#include "timer.hxx"
#include "valid.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace CarpetX {
using namespace std;

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
  ostringstream buf;
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
  ostringstream buf;
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
                      const function<string()> &msg) {
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
                     const function<string()> &msg) {
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
                      const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VERROR("%s: Array \"%s\" is invalid on time level %d; "
                "required %s, found %s",
                msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
                tl, required.explanation().c_str(),
                groupdata.valid.at(tl).at(vi).explanation().c_str());
}
void warn_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                     int tl, const valid_t &required,
                     const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VWARN(CCTK_WARN_ALERT,
               "%s: Array \"%s\" is invalid on time level %d; "
               "required %s, found %s",
               msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
               tl, required.explanation().c_str(),
               groupdata.valid.at(tl).at(vi).explanation().c_str());
}

////////////////////////////////////////////////////////////////////////////////

// Poison values to catch uninitialized variables

#if defined CCTK_REAL_PRECISION_4
constexpr std::uint32_t ipoison = 0xffc00000UL + 0xdead;
#elif defined CCTK_REAL_PRECISION_8
constexpr std::uint64_t ipoison = 0xfff8000000000000ULL + 0xdeadbeef;
#endif
static_assert(sizeof ipoison == sizeof(CCTK_REAL));

// Poison grid functions
void poison_invalid(const GHExt::PatchData::LevelData &leveldata,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = groupdata.valid.at(tl).at(vi).get();
  if (valid.valid_all())
    return;

  static Timer timer("poison_invalid<GF>");
  Interval interval(timer);

  CCTK_REAL poison;
  std::memcpy(&poison, &ipoison, sizeof poison);

  const active_levels_t active_levels(leveldata.level, leveldata.level + 1,
                                      leveldata.patch, leveldata.patch + 1);
  loop_over_blocks(active_levels, [&](int patch, int level, int index,
                                      int block, const cGH *cctkGH) {
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
  synchronize();
}

// Ensure grid functions are not poisoned
void check_valid(const GHExt::PatchData::LevelData &leveldata,
                 const GHExt::PatchData::LevelData::GroupData &groupdata,
                 int vi, int tl, const nan_handling_t nan_handling1,
                 const function<string()> &msg) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

#warning "TODO"
  const nan_handling_t nan_handling = nan_handling_t::forbid_nans;

  const valid_t &valid = groupdata.valid.at(tl).at(vi).get();
  if (!valid.valid_any())
    return;

  static Timer timer("check_valid<GF>");
  Interval interval(timer);

  CCTK_REAL poison;
  std::memcpy(&poison, &ipoison, sizeof poison);

  const auto nan_check = [&](const GridDescBase &grid,
                             const GF3D2<const CCTK_REAL> &gf,
                             const Loop::PointDesc &p) {
    return CCTK_BUILTIN_EXPECT(
        nan_handling == nan_handling_t::allow_nans
            ? std::memcmp(&gf(p.I), &poison, sizeof gf(p.I)) == 0
            : isnan(gf(p.I)),
        false);
  };

  size_t nan_count{0};
  array<int, 3> nan_imin, nan_imax;
  array<CCTK_REAL, 3> nan_xmin, nan_xmax;
  for (int d = 0; d < 3; ++d) {
    nan_imin[d] = numeric_limits<int>::max();
    nan_imax[d] = numeric_limits<int>::min();
    nan_xmin[d] = +1.0 / 0.0;
    nan_xmax[d] = -1.0 / 0.0;
  }
  const auto nan_count_update = [&](const GridDescBase &grid,
                                    const GF3D2<const CCTK_REAL> &gf,
                                    const Loop::PointDesc &p) {
    if (!nan_check(grid, gf, p))
      return;
    static mutex m;
    lock_guard<mutex> g(m);
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
  };

  const active_levels_t active_levels(leveldata.level, leveldata.level + 1,
                                      leveldata.patch, leveldata.patch + 1);
  loop_over_blocks(active_levels, [&](int patch, int level, int index,
                                      int block, const cGH *cctkGH) {
    const Loop::GridDescBaseDevice grid(cctkGH);
    const Loop::GF3D2layout layout(cctkGH, groupdata.indextype);
    const Loop::GF3D2<const CCTK_REAL> gf(
        layout, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, tl, groupdata.firstvarindex + vi)));

    if (valid.valid_all()) {
      grid.loop_idx(
          where_t::everywhere, groupdata.indextype, groupdata.nghostzones,
          [&](const Loop::PointDesc &p) { nan_count_update(grid, gf, p); });
    } else {
      if (valid.valid_int)
        grid.loop_idx(
            where_t::interior, groupdata.indextype, groupdata.nghostzones,
            [&](const Loop::PointDesc &p) { nan_count_update(grid, gf, p); });
      if (valid.valid_outer)
        grid.loop_idx(
            where_t::boundary, groupdata.indextype, groupdata.nghostzones,
            [&](const Loop::PointDesc &p) { nan_count_update(grid, gf, p); });
      if (valid.valid_ghosts)
        grid.loop_idx(
            where_t::ghosts, groupdata.indextype, groupdata.nghostzones,
            [&](const Loop::PointDesc &p) { nan_count_update(grid, gf, p); });
    }
  });
  synchronize();

  if (CCTK_BUILTIN_EXPECT(nan_count > 0, false)) {
#pragma omp critical
    {
      CCTK_VWARN(
          CCTK_WARN_ALERT,
          "%s: Grid function \"%s\" has %td nans, infinities, or poison on "
          "patch %d, "
          "refinement level %d, time level %d, in box [%d,%d,%d]:[%d,%d,%d] "
          "(%g,%g,%g):(%g,%g,%g); expected valid %s",
          msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
          size_t(nan_count), leveldata.patch, leveldata.level, tl, nan_imin[0],
          nan_imin[1], nan_imin[2], nan_imax[0], nan_imax[1], nan_imax[2],
          double(nan_xmin[0]), double(nan_xmin[1]), double(nan_xmin[2]),
          double(nan_xmax[0]), double(nan_xmax[1]), double(nan_xmax[2]),
          groupdata.valid.at(tl).at(vi).explanation().c_str());

      struct info_t {
        where_t where;
        vect<int, dim> I;
        vect<CCTK_REAL, dim> X;
        CCTK_REAL val;
      };
      std::vector<info_t> infos;

      loop_over_blocks(active_levels, [&](int patch, int level, int index,
                                          int block, const cGH *cctkGH) {
        const Loop::GridDescBaseDevice grid(cctkGH);
        const Loop::GF3D2layout layout(cctkGH, groupdata.indextype);
        const Loop::GF3D2<const CCTK_REAL> gf(
            layout, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                        cctkGH, tl, groupdata.firstvarindex + vi)));

        if (valid.valid_int)
          grid.loop_idx(
              where_t::interior, groupdata.indextype, groupdata.nghostzones,
              [&](const Loop::PointDesc &p) {
                if (nan_check(grid, gf, p))
#pragma omp critical(CarpetX_check_valid)
                  infos.push_back(info_t{where_t::interior, p.I, p.X, gf(p.I)});
              });
        if (valid.valid_outer)
          grid.loop_idx(
              where_t::boundary, groupdata.indextype, groupdata.nghostzones,
              [&](const Loop::PointDesc &p) {
                if (nan_check(grid, gf, p))
#pragma omp critical(CarpetX_check_valid)
                  infos.push_back(info_t{where_t::interior, p.I, p.X, gf(p.I)});
              });
        if (valid.valid_ghosts)
          grid.loop_idx(
              where_t::ghosts, groupdata.indextype, groupdata.nghostzones,
              [&](const Loop::PointDesc &p) {
                if (nan_check(grid, gf, p))
#pragma omp critical(CarpetX_check_valid)
                  infos.push_back(info_t{where_t::interior, p.I, p.X, gf(p.I)});
              });
      });
      synchronize();

      std::sort(infos.begin(), infos.end(),
                [](const info_t &a, const info_t &b) {
                  const std::less<vect<int, dim> > lt;
                  return lt(reversed(a.I), reversed(b.I));
                });

      std::ostringstream buf;
      buf << setprecision(std::numeric_limits<CCTK_REAL>::digits10 + 1);
      for (const auto &info : infos)
        buf << "\n"
            << info.where << " " << info.I << " " << info.X << " " << info.val;
      CCTK_WARN(CCTK_WARN_ALERT, buf.str().c_str());

      CCTK_VERROR("%s: Grid function \"%s\" contains nans, infinities, or "
                  "poison on patch %d, "
                  "refinement level %d, time level %d; expected valid %s",
                  msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
                  leveldata.patch, leveldata.level, tl,
                  groupdata.valid.at(tl).at(vi).explanation().c_str());
    }
  }
}

// Poison arrays
void poison_invalid(const GHExt::GlobalData::ArrayGroupData &arraygroupdata,
                    int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = arraygroupdata.valid.at(tl).at(vi).get();
  if (valid.valid_all())
    return;

  static Timer timer("poison_invalid<GA>");
  Interval interval(timer);

  CCTK_REAL poison;
  std::memcpy(&poison, &ipoison, sizeof poison);

  if (!valid.valid_int) {
    int dimension = arraygroupdata.dimension;
    CCTK_REAL *restrict const ptr =
        const_cast<CCTK_REAL *>(&arraygroupdata.data.at(tl).at(vi));
    const int *gsh = arraygroupdata.gsh;
    int n_elems = 1;
    for (int i = 0; i < dimension; i++)
      n_elems *= gsh[i];
    for (int i = 0; i < n_elems; i++)
      ptr[i] = poison;
  }
}

// Ensure arrays are not poisoned
void check_valid(const GHExt::GlobalData::ArrayGroupData &arraygroupdata,
                 int vi, int tl, const nan_handling_t nan_handling,
                 const function<string()> &msg) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = arraygroupdata.valid.at(tl).at(vi).get();
  if (!valid.valid_any())
    return;

  static Timer timer("check_valid<GA>");
  Interval interval(timer);

  CCTK_REAL poison;
  std::memcpy(&poison, &ipoison, sizeof poison);

  // arrays have no boundary so we expect them to alway be valid
  assert(valid.valid_outer && valid.valid_ghosts);

  size_t nan_count{0};
  if (valid.valid_int) {
    const CCTK_REAL *restrict const ptr = &arraygroupdata.data.at(tl).at(vi);
    int dimension = arraygroupdata.dimension;
    const int *gsh = arraygroupdata.gsh;
    int n_elems = 1;
    for (int i = 0; i < dimension; i++)
      n_elems *= gsh[i];
    for (int i = 0; i < n_elems; i++) {
      if (CCTK_BUILTIN_EXPECT(nan_handling == nan_handling_t::allow_nans
                                  ? std::memcmp(ptr, &poison, sizeof *ptr) == 0
                                  : isnan(*ptr),
                              false))
        ++nan_count;
    }
  }

  if (CCTK_BUILTIN_EXPECT(nan_count > 0, false))
    CCTK_VERROR("%s: Array \"%s\" has %td nans on time level %d; "
                "expected valid %s",
                msg().c_str(),
                CCTK_FullVarName(arraygroupdata.firstvarindex + vi), nan_count,
                tl, arraygroupdata.valid.at(tl).at(vi).explanation().c_str());
}

////////////////////////////////////////////////////////////////////////////////

// Checksums to catch illegal modifications

checksums_t
calculate_checksums(const vector<vector<vector<valid_t> > > &will_write) {
  DECLARE_CCTK_PARAMETERS;

  checksums_t checksums;

  if (!poison_undefined_values)
    return checksums;

  static Timer timer("calculate_checksums");
  Interval interval(timer);

  assert(active_levels);
  active_levels->loop([&](auto &restrict leveldata) {
    auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {

      for (const auto &groupdataptr : leveldata.groupdata) {
        if (groupdataptr == nullptr)
          continue;

        auto &restrict groupdata = *groupdataptr;
        const GridPtrDesc1 grid(leveldata, groupdata, mfi);

        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          for (int tl = 0; tl < int(groupdata.valid.size()); ++tl) {
            const tiletag_t tiletag{leveldata.patch,
                                    leveldata.level,
                                    mfi.tilebox(),
                                    groupdata.groupindex,
                                    vi,
                                    tl};

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

            const amrex::Array4<const CCTK_REAL> &vars =
                groupdata.mfab.at(tl)->array(mfi);
            const GF3D1<const CCTK_REAL> var_ = grid.gf3d(vars, vi);

            checksum_t checksum(to_check);
            checksum.add(tiletag);
            const auto add_point{
                [&](const Loop::PointDesc &p) { checksum.add(var_(p.I)); }};

            if (to_check.valid_int)
              grid.loop_idx(where_t::interior, groupdata.indextype,
                            groupdata.nghostzones, add_point);

            if (to_check.valid_outer)
              grid.loop_idx(where_t::boundary, groupdata.indextype,
                            groupdata.nghostzones, add_point);

            if (to_check.valid_ghosts)
              grid.loop_idx(where_t::ghosts, groupdata.indextype,
                            groupdata.nghostzones, add_point);

#pragma omp critical
            checksums[tiletag] = checksum;
          }
        }
      }
    }
  });

  return checksums;
}

void check_checksums(const checksums_t &checksums,
                     const std::function<string()> &where) {
  DECLARE_CCTK_PARAMETERS;

  if (!poison_undefined_values)
    return;
  if (checksums.empty())
    return;

  static Timer timer("check_checksums");
  Interval interval(timer);

  assert(active_levels);
  active_levels->loop([&](auto &restrict leveldata) {
    auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {

      for (const auto &groupdataptr : leveldata.groupdata) {
        if (groupdataptr == nullptr)
          continue;

        auto &restrict groupdata = *groupdataptr;
        const GridPtrDesc1 grid(leveldata, groupdata, mfi);

        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          for (int tl = 0; tl < int(groupdata.valid.size()); ++tl) {
            const tiletag_t tiletag{leveldata.patch,
                                    leveldata.level,
                                    mfi.tilebox(),
                                    groupdata.groupindex,
                                    vi,
                                    tl};

            if (!checksums.count(tiletag))
              continue;

            const auto &old_checksum = checksums.at(tiletag);
            const auto &did_check = old_checksum.where;
            assert(did_check.valid_any());

            const amrex::Array4<const CCTK_REAL> &vars =
                groupdata.mfab.at(tl)->array(mfi);
            const GF3D1<const CCTK_REAL> var_ = grid.gf3d(vars, vi);

            checksum_t checksum(did_check);
            checksum.add(tiletag);
            const auto add_point{
                [&](const Loop::PointDesc &p) { checksum.add(var_(p.I)); }};

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
              CCTK_VERROR(
                  "%s: Checksum mismatch: variable %s, tile %s, "
                  "int:%d,outer:%d,ghosts:%d, old checksum %s, new checksum %s",
                  where().c_str(),
                  CCTK_FullVarName(groupdata.firstvarindex + tiletag.vi),
                  string(tiletag).c_str(), int(did_check.valid_int),
                  int(did_check.valid_outer), int(did_check.valid_ghosts),
                  string(old_checksum).c_str(), string(checksum).c_str());
          }
        }
      }
    }
  });
}

} // namespace CarpetX
