#ifndef CARPETX_CARPETX_SCHEDULE_HXX
#define CARPETX_CARPETX_SCHEDULE_HXX

#include "driver.hxx"
#include "loop.hxx"

#include <cctk.h>
#include <cctk_Schedule.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <optional>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace std;
using namespace Loop;

int Initialise(tFleshConfig *config);
int Evolve(tFleshConfig *config);
int Shutdown(tFleshConfig *config);

int SyncGroupsByDirI(const cGH *restrict cctkGH, int numgroups,
                     const int *groups, const int *directions);

int CallFunction(void *function, cFunctionData *attribute, void *data);

int GroupStorageIncrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status);
int GroupStorageDecrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status);
int EnableGroupStorage(const cGH *cctkGH, const char *groupname);
int DisableGroupStorage(const cGH *cctkGH, const char *groupname);

////////////////////////////////////////////////////////////////////////////////

struct active_levels_t {
  int min_level, max_level;
  int min_patch, max_patch;

  // active_levels_t() = delete;

  active_levels_t(const int min_level, const int max_level, const int min_patch,
                  const int max_patch)
      : min_level(min_level), max_level(max_level), min_patch(min_patch),
        max_patch(max_patch) {
    assert(min_level >= 0);
    assert(max_level <= ghext->num_levels());
    assert(min_patch >= 0);
    assert(max_patch <= ghext->num_patches());
  }
  active_levels_t(const int min_level, const int max_level)
      : active_levels_t(min_level, max_level, 0, ghext->num_patches()) {}
  active_levels_t() : active_levels_t(0, ghext->num_levels()) {}

private:
  void assert_consistent_iterations() const {
    rat64 good_iteration = -1;
    for (int level = min_level; level < max_level; ++level) {
      for (int patch = min_patch; patch < max_patch; ++patch) {
        const auto &patchdata = ghext->patchdata.at(patch);
        if (level < int(patchdata.leveldata.size())) {
          const auto &leveldata = patchdata.leveldata.at(level);
          const auto &iteration = leveldata.iteration;
          if (good_iteration == -1)
            good_iteration = iteration;
          assert(iteration == good_iteration);
        }
      }
    }
  }

public:
  // Loop over all active patches of all active levels
  template <typename F> void loop(F f) const {
    assert_consistent_iterations();
    for (int level = min_level; level < max_level; ++level) {
      for (int patch = min_patch; patch < max_patch; ++patch) {
        auto &patchdata = ghext->patchdata.at(patch);
        if (level < int(patchdata.leveldata.size()))
          f(patchdata.leveldata.at(level));
      }
    }
  }

  // Loop over all active patches of all active levels
  template <typename F> void loop_reverse(F f) const {
    assert_consistent_iterations();
    for (int level = max_level - 1; level >= min_level; --level) {
      for (int patch = min_patch; patch < max_patch; ++patch) {
        auto &patchdata = ghext->patchdata.at(patch);
        if (level < int(patchdata.leveldata.size()))
          f(patchdata.leveldata.at(level));
      }
    }
  }
};

// The levels CallFunction should traverse
// TODO: Move this into ghext
extern optional<active_levels_t> active_levels;

////////////////////////////////////////////////////////////////////////////////

// Like an MFIter, but does not support iteration, instead it can be copied
struct MFPointer {
  int m_index;
  amrex::Box m_validbox; // interior of component
  amrex::Box m_tilebox;  // interior of tile

  MFPointer() = delete;
  MFPointer(const MFPointer &) = default;
  MFPointer(MFPointer &&) = default;
  MFPointer &operator=(const MFPointer &) = default;
  MFPointer &operator=(MFPointer &&) = default;
  MFPointer(const amrex::MFIter &mfi)
      : m_index((assert(mfi.isValid()), mfi.index())),
        m_validbox(mfi.validbox()), m_tilebox(mfi.tilebox()) {}

  constexpr int index() const noexcept { return m_index; }

  constexpr amrex::Box validbox() const noexcept { return m_validbox; }

  constexpr amrex::Box tilebox() const noexcept { return m_tilebox; }

  const amrex::Box fabbox(const amrex::IntVect &ng) const noexcept {
    return amrex::grow(validbox(), ng);
  }

  const amrex::Box growntilebox(const amrex::IntVect &ng) const noexcept {
    // return m_growntilebox;
    amrex::Box bx = tilebox();
    const amrex::Box &vbx = validbox();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      if (bx.smallEnd(d) == vbx.smallEnd(d))
        bx.growLo(d, ng[d]);
      if (bx.bigEnd(d) == vbx.bigEnd(d))
        bx.growHi(d, ng[d]);
    }
    return bx;
  }
};

struct GridDesc : GridDescBase {

  GridDesc() = delete;
  GridDesc(const GHExt::PatchData::LevelData &leveldata, const MFPointer &mfp);
  GridDesc(const GHExt::PatchData::LevelData &leveldata,
           const int global_component);
  GridDesc(const cGH *cctkGH) : GridDescBase(cctkGH) {}
};

// TODO: remove this
struct GridPtrDesc : GridDesc {
  amrex::Dim3 cactus_offset;

  GridPtrDesc() = delete;
  GridPtrDesc(const GHExt::PatchData::LevelData &leveldata,
              const MFPointer &mfp);

  template <typename T> T *ptr(const amrex::Array4<T> &vars, int vi) const {
    return vars.ptr(cactus_offset.x, cactus_offset.y, cactus_offset.z, vi);
  }
  template <typename T>
  T &idx(const amrex::Array4<T> &vars, int i, int j, int k, int vi) const {
    return vars(cactus_offset.x + i, cactus_offset.y + i, cactus_offset.z + j,
                vi);
  }
};

struct GridPtrDesc1 : GridDesc {
  amrex::Dim3 cactus_offset;
  array<int, dim> gimin, gimax;
  array<int, dim> gash;

  GridPtrDesc1() = delete;
  GridPtrDesc1(const GridPtrDesc1 &) = delete;
  GridPtrDesc1 &operator=(const GridPtrDesc1 &) = delete;

  GridPtrDesc1(const GHExt::PatchData::LevelData &leveldata,
               const GHExt::PatchData::LevelData::GroupData &groupdata,
               const MFPointer &mfp);

  template <typename T> T *ptr(const amrex::Array4<T> &vars, int vi) const {
    return vars.ptr(cactus_offset.x + gimin[0], cactus_offset.y + gimin[1],
                    cactus_offset.z + gimin[2], vi);
  }
  template <typename T>
  T &idx(const amrex::Array4<T> &vars, int i, int j, int k, int vi) const {
    return vars(cactus_offset.x + gimin[0] + i, cactus_offset.y + gimin[1] + j,
                cactus_offset.z + gimin[2] + k, vi);
  }

  template <typename T>
  GF3D1<T> gf3d(const amrex::Array4<T> &vars, int vi) const {
    return GF3D1<T>(ptr(vars, vi), gimin, gimax, gash);
  }

  friend ostream &operator<<(ostream &os, const GridPtrDesc1 &p) {
    os << "GridPtrDesc1{" << (const GridDescBase &)p << ", "
       << "cactus_offset:"
       << "{" << p.cactus_offset.x << "," << p.cactus_offset.y << ","
       << p.cactus_offset.z << "}, "
       << "gimin:"
       << "{" << p.gimin[0] << "," << p.gimin[1] << "," << p.gimin[2] << "}, "
       << "gimax:"
       << "{" << p.gimax[0] << "," << p.gimax[1] << "," << p.gimax[2] << "}, "
       << "gash:"
       << "{" << p.gash[0] << "," << p.gash[1] << "," << p.gash[2] << "}";
    return os;
  }
};

////////////////////////////////////////////////////////////////////////////////

cGH *copy_cctkGH(const cGH *restrict const sourceGH);
void delete_cctkGH(cGH *cctkGH);

bool in_local_mode(const cGH *restrict cctkGH);
bool in_patch_mode(const cGH *restrict cctkGH);
bool in_level_mode(const cGH *restrict cctkGH);
bool in_global_mode(const cGH *restrict cctkGH);
bool in_meta_mode(const cGH *restrict cctkGH);

void update_cctkGH(cGH *restrict cctkGH, const cGH *restrict sourceGH);
void enter_global_mode(cGH *restrict cctkGH);
void leave_global_mode(cGH *restrict cctkGH);
void enter_level_mode(cGH *restrict cctkGH, int level);
void leave_level_mode(cGH *restrict cctkGH, int level);
void enter_patch_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata);
void leave_patch_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata);
void enter_local_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata,
                      const MFPointer &mfp);
void leave_local_mode(cGH *restrict cctkGH,
                      const GHExt::PatchData::LevelData &restrict leveldata,
                      const MFPointer &mfp);

// Loop over all components of a single patch and level
void loop_over_components(
    amrex::FabArrayBase &fab,
    const std::function<void(int index, int component)> &component_kernel);
// Loop over all components of several patches and levels
void loop_over_components(
    const active_levels_t &active_levels,
    const std::function<void(int patch, int level, int index, int component,
                             const cGH *cctkGH)> &component_kernel);
void synchronize();

// These functions are defined in valid.cxx. These prototypes should
// be moved to valid.hxx. Unfortunately, they depend on GHExt, which is declared
// in driver.hxx, which includes valid.hxx. Declaring the prorotypes here avoids
// that circular reference. This should be fixed.

void error_if_invalid(const GHExt::PatchData::LevelData::GroupData &grouppdata,
                      int vi, int tl, const valid_t &required,
                      const function<string()> &msg);
void warn_if_invalid(const GHExt::PatchData::LevelData::GroupData &grouppdata,
                     int vi, int tl, const valid_t &required,
                     const function<string()> &msg);

void error_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata,
                      int vi, int tl, const valid_t &required,
                      const function<string()> &msg);
void warn_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                     int tl, const valid_t &required,
                     const function<string()> &msg);

enum class nan_handling_t { allow_nans, forbid_nans };

void poison_invalid(const GHExt::PatchData::LevelData &leveldata,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    int vi, int tl);
void check_valid(const GHExt::PatchData::LevelData &leveldata,
                 const GHExt::PatchData::LevelData::GroupData &groupdata,
                 int vi, int tl, nan_handling_t nan_handling,
                 const function<string()> &msg);

void poison_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                    int tl);
void check_valid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                 int tl, nan_handling_t nan_handling,
                 const function<string()> &msg);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_SCHEDULE_HXX
