#ifndef CARPETX_CARPETX_DRIVER_HXX
#define CARPETX_CARPETX_DRIVER_HXX

#include "loop.hxx"
#include "valid.hxx"

#include <rational.hxx>
#include <tuple.hxx>

#include <cctk.h>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace Arith;

using Loop::dim;

using rat64 = rational<int64_t>;

// Compile-time capacity of the subcycling band arrays. The runtime-active
// count lives in GHExt::num_rk_stages (RK4 uses 4, SSPRK3 uses 3).
inline constexpr int max_num_rk_stages = 4;

// TODO: It seems that AMReX now also has `RB90`, `RB180`, and
// `PolarB` boundary conditions. Make these available as well.

// Symmetries are domain properties
enum class symmetry_t {
  none,
  interpatch,
  periodic,
  reflection,
};
std::ostream &operator<<(std::ostream &os, const symmetry_t symmetry);

// Boundary conditions are group properties. They are valid only for faces where
// the domain symmetry is `none`.
enum class boundary_t {
  none,
  symmetry_boundary,
  dirichlet,
  linear_extrapolation,
  neumann,
  robin,
};
std::ostream &operator<<(std::ostream &os, const boundary_t boundary);

static_assert(AMREX_SPACEDIM == dim,
              "AMReX's AMREX_SPACEDIM must be the same as Cactus's cctk_dim");

static_assert(is_same<amrex::Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

////////////////////////////////////////////////////////////////////////////////

// AMR driver
class CactusAmrCore final : public amrex::AmrCore {
  int patch;

public:
  bool cactus_is_initialized = false;
  std::vector<bool> level_modified;

  CactusAmrCore();
  CactusAmrCore(int patch, const amrex::RealBox *rb, int max_level_in,
                const amrex::Vector<int> &n_cell_in, int coord = -1,
                amrex::Vector<amrex::IntVect> ref_ratios =
                    amrex::Vector<amrex::IntVect>(),
                const int *is_per = nullptr);
  CactusAmrCore(int patch, const amrex::RealBox &rb, int max_level_in,
                const amrex::Vector<int> &n_cell_in, int coord,
                amrex::Vector<amrex::IntVect> const &ref_ratios,
                amrex::Array<int, AMREX_SPACEDIM> const &is_per);
  CactusAmrCore(const amrex::AmrCore &rhs) = delete;
  CactusAmrCore &operator=(const amrex::AmrCore &rhs) = delete;

  virtual ~CactusAmrCore() override;

  virtual void ErrorEst(int level, amrex::TagBoxArray &tags, amrex::Real time,
                        int ngrow) override;
  void SetupLevel(int level, const amrex::BoxArray &ba,
                  const amrex::DistributionMapping &dm,
                  const function<string()> &why);
  // Re-partition a recovered level's covered region into a fresh box
  // decomposition matching the current max_grid_size / node count. Pure
  // geometry: the returned BoxArray covers exactly the same cell union as the
  // input, so no field data is changed. Level 0 regenerates the canonical
  // full-domain decomposition; levels > 0 coalesce then re-tile + load-balance.
  amrex::BoxArray RechopLevel(int level, amrex::BoxArray ba) const;
  virtual void
  MakeNewLevelFromScratch(int level, amrex::Real time,
                          const amrex::BoxArray &ba,
                          const amrex::DistributionMapping &dm) override;
  virtual void
  MakeNewLevelFromCoarse(int level, amrex::Real time, const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm) override;
  virtual void RemakeLevel(int level, amrex::Real time,
                           const amrex::BoxArray &ba,
                           const amrex::DistributionMapping &dm) override;
  virtual void ClearLevel(int level) override;
};

// Cactus grid hierarchy extension
struct GHExt {

  GHExt() = default;
  GHExt(const GHExt &) = delete;
  GHExt(GHExt &&) = delete;
  GHExt &operator=(const GHExt &) = delete;
  GHExt &operator=(GHExt &&) = delete;

  struct cctkGHptr {
    cGH *cctkGH;
    cctkGHptr(const cctkGHptr &) = delete;
    cctkGHptr(cctkGHptr &&ptr) : cctkGH(ptr.cctkGH) { ptr.cctkGH = nullptr; }
    cctkGHptr &operator=(const cctkGHptr &) = delete;
    cctkGHptr &operator=(cctkGHptr &&ptr);
    cctkGHptr() : cctkGH(nullptr) {}
    cctkGHptr(cGH *&&cctkGH) : cctkGH(cctkGH) {}
    cctkGHptr &operator=(cGH *&&cctkGH);
    ~cctkGHptr();
    operator bool() const { return bool(cctkGH); }
    cGH *get() const { return cctkGH; }
  };

  cctkGHptr global_cctkGH;
  std::vector<cctkGHptr> level_cctkGHs; // [reflevel]

  struct CommonGroupData {
    std::string groupname;
    int groupindex;
    int firstvarindex;
    int numvars;

    bool do_checkpoint; // whether to checkpoint
    bool do_evolve;     // whether this is an evolved state variable
    bool do_restrict;   // whether to restrict

    std::vector<std::vector<why_valid_t> > valid; // [time level][var index]

    // TODO: add poison_invalid and check_valid functions

    friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                     const CommonGroupData &commongroupdata);
  };

  struct GlobalData {
    // all data that exists on all levels

    class AnyTypeVector {

    public:
      // access to a single element of a AnyTypeVector
      class AnyTypeScalarRef {
      public:
        AnyTypeScalarRef() = delete;
        AnyTypeScalarRef(const AnyTypeVector &vect_, size_t idx_)
            : _vect(vect_), _idx(idx_) {}

      private:
        const AnyTypeVector &_vect;
        const size_t _idx;

        friend YAML::Emitter &
        operator<<(YAML::Emitter &yaml,
                   const AnyTypeScalarRef &anytypescalarref);
        friend std::ostream &operator<<(std::ostream &os,
                                        const AnyTypeScalarRef &scalar);
      };

      AnyTypeVector() : _type(-1), _typesize(-1), _count(0), _data(nullptr) {};
      AnyTypeVector(int type_, size_t count_)
          : _type(-1), _typesize(-1), _count(0), _data(nullptr) {
        alloc(type_, count_);
        assert(_type == type_);
        assert(_typesize != -1);
        assert(_count == count_);
        assert(_data != nullptr);
      };
      // Noncopyable for now
      AnyTypeVector(const AnyTypeVector &) = delete;
      AnyTypeVector &operator=(const AnyTypeVector &) = delete;
      AnyTypeVector &operator=(AnyTypeVector &&other) {
        swap(other);
        return *this;
      }
      AnyTypeVector(AnyTypeVector &&other)
          : _type(other._type), _typesize(other._typesize),
            _count(other._count), _data(other._data) {
        other._type = -1;
        other._typesize = -1;
        other._count = 0;
        other._data = nullptr;
      }
      void swap(AnyTypeVector &other) {
        std::swap(this->_type, other._type);
        std::swap(this->_typesize, other._typesize);
        std::swap(this->_count, other._count);
        std::swap(this->_data, other._data);
      }

      ~AnyTypeVector() {
        if (_data != nullptr) {
          assert(_type != -1);
          assert(_typesize != -1);
          amrex::The_Arena()->free(_data);
          _type = -1;
          _typesize = -1;
          _count = 0;
          _data = nullptr;
        }
        assert(_type == -1);
        assert(_typesize == -1);
        assert(_count == 0);
        assert(_data == nullptr);
      };

      void alloc(int type_, size_t count_) {
        assert(type_ == CCTK_VARIABLE_INT || type_ == CCTK_VARIABLE_REAL ||
               type_ == CCTK_VARIABLE_COMPLEX);

        assert(_type == -1);
        assert(_typesize == -1);
        assert(_count == 0);
        assert(_data == nullptr);

        _type = type_;
        _typesize = CCTK_VarTypeSize(_type);
        assert(_typesize > 0);
        _count = count_;
        _data = amrex::The_Arena()->alloc(_typesize * _count);
      }

      void free() {
        assert(_type != -1);
        assert(_typesize != -1);
        assert(_data != nullptr);
        amrex::The_Arena()->free(_data);
        _type = -1;
        _typesize = -1;
        _count = 0;
        _data = nullptr;
      }

      int type() const { return _type; };
      int typesize() const { return _typesize; };

      const void *data_at(size_t i) const {
#ifdef CCTK_DEBUG
        if (i >= _count) {
          CCTK_VERROR("invalid index %zd exceeds %zd", i, _count);
        }
#endif
        assert(i < _count);
        return (char *)_data + i * _typesize;
      };

      void *data_at(size_t i) {
#ifdef CCTK_DEBUG
        if (i >= _count) {
          CCTK_VERROR("invalid index %zu exceeds %zu", i, _count);
        }
#endif
        assert(i < _count);
        return (char *)_data + i * _typesize;
      };

      AnyTypeScalarRef operator[](size_t idx) const {
        return AnyTypeScalarRef(*this, idx);
      }

      size_t size() const { return _count; };

      friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                       const AnyTypeVector &anytypevector);

    private:
      int _type, _typesize;
      size_t _count;
      void *_data;
    };

    // For subcycling in time, there really should be one copy of each
    // integrated grid scalar per level. We don't do that yet; instead,
    // we assume that grid scalars only hold "analysis" data.

    struct ArrayGroupData : public CommonGroupData {
      vector<AnyTypeVector> data; // [time level][var index + grid point index]
      int array_size;
      int dimension;
      int activetimelevels;
      int lsh[dim];
      int ash[dim];
      int gsh[dim];
      int lbnd[dim];
      int ubnd[dim];
      int bbox[2 * dim];
      int nghostzones[dim];

      ArrayGroupData() {
        array_size = -1;
        dimension = -1;
        activetimelevels = -1;
        for (int d = 0; d < dim; d++) {
          lsh[d] = -1;
          ash[d] = -1;
          gsh[d] = -1;
          lbnd[d] = -1;
          ubnd[d] = -1;
          bbox[2 * d] = bbox[2 * d + 1] = -1;
          nghostzones[d] = -1;
        }
      }

      friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                       const ArrayGroupData &arraygroupdata);
    };
    // TODO: right now this is sized for the total number of groups
    std::vector<std::unique_ptr<ArrayGroupData> >
        arraygroupdata; // [group index]

    friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                     const GlobalData &globaldata);
  };
  GlobalData globaldata;

  struct PatchData {
    PatchData() = delete;
    PatchData(const PatchData &) = delete;
    PatchData &operator=(const PatchData &) = delete;
    PatchData(PatchData &&) = default;
    PatchData &operator=(PatchData &&) = default;

    PatchData(int patch);

    int patch;

    bool is_cartesian;

    std::array<std::array<symmetry_t, dim>, 2> symmetries;

    // AMReX grid structure
    // TODO: convert this from unique_ptr to optional
    std::unique_ptr<CactusAmrCore> amrcore;

    struct LevelData {
      LevelData() = delete;
      LevelData(const LevelData &) = delete;
      LevelData &operator=(const LevelData &) = delete;
      LevelData(LevelData &&) = default;
      LevelData &operator=(LevelData &&) = default;

      LevelData(const int patch, const int level, const amrex::BoxArray &ba,
                const amrex::DistributionMapping &dm,
                const std::function<std::string()> &why);

      int patch, level;

      // This level uses subcycling with respect to the next coarser
      // level. (Ignored for the coarsest level.)
      bool is_subcycling_level;

      // Iteration and time at which this cycle level is valid
      rat64 iteration, delta_iteration;

      // Fabamrex::ArrayBase object holding a cell-centred BoxArray for
      // iterating over grid functions. This stores the grid structure
      // and its distribution over all processes, but holds no data.
      std::unique_ptr<amrex::FabArrayBase> fab;

      // Per-centering coarse-fine ghost masks, built single-threaded by
      // build_cf_mask and read by get_cf_mask.
      // Indexed by (indextype[0]<<2)|(indextype[1]<<1)|indextype[2].
      mutable std::array<std::unique_ptr<amrex::iMultiFab>, 8> cf_masks;

      // Per-centering coarse-fine boundary-band geometry for subcycling RK
      // k-stages, built lazily by build_bands and used to allocate the
      // per-group band MultiFabs. Indexed by centering s, mirroring cf_masks.
      //   source_band_*   : coarse cells under this level's children's cf-ghost
      //                     footprint (child fpc.ba_crse_patch). Empty on the
      //                     finest level.
      //   consumer_band_* : this level's own cf-ghost region (fpc.ba_fine_patch
      //                     w.r.t. the parent). Empty at level 0.
      // A non-null BoxArray slot (even if the BoxArray itself is empty) means
      // the geometry for that centering has been built; both slots are set
      // together. The DistributionMapping is the fpc.dm_patch the consumer band
      // must share for the band->band FillPatchInterp to stay local. Shared by
      // both band families (ks_* and old_*).
      mutable std::array<std::unique_ptr<amrex::BoxArray>, 8> source_band_ba;
      mutable std::array<std::unique_ptr<amrex::DistributionMapping>, 8>
          source_band_dm;
      mutable std::array<std::unique_ptr<amrex::BoxArray>, 8> consumer_band_ba;
      mutable std::array<std::unique_ptr<amrex::DistributionMapping>, 8>
          consumer_band_dm;

      // The child (level+1) BoxArray the source band was last built against.
      // A mismatch with the current child layout (which AMReX may have changed
      // without re-making this coarser level) rebuilds the source band. null
      // means not yet built; an empty BoxArray means there was no child.
      mutable std::array<std::unique_ptr<amrex::BoxArray>, 8>
          source_band_child_ba;

      // Returns the coarse-fine ghost mask for this (level, centering), or
      // nullptr at level 0 / when subcycling is disabled. Pure reader,
      // side-effect-free and safe to call from a parallel consume; callers MUST
      // warm the centerings single-threaded via build_cf_mask first.
      amrex::iMultiFab *
      get_cf_mask(const std::array<int, dim> &indextype,
                  const std::array<int, dim> &nghostzones) const;

      // Build and cache the coarse-fine ghost mask for this (level,
      // centering). Idempotent; a no-op at level 0 / when subcycling is
      // disabled. The only caller of iMultiFab::BuildMask, which opens its own
      // OpenMP region, so this must run single-threaded (no active MFIter, no
      // enclosing parallel region).
      void build_cf_mask(const std::array<int, dim> &indextype,
                         const std::array<int, dim> &nghostzones) const;

      cctkGHptr patch_cctkGH;
      std::vector<cctkGHptr> local_cctkGHs; // [component]

      cGH *get_patch_cctkGH() const { return patch_cctkGH.get(); }
      cGH *get_local_cctkGH(const int component) const {
        return local_cctkGHs.at(component).get();
      }

      struct GroupData : public CommonGroupData {
        GroupData() = delete;
        GroupData(const GroupData &) = delete;
        GroupData &operator=(const GroupData &) = delete;
        GroupData(GroupData &&) = delete;
        GroupData &operator=(GroupData &&) = delete;

        GroupData(int patch, int level, int gi, const amrex::BoxArray &ba,
                  const amrex::DistributionMapping &dm,
                  const std::function<std::string()> &why);

        int patch, level;

        std::array<int, dim> indextype;
        std::array<int, dim> nghostzones;

        amrex::Interpolater *interpolator;

        std::array<std::array<boundary_t, dim>, 2> boundaries;
        bool all_faces_have_symmetries_or_boundaries() const;
        std::vector<array<int, dim> > parities;
        std::vector<CCTK_REAL> dirichlet_values;
        std::vector<CCTK_REAL> robin_values;
        amrex::Vector<amrex::BCRec> bcrecs;

        // Apply outer (physical) boundary conditions to a MultiFab
        void apply_boundary_conditions(amrex::MultiFab &mfab) const;

        // each amrex::MultiFab has numvars components
        std::vector<std::unique_ptr<amrex::MultiFab> > mfab; // [time level]

        // Coarse-fine boundary bands holding the subcycling RK k-stages
        // (zero-ghost, numvars comps), allocated lazily by build_bands only
        // under subcycling for evolved groups. Indexed by RK stage.
        //   ks_source_band[s]  : coarse cells under children's cf-ghost
        //                        footprint, written by setks and read as the
        //                        prolongation source. Empty on the finest
        //                        level.
        //   ks_consumer_band[s]: this level's own cf-ghost region, filled by
        //                        prolongation from the parent and read by the
        //                        dense-output kernel. Empty at level 0.
        mutable std::array<std::unique_ptr<amrex::MultiFab>, max_num_rk_stages>
            ks_source_band;
        mutable std::array<std::unique_ptr<amrex::MultiFab>, max_num_rk_stages>
            ks_consumer_band;

        // Coarse-fine boundary bands holding the subcycling old state u(t_n), a
        // single snapshot (not RK-stage indexed) sharing the ks bands' geometry
        // and lifecycle above. old_source_band is filled from var(tl=0) at
        // solve start; old_consumer_band is its prolongation, read by the
        // dense-output kernel as the u(t_n) base.
        mutable std::unique_ptr<amrex::MultiFab> old_source_band;
        mutable std::unique_ptr<amrex::MultiFab> old_consumer_band;

        // flux register between this and the next coarser level
        std::unique_ptr<amrex::FluxRegister> freg;
        // associated flux group indices
        std::array<int, dim> fluxes; // [dir]

        // CarpetX can allocate and free (temporary) multifabs that
        // are associated with a Cactus grid function group. These
        // multifabs remain allocated when they are freed, which makes
        // it efficient when they are re-allocated later. However,
        // they are freed when the current level changes during
        // regridding (and the shape of the multifab presumably
        // changes). This is used e.g. by ODESolvers for its
        // temporaries.
      private:
        mutable std::vector<std::unique_ptr<amrex::MultiFab> > tmp_mfabs;
        mutable std::size_t next_tmp_mfab;

      public:
        void init_tmp_mfabs() const;
        amrex::MultiFab *alloc_tmp_mfab() const;
        void free_tmp_mfabs() const;

        friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                         const GroupData &groupdata);
      };
      // TODO: right now this is sized for the total number of groups
      std::vector<unique_ptr<GroupData> > groupdata; // [group index]

      // Build (lazily, idempotently) the coarse-fine band geometry for this
      // group's centering and allocate the group's ks_source_band/
      // ks_consumer_band MultiFabs (zero ghost, numvars comps). A no-op when
      // subcycling is disabled or the group is not evolved. Computes the
      // source-band geometry from the next-finer level's fpc, so it must run
      // after all levels exist; like build_cf_mask it warms a cache and must
      // run single-threaded.
      void build_bands(const GroupData &groupdata) const;

      friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                       const LevelData &leveldata);
    };
    std::vector<LevelData> leveldata; // [reflevel]

    friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                     const PatchData &patchdata);
  };
  std::vector<PatchData> patchdata; // [patch]

  // Number of active timelevels per group, populated from schedule.ccl STORAGE
  // statements during CCTKi_ScheduleGHInit. Indexed by group index; 0 means no
  // storage. Frozen between schedule init and the first allocation pass.
  std::vector<int> active_timelevels; // [group index]
  bool storage_frozen = false;

  bool use_subcycling = false;

  // Active number of RK stages for subcycling, set from ODESolvers::method at
  // WRAGH (SSPRK3 -> 3, else 4). Must be <= max_num_rk_stages.
  int num_rk_stages = 4;

  // Per-level iteration values read from checkpoint; consumed by recovery fixup
  // in schedule.cxx. Indexed [patch][level]. Empty outside of recovery window.
  std::vector<std::vector<std::optional<rat64> > > recovered_level_iterations;

  int num_patches() const { return patchdata.size(); }
  int num_levels(const int patch) const {
    return patchdata.at(patch).leveldata.size();
  }
  int num_levels() const {
    int nlevels = 0;
    using std::max;
    for (const auto &pd : patchdata)
      nlevels = max(nlevels, int(pd.leveldata.size()));
    return nlevels;
  }

  cGH *get_global_cctkGH() const { return global_cctkGH.get(); }
  cGH *get_level_cctkGH(const int level) const {
    return level_cctkGHs.at(level).get();
  }
  cGH *get_patch_cctkGH(const int level, const int patch) const {
    return patchdata.at(patch).leveldata.at(level).patch_cctkGH.get();
  }
  cGH *get_local_cctkGH(const int level, const int patch,
                        const int component) const {
    return patchdata.at(patch)
        .leveldata.at(level)
        .local_cctkGHs.at(component)
        .get();
  }

  friend YAML::Emitter &operator<<(YAML::Emitter &yaml, const GHExt &ghext);
  friend std::ostream &operator<<(std::ostream &os, const GHExt &ghext);
};

extern std::unique_ptr<GHExt> ghext;

// True iff every level of every patch sits at the same subcycling iteration,
// i.e. the checkpoint is time-aligned. Always true without subcycling. When
// false, the fine consumer bands hold mid-cycle state that must be serialized.
bool all_levels_synchronized();

// Subcycling consumer-band kinds serialized at unsynchronized checkpoints:
// ks_consumer is the RK stages 0..max_num_rk_stages-1, old_consumer the u(t_n)
// snapshot.
enum class band_kind { ks_consumer, old_consumer };

// Token identifying a consumer band in checkpoint names: "ksc_s00".."ksc_s03"
// for ks_consumer, "oldc" for old_consumer. Shared by both IO backends.
std::string subcycling_band_tag(band_kind kind, int stage = -1);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_DRIVER_HXX
