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

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace Arith;

using Loop::dim;

using rat64 = rational<int64_t>;

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

      friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                       const LevelData &leveldata);
    };
    std::vector<LevelData> leveldata; // [reflevel]

    friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                     const PatchData &patchdata);
  };
  std::vector<PatchData> patchdata; // [patch]

  int num_patches() const { return patchdata.size(); }
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

amrex::Interpolater *get_interpolator(const std::array<int, dim> indextype);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_DRIVER_HXX
