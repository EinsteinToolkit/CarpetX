#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace amrex {
struct IntVect {
  std::array<int, 3> v{};
  int operator[](int d) const { return v.at(static_cast<std::size_t>(d)); }
};
struct IndexType { IntVect value; IntVect toIntVect() const { return value; } };
struct Box {
  IntVect lo, hi;
  IntVect smallEnd() const { return lo; }
  IntVect bigEnd() const { return hi; }
};
class BoxArray {
public:
  std::vector<Box> boxes;
  IndexType type{{{{0, 0, 0}}}};
  int size() const { return static_cast<int>(boxes.size()); }
  const Box &operator[](int i) const { return boxes.at(static_cast<std::size_t>(i)); }
  IndexType ixType() const { return type; }
};
class DistributionMapping {
public:
  std::vector<int> processors;
  const std::vector<int> &ProcessorMap() const { return processors; }
};
struct TileRecord { int index; Box valid, tile; };
class FabArrayBase { public: std::vector<TileRecord> tiles; };
class MFItInfo { public: MFItInfo &EnableTiling() { return *this; } };
class MFIter {
public:
  MFIter(const FabArrayBase &fab, const MFItInfo &) : fab_(&fab) {}
  bool isValid() const { return pos_ < fab_->tiles.size(); }
  MFIter &operator++() { ++pos_; return *this; }
  int index() const { return fab_->tiles.at(pos_).index; }
  Box validbox() const { return fab_->tiles.at(pos_).valid; }
  Box tilebox() const { return fab_->tiles.at(pos_).tile; }
private:
  const FabArrayBase *fab_;
  std::size_t pos_{0};
};
template <typename T> class Array4 {
public:
  Array4(T *data, Box allocation, int components)
      : data_(data), box_(allocation), components_(components) {}
  T *ptr(int i, int j, int k, int component) const {
    const auto lo = box_.smallEnd();
    const auto hi = box_.bigEnd();
    if (i < lo[0] || i > hi[0] || j < lo[1] || j > hi[1] ||
        k < lo[2] || k > hi[2] || component < 0 || component >= components_)
      return nullptr;
    const std::size_t nx = static_cast<std::size_t>(hi[0] - lo[0] + 1);
    const std::size_t ny = static_cast<std::size_t>(hi[1] - lo[1] + 1);
    const std::size_t nz = static_cast<std::size_t>(hi[2] - lo[2] + 1);
    const std::size_t point = static_cast<std::size_t>(i - lo[0]) +
        nx * (static_cast<std::size_t>(j - lo[1]) +
              ny * static_cast<std::size_t>(k - lo[2]));
    return data_ + point + nx * ny * nz * static_cast<std::size_t>(component);
  }
private:
  T *data_;
  Box box_;
  int components_;
};
class MultiFab {
public:
  MultiFab(BoxArray ba, DistributionMapping dm, IntVect grow, int components,
           bool defined = true, bool alias_copy = false)
      : ba_(std::move(ba)), dm_(std::move(dm)), grow_(grow),
        components_(components), defined_(defined), alias_copy_(alias_copy),
        data_(std::make_shared<std::vector<std::vector<double>>>()) {
    data_->resize(static_cast<std::size_t>(ba_.size()));
    for (int i = 0; i < ba_.size(); ++i) {
      const auto lo = ba_[i].smallEnd(); const auto hi = ba_[i].bigEnd();
      std::size_t points = 1;
      for (int d = 0; d < 3; ++d)
        points *= static_cast<std::size_t>(hi[d] - lo[d] + 1 + 2 * grow_[d]);
      (*data_)[static_cast<std::size_t>(i)].assign(
          points * static_cast<std::size_t>(components_), 1000.0 + i);
    }
  }
  MultiFab deepCopy() const {
    MultiFab result(ba_, dm_, grow_, components_, defined_, false);
    if (alias_copy_) result.data_ = data_;
    else *result.data_ = *data_;
    return result;
  }
  bool isDefined() const { return defined_; }
  int nComp() const { return components_; }
  const BoxArray &boxArray() const { return ba_; }
  const DistributionMapping &DistributionMap() const { return dm_; }
  IntVect nGrowVect() const { return grow_; }
  Array4<double> array(int index) {
    auto box = ba_[index];
    for (int d = 0; d < 3; ++d) { box.lo.v[d] -= grow_[d]; box.hi.v[d] += grow_[d]; }
    return Array4<double>((*data_)[static_cast<std::size_t>(index)].data(), box,
                          components_);
  }
  void set_defined(bool value) { defined_ = value; }
  void set_components(int value) { components_ = value; }
  void set_grow(IntVect value) { grow_ = value; }
  void set_processors(std::vector<int> value) { dm_.processors = std::move(value); }
  void set_index_type(IntVect value) { ba_.type.value = value; }
  void set_box_hi(int index, IntVect value) { ba_.boxes.at(index).hi = value; }
private:
  BoxArray ba_;
  DistributionMapping dm_;
  IntVect grow_;
  int components_;
  bool defined_;
  bool alias_copy_;
  std::shared_ptr<std::vector<std::vector<double>>> data_;
};
} // namespace amrex

using CCTK_REAL = double;
constexpr int CCTK_GF = 1;
struct cGroup { int grouptype; int numvars; };
struct cGH {
  int cctk_dim{}, cctk_iteration{};
  int *cctk_gsh{}, *cctk_lsh{}, *cctk_lbnd{}, *cctk_ubnd{};
  int cctk_level{}, cctk_patch{}, cctk_npatches{}, cctk_component{};
  int *cctk_tile_min{}, *cctk_tile_max{}, *cctk_ash{};
  int cctk_alignment{}, cctk_alignment_offset{};
  int *cctk_to{}, *cctk_from{};
  CCTK_REAL cctk_delta_time{};
  CCTK_REAL *cctk_delta_space{}, *cctk_origin_space{};
  int *cctk_bbox{}, *cctk_levfac{}, *cctk_levoff{}, *cctk_levoffdenom{};
  int cctk_timefac{}, cctk_convlevel{}, cctk_convfac{};
  int *cctk_nghostzones{};
  CCTK_REAL cctk_time{};
  const void *current_scheduled_function{};
  const char *identity{};
  void ***data{};
  void **extensions{};
  void *GroupData{};
};
static std::int64_t current_epoch = 7;
static std::array<int, 8> declared_tls{{1, 2, 2, 1, 1, 2, 2, 1}};
int CCTK_NumVars() { return static_cast<int>(declared_tls.size()); }
int CCTK_DeclaredTimeLevelsVI(int vi) { return declared_tls.at(vi); }
int CCTK_FirstVarIndexI(int gi) { return gi == 0 ? 1 : (gi == 2 ? 5 : -1); }
int CCTK_GroupData(int gi, cGroup *group) {
  if ((gi == 0 || gi == 2) && group != nullptr) {
    group->grouptype = CCTK_GF; group->numvars = 2; return 0;
  }
  return 1;
}
extern "C" std::int64_t CarpetX_GetEpoch() { return current_epoch; }

namespace CarpetX {
inline constexpr int dim = 3;
struct TemplateGeometry {
  std::array<int, 3> gsh{{20, 20, 20}}, lsh{{8, 8, 8}}, lbnd{{0, 0, 0}},
      ubnd{{7, 7, 7}}, tile_min{{0, 0, 0}}, tile_max{{3, 7, 7}},
      ash{{12, 12, 12}}, to{{0, 0, 0}}, from{{0, 0, 0}},
      levfac{{2, 2, 2}}, levoff{{0, 0, 0}}, levoffdenom{{1, 1, 1}},
      nghostzones{{2, 2, 2}};
  std::array<int, 6> bbox{{1, 1, 1, 1, 1, 1}};
  std::array<double, 3> delta{{0.25, 0.25, 0.25}}, origin{{0, 0, 0}};
};
class GHExt {
public:
  struct LocalGHPtr {
    LocalGHPtr() = default;
    LocalGHPtr(std::unique_ptr<cGH> value) : value_(std::move(value)) {}
    LocalGHPtr(LocalGHPtr &&) noexcept = default;
    LocalGHPtr &operator=(LocalGHPtr &&) noexcept = default;
    LocalGHPtr(const LocalGHPtr &) = delete;
    LocalGHPtr &operator=(const LocalGHPtr &) = delete;
    cGH *get() const { return value_.get(); }
  private:
    std::unique_ptr<cGH> value_;
  };
  struct PatchData {
    struct LevelData {
      struct GroupData {
        int groupindex{}, firstvarindex{}, numvars{}, patch{}, level{};
        std::array<int, 3> indextype{{0, 0, 0}}, nghostzones{{1, 1, 1}};
        std::vector<std::unique_ptr<amrex::MultiFab>> mfab;
        std::vector<std::vector<int>> valid;
      };
      int patch{}, level{};
      std::unique_ptr<amrex::FabArrayBase> fab;
      std::vector<LocalGHPtr> local_cctkGHs;
      std::vector<std::unique_ptr<TemplateGeometry>> geometry;
      std::vector<std::unique_ptr<GroupData>> groupdata;
    };
    std::vector<LevelData> leveldata;
  };
  std::vector<PatchData> patchdata;
};
struct MFPointer {
  explicit MFPointer(const amrex::MFIter &mfi)
      : index_(mfi.index()), valid_(mfi.validbox()), tile_(mfi.tilebox()) {}
  int index() const { return index_; }
  amrex::Box validbox() const { return valid_; }
  amrex::Box tilebox() const { return tile_; }
  int index_; amrex::Box valid_, tile_;
};
struct GridPtrDesc1 {
  GridPtrDesc1(const GHExt::PatchData::LevelData &,
               const GHExt::PatchData::LevelData::GroupData &group,
               const MFPointer &mfp) : mfp_(mfp), group_(group) {}
  template <typename T> T *ptr(const amrex::Array4<T> &vars, int vi) const {
    const auto lo = mfp_.validbox().smallEnd();
    return vars.ptr(lo[0] - group_.nghostzones[0],
                    lo[1] - group_.nghostzones[1],
                    lo[2] - group_.nghostzones[2], vi);
  }
  MFPointer mfp_;
  const GHExt::PatchData::LevelData::GroupData &group_;
};
} // namespace CarpetX

#define CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
#define CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_STANDALONE
#define private public
#include "../src/subcycling_scratch_state.cxx"
#undef private
#include "../src/subcycling_scratch_adapter.hxx"

namespace CarpetX {
CertifiedScratchLevelFrame::CertifiedScratchLevelFrame(int patch,
                                                       ScratchLevelFrame frame)
    : patch_(patch), frame_(std::move(frame)) {}
CertifiedScratchLevelFrame::CertifiedScratchLevelFrame(
    CertifiedScratchLevelFrame &&other) noexcept
    : patch_(other.patch_), frame_(std::move(other.frame_)) { other.patch_ = -1; }
CertifiedScratchLevelFrame &CertifiedScratchLevelFrame::operator=(
    CertifiedScratchLevelFrame &&other) noexcept {
  if (this != &other) { patch_ = other.patch_; frame_ = std::move(other.frame_); other.patch_ = -1; }
  return *this;
}
int CertifiedScratchLevelFrame::patch() const noexcept { return patch_; }
const ScratchLevelFrame &CertifiedScratchLevelFrame::frame() const noexcept { return frame_; }
} // namespace CarpetX

#include "../src/subcycling_scratch_local_gh.cxx"

namespace CarpetX {
class ScratchLocalLevelBindingTestAccess {
public:
  static CertifiedScratchLevelFrame certify(int patch, ScratchLevelFrame frame) {
    return CertifiedScratchLevelFrame(patch, std::move(frame));
  }
  static ScratchLevelFrame &frame(CertifiedScratchLevelFrame &certified) {
    return certified.frame_;
  }
  static ScratchLevelFrame &frame(ScratchLocalLevelBinding &binding) {
    return binding.frame_.frame_;
  }
  static const cGH &gh(const ScratchLocalLevelBinding &binding, std::size_t i) {
    return binding.contexts_.at(i)->gh;
  }
  static int live_storage_count() {
    return ScratchLocalLevelBinding::Storage::live_count;
  }
};
} // namespace CarpetX

namespace {
using namespace CarpetX;
amrex::Box box(std::array<int, 3> lo, std::array<int, 3> hi) {
  return amrex::Box{amrex::IntVect{lo}, amrex::IntVect{hi}};
}
amrex::BoxArray boxes() {
  amrex::BoxArray ba;
  for (int i = 0; i < 6; ++i) ba.boxes.push_back(box({i * 8, 0, 0}, {i * 8 + 7, 7, 7}));
  return ba;
}
std::unique_ptr<cGH> make_template(TemplateGeometry &g, int component,
                                   std::array<int, 3> tile_lo,
                                   std::array<int, 3> tile_hi) {
  g.tile_min = tile_lo; g.tile_max = tile_hi;
  auto gh = std::make_unique<cGH>();
  gh->cctk_dim = 3; gh->cctk_level = 0; gh->cctk_patch = 0;
  gh->cctk_npatches = 1; gh->cctk_component = component;
  gh->cctk_alignment = 64; gh->cctk_alignment_offset = 0;
  gh->cctk_convlevel = 2; gh->cctk_convfac = 4;
  gh->cctk_gsh=g.gsh.data(); gh->cctk_lsh=g.lsh.data(); gh->cctk_lbnd=g.lbnd.data();
  gh->cctk_ubnd=g.ubnd.data(); gh->cctk_tile_min=g.tile_min.data();
  gh->cctk_tile_max=g.tile_max.data(); gh->cctk_ash=g.ash.data();
  gh->cctk_to=g.to.data(); gh->cctk_from=g.from.data(); gh->cctk_delta_space=g.delta.data();
  gh->cctk_origin_space=g.origin.data(); gh->cctk_bbox=g.bbox.data();
  gh->cctk_levfac=g.levfac.data(); gh->cctk_levoff=g.levoff.data();
  gh->cctk_levoffdenom=g.levoffdenom.data(); gh->cctk_nghostzones=g.nghostzones.data();
  gh->identity="live"; gh->extensions=reinterpret_cast<void **>(0x1);
  gh->GroupData=reinterpret_cast<void *>(0x2); gh->current_scheduled_function=reinterpret_cast<void *>(0x3);
  return gh;
}
GHExt make_ghext(bool alias_copy = false) {
  GHExt result; result.patchdata.resize(1); result.patchdata[0].leveldata.resize(1);
  auto &level = result.patchdata[0].leveldata[0]; level.patch=0; level.level=0;
  level.fab=std::make_unique<amrex::FabArrayBase>();
  level.fab->tiles = {{2, box({16,0,0},{23,7,7}), box({16,0,0},{19,7,7})},
                      {2, box({16,0,0},{23,7,7}), box({20,0,0},{23,7,7})},
                      {5, box({40,0,0},{47,7,7}), box({40,0,0},{47,7,7})}};
  for (const auto &tile : level.fab->tiles) {
    auto geometry=std::make_unique<TemplateGeometry>();
    auto gh=make_template(*geometry, tile.index, tile.tile.lo.v, tile.tile.hi.v);
    level.geometry.push_back(std::move(geometry)); level.local_cctkGHs.push_back(std::move(gh));
  }
  level.groupdata.resize(3);
  const amrex::DistributionMapping dm{{0,0,0,0,0,0}}; const amrex::IntVect grow{{2,2,2}};
  for (int gi : {0,2}) {
    auto group=std::make_unique<GHExt::PatchData::LevelData::GroupData>();
    group->groupindex=gi; group->firstvarindex=(gi==0?1:5); group->numvars=2;
    group->patch=0; group->level=0; group->valid={{1,1}};
    group->mfab.push_back(std::make_unique<amrex::MultiFab>(boxes(),dm,grow,2,true,alias_copy));
    level.groupdata[static_cast<std::size_t>(gi)]=std::move(group);
  }
  return result;
}
CertifiedScratchLevelFrame make_certified(const GHExt &ghext) {
  const auto &level=ghext.patchdata[0].leveldata[0];
  UncertifiedScratchLevelManifest manifest{current_epoch,0,{}};
  std::vector<ScratchGFView> views; std::vector<std::vector<ScratchValidity>> validity(2);
  int ordinal=0;
  for (int gi : {0,2}) {
    auto &group=*level.groupdata[static_cast<std::size_t>(gi)];
    validity[ordinal].assign(2,ScratchValidity{true,true,true});
    ScratchGFKey key{current_epoch,0,0,gi}; const void *identity=group.mfab[0].get();
    manifest.entries.push_back({key,identity});
    views.push_back({key,0,identity,group.mfab[0].get(),&validity[ordinal]}); ++ordinal;
  }
  return ScratchLocalLevelBindingTestAccess::certify(
      0,ScratchLevelFrame::copy_tl0(manifest,views));
}
template <typename F> void expect_reject(F &&f) {
  bool rejected=false; try { f(); } catch (const std::invalid_argument &) { rejected=true; }
  assert(rejected);
}
void test_two_noncontiguous_groups_and_complete_null_table() {
  auto ghext=make_ghext(); auto certified=make_certified(ghext);
  auto binding=ScratchLocalLevelBinding::bind(ghext,std::move(certified));
  assert(binding->local_context_count()==3); const auto &gh=ScratchLocalLevelBindingTestAccess::gh(*binding,0);
  for (int vi=0;vi<CCTK_NumVars();++vi) for(int tl=0;tl<CCTK_DeclaredTimeLevelsVI(vi);++tl)
    assert(((vi==1||vi==2||vi==5||vi==6)&&tl==0) ? gh.data[vi][tl]!=nullptr : gh.data[vi][tl]==nullptr);
}
void test_local_ordinal_is_not_fab_index_and_repeated_tiles() {
  auto ghext=make_ghext(); auto certified=make_certified(ghext);
  auto binding=ScratchLocalLevelBinding::bind(ghext,std::move(certified));
  assert(ScratchLocalLevelBindingTestAccess::gh(*binding,0).cctk_component==2);
  assert(ScratchLocalLevelBindingTestAccess::gh(*binding,1).cctk_component==2);
  assert(ScratchLocalLevelBindingTestAccess::gh(*binding,2).cctk_component==5);
  assert(ScratchLocalLevelBindingTestAccess::gh(*binding,0).cctk_tile_min[0]!=
         ScratchLocalLevelBindingTestAccess::gh(*binding,1).cctk_tile_min[0]);
}
void test_grid_ptr_mapping_and_scratch_isolation() {
  auto ghext=make_ghext(); auto certified=make_certified(ghext);
  auto binding=ScratchLocalLevelBinding::bind(ghext,std::move(certified));
  auto &level=ghext.patchdata[0].leveldata[0];
  auto &scratch_frame=ScratchLocalLevelBindingTestAccess::frame(*binding);
  std::size_t ordinal=0;
  const auto mfitinfo=amrex::MFItInfo().EnableTiling();
  for (amrex::MFIter mfi(*level.fab,mfitinfo);mfi.isValid();++mfi,++ordinal) {
    const auto &gh=ScratchLocalLevelBindingTestAccess::gh(*binding,ordinal);
    std::size_t frame_index=0;
    for (int gi : {0,2}) {
      auto &group=*level.groupdata[static_cast<std::size_t>(gi)];
      auto scratch_vars=scratch_frame.mutable_multifab(frame_index).array(mfi.index());
      auto live_vars=group.mfab[0]->array(mfi.index());
      GridPtrDesc1 grid(level,group,MFPointer(mfi));
      for(int vi=0;vi<group.numvars;++vi) {
        auto *expected=grid.ptr(scratch_vars,vi);
        auto *live=grid.ptr(live_vars,vi);
        assert(gh.data[group.firstvarindex+vi][0]==expected);
        assert(expected!=nullptr && expected!=live);
      }
      ++frame_index;
    }
  }
  auto *live=level.groupdata[0]->mfab[0]->array(2).ptr(15,-1,-1,0);
  auto *scratch=static_cast<double *>(
      ScratchLocalLevelBindingTestAccess::gh(*binding,0).data[1][0]);
  const double before=*live; *scratch=42.0; assert(*live==before);
}
void test_geometry_is_deep_owned_and_aliases_are_null() {
  auto ghext=make_ghext(); auto certified=make_certified(ghext);
  auto binding=ScratchLocalLevelBinding::bind(ghext,std::move(certified));
  const auto &gh=ScratchLocalLevelBindingTestAccess::gh(*binding,0); int copied=gh.cctk_gsh[0];
  ghext.patchdata[0].leveldata[0].local_cctkGHs.at(0).get()->cctk_gsh[0]=999;
  assert(gh.cctk_gsh[0]==copied && gh.identity==nullptr && gh.extensions==nullptr &&
         gh.GroupData==nullptr && gh.current_scheduled_function==nullptr);
  assert(gh.cctk_iteration==-1 && gh.cctk_timefac==0 && std::isnan(gh.cctk_time) && std::isnan(gh.cctk_delta_time));
}
void test_stale_and_incomplete_frames_reject() {
  { auto g=make_ghext(); auto c=make_certified(g); ++current_epoch;
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); --current_epoch; }
  { auto g=make_ghext(); auto c=make_certified(g);
    ScratchLocalLevelBindingTestAccess::frame(c).entries_.pop_back();
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
}
void test_duplicate_and_extra_frames_reject() {
  { auto g=make_ghext(); auto c=make_certified(g); auto &frame=ScratchLocalLevelBindingTestAccess::frame(c);
    frame.entries_[1].key=frame.entries_[0].key;
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); auto &frame=ScratchLocalLevelBindingTestAccess::frame(c);
    const auto &source=frame.entries_.back();
    auto owned=std::make_unique<amrex::MultiFab>(source.multifab->deepCopy());
    frame.entries_.push_back(ScratchLevelFrame::Entry{
        ScratchGFKey{current_epoch,0,0,3},std::move(owned),source.validity});
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
}
void test_group_ranges_and_time_levels_reject() {
  { auto g=make_ghext(); auto c=make_certified(g); g.patchdata[0].leveldata[0].groupdata[2]->firstvarindex=2;
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); declared_tls[5]=0;
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); declared_tls[5]=2; }
}
void test_layout_distribution_grow_component_and_validity_reject() {
  { auto g=make_ghext(); auto c=make_certified(g); ScratchLocalLevelBindingTestAccess::frame(c).mutable_multifab(0).set_grow({{{3,2,2}}});
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); g.patchdata[0].leveldata[0].groupdata[0]->valid[0].pop_back();
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); ScratchLocalLevelBindingTestAccess::frame(c).mutable_validity(0).pop_back();
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); ScratchLocalLevelBindingTestAccess::frame(c).mutable_multifab(0).set_processors({0,0,1,0,0,0});
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); ScratchLocalLevelBindingTestAccess::frame(c).mutable_multifab(0).set_box_hi(0,{{{8,7,7}}});
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); ScratchLocalLevelBindingTestAccess::frame(c).mutable_multifab(0).set_index_type({{{1,0,0}}});
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); ScratchLocalLevelBindingTestAccess::frame(c).mutable_multifab(0).set_components(3);
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(); auto c=make_certified(g); ScratchLocalLevelBindingTestAccess::frame(c).mutable_multifab(0).set_defined(false);
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
}
void test_cached_context_mismatch_and_pointer_alias_reject() {
  { auto g=make_ghext(); auto c=make_certified(g); g.patchdata[0].leveldata[0].local_cctkGHs.at(0).get()->cctk_component=0;
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
  { auto g=make_ghext(true); auto c=make_certified(g);
    expect_reject([&]{ (void)ScratchLocalLevelBinding::bind(g,std::move(c)); }); }
}
void test_later_context_failure_destroys_partial_storage_without_publication() {
  assert(ScratchLocalLevelBindingTestAccess::live_storage_count()==0);
  auto g=make_ghext(); auto c=make_certified(g);
  g.patchdata[0].leveldata[0].local_cctkGHs.at(1).get()->cctk_gsh=nullptr;
  std::unique_ptr<ScratchLocalLevelBinding> published;
  expect_reject([&]{ published=ScratchLocalLevelBinding::bind(g,std::move(c)); });
  assert(!published);
  assert(ScratchLocalLevelBindingTestAccess::live_storage_count()==0);
}
} // namespace

int main() {
  static_assert(!std::is_copy_constructible_v<ScratchLocalLevelBinding>);
  static_assert(!std::is_move_constructible_v<ScratchLocalLevelBinding>);
  test_two_noncontiguous_groups_and_complete_null_table();
  test_local_ordinal_is_not_fab_index_and_repeated_tiles();
  test_grid_ptr_mapping_and_scratch_isolation();
  test_geometry_is_deep_owned_and_aliases_are_null();
  test_stale_and_incomplete_frames_reject();
  test_duplicate_and_extra_frames_reject();
  test_group_ranges_and_time_levels_reject();
  test_layout_distribution_grow_component_and_validity_reject();
  test_cached_context_mismatch_and_pointer_alias_reject();
  test_later_context_failure_destroys_partial_storage_without_publication();
  std::cout << "Phase 8B2 scratch local cGH tests passed\n";
}
