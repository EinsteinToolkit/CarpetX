#include "subcycling_scratch_local_gh.hxx"

#ifndef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_STANDALONE
#include "driver.hxx"
#include "schedule.hxx"
#include <cctk.h>
#endif

#include <limits>
#include <stdexcept>
#include <utility>

namespace CarpetX {
namespace {

[[noreturn]] void reject(const char *message) {
  throw std::invalid_argument(message);
}

template <typename T>
void copy_required(std::vector<T> &dst, const T *src, std::size_t count) {
  if (src == nullptr)
    reject("scratch local cGH geometry source is null");
  dst.assign(src, src + count);
}

bool same_intvect(const amrex::IntVect &a,
                  const amrex::IntVect &b) noexcept {
  for (int d = 0; d < dim; ++d)
    if (a[d] != b[d])
      return false;
  return true;
}

bool same_box(const amrex::Box &a, const amrex::Box &b) noexcept {
  return same_intvect(a.smallEnd(), b.smallEnd()) &&
         same_intvect(a.bigEnd(), b.bigEnd());
}

bool same_layout(const amrex::MultiFab &a, const amrex::MultiFab &b) {
  if (a.nComp() != b.nComp() ||
      !same_intvect(a.nGrowVect(), b.nGrowVect()))
    return false;
  const auto &aba = a.boxArray();
  const auto &abb = b.boxArray();
  if (aba.size() != abb.size() ||
      !same_intvect(aba.ixType().toIntVect(), abb.ixType().toIntVect()))
    return false;
  for (int i = 0; i < aba.size(); ++i)
    if (!same_box(aba[i], abb[i]))
      return false;
  return a.DistributionMap().ProcessorMap() ==
         b.DistributionMap().ProcessorMap();
}

struct ValidatedEntry {
  std::size_t frame_index;
  int group_index;
  int first_var;
  int num_vars;
};

struct ValidatedBinding {
  const GHExt::PatchData::LevelData *leveldata;
  std::vector<ValidatedEntry> entries;
  std::size_t local_context_count;
};

ValidatedBinding validate_binding(const GHExt &ghext,
                                  const CertifiedScratchLevelFrame &certified) {
  const int patch = certified.patch();
  const auto &frame = certified.frame();
  const int level = frame.level();
  if (patch < 0 || static_cast<std::size_t>(patch) >= ghext.patchdata.size())
    reject("scratch local cGH patch is out of range");
  const auto &patchdata = ghext.patchdata[static_cast<std::size_t>(patch)];
  if (level < 0 ||
      static_cast<std::size_t>(level) >= patchdata.leveldata.size())
    reject("scratch local cGH level is out of range");
  if (frame.hierarchy_epoch() < 0 ||
      frame.hierarchy_epoch() !=
          static_cast<std::int64_t>(CarpetX_GetEpoch()))
    reject("scratch local cGH frame epoch is stale");
  const auto &leveldata =
      patchdata.leveldata[static_cast<std::size_t>(level)];
  if (leveldata.patch != patch || leveldata.level != level ||
      leveldata.fab == nullptr)
    reject("scratch local cGH level metadata is inconsistent");

  const int total_vars = CCTK_NumVars();
  if (total_vars < 0)
    reject("scratch local cGH variable count is invalid");
  std::vector<bool> claimed(static_cast<std::size_t>(total_vars), false);
  std::vector<ValidatedEntry> entries;
  entries.reserve(frame.entry_count());
  std::size_t frame_index = 0;
  for (std::size_t gi = 0; gi < leveldata.groupdata.size(); ++gi) {
    const auto &owner = leveldata.groupdata[gi];
    if (owner == nullptr)
      continue;
    if (frame_index >= frame.entry_count())
      reject("scratch local cGH frame omits a canonical GF group");
    const auto &groupdata = *owner;
    const auto &key = frame.key(frame_index);
    if (gi > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        key.hierarchy_epoch != frame.hierarchy_epoch() ||
        key.patch != patch || key.level != level ||
        key.group_index != static_cast<int>(gi) ||
        groupdata.groupindex != static_cast<int>(gi) ||
        groupdata.patch != patch || groupdata.level != level)
      reject("scratch local cGH frame key or group metadata is inconsistent");
    cGroup group;
    if (CCTK_GroupData(static_cast<int>(gi), &group) != 0 ||
        group.grouptype != CCTK_GF || group.numvars != groupdata.numvars ||
        CCTK_FirstVarIndexI(static_cast<int>(gi)) !=
            groupdata.firstvarindex)
      reject("scratch local cGH canonical entry is not a grid function");
    if (groupdata.firstvarindex < 0 || groupdata.numvars <= 0 ||
        groupdata.numvars > total_vars ||
        groupdata.firstvarindex > total_vars - groupdata.numvars)
      reject("scratch local cGH group variable range is invalid");
    for (int vi = 0; vi < groupdata.numvars; ++vi) {
      const int var = groupdata.firstvarindex + vi;
      if (claimed[static_cast<std::size_t>(var)])
        reject("scratch local cGH group variable ranges overlap");
      if (CCTK_DeclaredTimeLevelsVI(var) <= 0)
        reject("scratch local cGH variable has no declared time level");
      claimed[static_cast<std::size_t>(var)] = true;
    }
    if (groupdata.mfab.empty() || groupdata.valid.empty() ||
        groupdata.mfab[0] == nullptr)
      reject("scratch local cGH canonical TL0 is unavailable");
    const auto &live = *groupdata.mfab[0];
    const auto &scratch = frame.multifab(frame_index);
    if (!live.isDefined() || !scratch.isDefined() ||
        live.nComp() != groupdata.numvars ||
        scratch.nComp() != groupdata.numvars ||
        groupdata.valid[0].size() !=
            static_cast<std::size_t>(groupdata.numvars) ||
        frame.validity(frame_index).size() !=
            static_cast<std::size_t>(groupdata.numvars) ||
        !same_layout(live, scratch))
      reject("scratch local cGH TL0 layout or validity is inconsistent");
    entries.push_back(ValidatedEntry{frame_index, static_cast<int>(gi),
                                     groupdata.firstvarindex,
                                     groupdata.numvars});
    ++frame_index;
  }
  if (entries.empty() || frame_index != frame.entry_count())
    reject("scratch local cGH frame has an extra or empty canonical set");

  std::size_t local_count = 0;
  const auto mfitinfo = amrex::MFItInfo().EnableTiling();
  for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {
    if (local_count >= leveldata.local_cctkGHs.size())
      reject("scratch local cGH cache omits a tiled context");
    const cGH *local = leveldata.local_cctkGHs[local_count].get();
    if (local == nullptr || local->cctk_patch != patch ||
        local->cctk_level != level || local->cctk_component != mfi.index())
      reject("scratch local cGH cached context does not match MFIter");
    ++local_count;
  }
  if (local_count != leveldata.local_cctkGHs.size())
    reject("scratch local cGH cache contains extra tiled contexts");
  return ValidatedBinding{&leveldata, std::move(entries), local_count};
}

} // namespace

struct ScratchLocalLevelBinding::Storage {
#ifdef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
  inline static int live_count = 0;
#endif
  cGH gh{};
  std::vector<int> gsh, lsh, lbnd, ubnd, tile_min, tile_max, ash;
  std::vector<int> to, from, bbox, levfac, levoff, levoffdenom, nghostzones;
  std::vector<CCTK_REAL> delta_space, origin_space;
  std::vector<std::vector<void *>> data_rows;
  std::vector<void **> data_row_pointers;

  Storage(const cGH &src, const int total_vars) {
    if (src.cctk_dim != dim || total_vars < 0)
      reject("scratch local cGH source dimension is invalid");
    gh.cctk_dim = src.cctk_dim;
    gh.cctk_level = src.cctk_level;
    gh.cctk_patch = src.cctk_patch;
    gh.cctk_npatches = src.cctk_npatches;
    gh.cctk_component = src.cctk_component;
    gh.cctk_alignment = src.cctk_alignment;
    gh.cctk_alignment_offset = src.cctk_alignment_offset;
    gh.cctk_convlevel = src.cctk_convlevel;
    gh.cctk_convfac = src.cctk_convfac;
    gh.cctk_iteration = -1;
    gh.cctk_delta_time = std::numeric_limits<CCTK_REAL>::quiet_NaN();
    gh.cctk_time = std::numeric_limits<CCTK_REAL>::quiet_NaN();
    gh.cctk_timefac = 0;
    gh.current_scheduled_function = nullptr;
    gh.identity = nullptr;
    gh.extensions = nullptr;
    gh.GroupData = nullptr;

    copy_required(gsh, src.cctk_gsh, dim);
    copy_required(lsh, src.cctk_lsh, dim);
    copy_required(lbnd, src.cctk_lbnd, dim);
    copy_required(ubnd, src.cctk_ubnd, dim);
    copy_required(tile_min, src.cctk_tile_min, dim);
    copy_required(tile_max, src.cctk_tile_max, dim);
    copy_required(ash, src.cctk_ash, dim);
    copy_required(to, src.cctk_to, dim);
    copy_required(from, src.cctk_from, dim);
    copy_required(delta_space, src.cctk_delta_space, dim);
    copy_required(origin_space, src.cctk_origin_space, dim);
    copy_required(bbox, src.cctk_bbox, 2 * dim);
    copy_required(levfac, src.cctk_levfac, dim);
    copy_required(levoff, src.cctk_levoff, dim);
    copy_required(levoffdenom, src.cctk_levoffdenom, dim);
    copy_required(nghostzones, src.cctk_nghostzones, dim);
    gh.cctk_gsh = gsh.data(); gh.cctk_lsh = lsh.data();
    gh.cctk_lbnd = lbnd.data(); gh.cctk_ubnd = ubnd.data();
    gh.cctk_tile_min = tile_min.data(); gh.cctk_tile_max = tile_max.data();
    gh.cctk_ash = ash.data(); gh.cctk_to = to.data(); gh.cctk_from = from.data();
    gh.cctk_delta_space = delta_space.data();
    gh.cctk_origin_space = origin_space.data(); gh.cctk_bbox = bbox.data();
    gh.cctk_levfac = levfac.data(); gh.cctk_levoff = levoff.data();
    gh.cctk_levoffdenom = levoffdenom.data();
    gh.cctk_nghostzones = nghostzones.data();

    data_rows.resize(static_cast<std::size_t>(total_vars));
    data_row_pointers.resize(static_cast<std::size_t>(total_vars), nullptr);
    for (int var = 0; var < total_vars; ++var) {
      const int tls = CCTK_DeclaredTimeLevelsVI(var);
      if (tls <= 0)
        reject("scratch local cGH variable has no declared time level");
      auto &row = data_rows[static_cast<std::size_t>(var)];
      row.assign(static_cast<std::size_t>(tls), nullptr);
      data_row_pointers[static_cast<std::size_t>(var)] = row.data();
    }
    gh.data = data_row_pointers.data();
#ifdef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
    ++live_count;
#endif
  }

  ~Storage() {
#ifdef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
    --live_count;
#endif
  }
};

std::unique_ptr<ScratchLocalLevelBinding>
ScratchLocalLevelBinding::bind(const GHExt &ghext,
                               CertifiedScratchLevelFrame &&frame) {
  (void)validate_binding(ghext, frame);
  std::unique_ptr<ScratchLocalLevelBinding> result(
      new ScratchLocalLevelBinding(std::move(frame)));
  result->build_contexts(ghext);
  return result;
}

ScratchLocalLevelBinding::ScratchLocalLevelBinding(
    CertifiedScratchLevelFrame &&frame)
    : frame_(std::move(frame)) {}
ScratchLocalLevelBinding::~ScratchLocalLevelBinding() = default;

void ScratchLocalLevelBinding::build_contexts(const GHExt &ghext) {
  const auto validated = validate_binding(ghext, frame_);
  const auto &leveldata = *validated.leveldata;
  contexts_.reserve(validated.local_context_count);
  std::size_t local_ordinal = 0;
  const auto mfitinfo = amrex::MFItInfo().EnableTiling();
  for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {
    cGH *const source = leveldata.local_cctkGHs.at(local_ordinal).get();
    if (source == nullptr)
      reject("scratch local cGH cached context is null during binding");
    auto storage = std::make_unique<Storage>(*source, CCTK_NumVars());
    for (const auto &entry : validated.entries) {
      const auto &groupdata =
          *leveldata.groupdata[static_cast<std::size_t>(entry.group_index)];
      auto scratch_vars =
          frame_.frame_.mutable_multifab(entry.frame_index).array(mfi.index());
      auto live_vars = groupdata.mfab[0]->array(mfi.index());
      GridPtrDesc1 grid(leveldata, groupdata, MFPointer(mfi));
      for (int vi = 0; vi < entry.num_vars; ++vi) {
        void *scratch_ptr = grid.ptr(scratch_vars, vi);
        void *live_ptr = grid.ptr(live_vars, vi);
        if (scratch_ptr == nullptr || scratch_ptr == live_ptr)
          reject("scratch local cGH pointer aliases canonical TL0");
        storage->gh.data[entry.first_var + vi][0] = scratch_ptr;
      }
    }
    contexts_.push_back(std::move(storage));
    ++local_ordinal;
  }
  if (frame_.frame_.hierarchy_epoch() !=
      static_cast<std::int64_t>(CarpetX_GetEpoch()))
    reject("hierarchy epoch changed while building scratch local cGH views");
}

int ScratchLocalLevelBinding::patch() const noexcept { return frame_.patch_; }
int ScratchLocalLevelBinding::level() const noexcept {
  return frame_.frame_.level();
}
std::int64_t ScratchLocalLevelBinding::hierarchy_epoch() const noexcept {
  return frame_.frame_.hierarchy_epoch();
}
std::size_t ScratchLocalLevelBinding::local_context_count() const noexcept {
  return contexts_.size();
}

cGH *ScratchLocalLevelBinding::context_for_executor(
    const std::size_t ordinal) noexcept {
  return ordinal < contexts_.size() ? &contexts_[ordinal]->gh : nullptr;
}

ScratchLevelFrame &ScratchLocalLevelBinding::frame_for_executor() noexcept {
  return frame_.frame_;
}

int ScratchLocalLevelBinding::level_for_executor() const noexcept {
  return frame_.frame_.level();
}

void ScratchLocalLevelBinding::claim_for_executor() {
  if (executor_faulted_.load())
    throw std::logic_error("scratch binding is faulted");
  bool expected = false;
  if (!executor_busy_.compare_exchange_strong(expected, true))
    throw std::logic_error("scratch binding is busy");
  if (executor_faulted_.load()) {
    executor_busy_.store(false);
    throw std::logic_error("scratch binding is faulted");
  }
}

void ScratchLocalLevelBinding::release_for_executor() noexcept {
  executor_busy_.store(false);
}

void ScratchLocalLevelBinding::fault_for_executor() noexcept {
  executor_faulted_.store(true);
}

} // namespace CarpetX

#ifndef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_STANDALONE
#include "subcycling_scratch_schedule_executor.hxx"

namespace CarpetX {

void ScratchLocalLevelBinding::install_stage_coordinates_for_executor(
    const ScratchStageCoordinates &coordinates) noexcept {
  for (auto &context : contexts_) {
    context->gh.cctk_iteration = coordinates.level_iteration;
    context->gh.cctk_time = coordinates.stage_time;
    context->gh.cctk_delta_time = coordinates.base_delta_time;
    context->gh.cctk_timefac = coordinates.time_refinement_factor;
  }
}

} // namespace CarpetX
#endif
