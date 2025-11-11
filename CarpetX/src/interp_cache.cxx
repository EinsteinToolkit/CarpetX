#include "interp_cache.hxx"
#include "driver.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cstring>

namespace CarpetX {

namespace {

// Simple hash combiner
inline std::uint64_t hash_combine(std::uint64_t h1, std::uint64_t h2) {
  return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
}

// Hash an array of integers
std::uint64_t hash_int_array(const int *data, std::size_t size) {
  std::uint64_t hash = 0x9e3779b97f4a7c15ULL;
  for (std::size_t i = 0; i < size; ++i) {
    hash = hash_combine(hash, static_cast<std::uint64_t>(data[i]));
  }
  return hash;
}

} // anonymous namespace

std::uint64_t InterpTargetCache::compute_grid_state_hash() const {
  if (!ghext) {
    return 0;
  }

  std::uint64_t hash = 0x9e3779b97f4a7c15ULL;

  // Hash number of patches
  hash = hash_combine(hash, static_cast<std::uint64_t>(ghext->num_patches()));

  // Hash each patch's level structure
  for (const auto &patchdata : ghext->patchdata) {
    const int patch = patchdata.patch;
    hash = hash_combine(hash, static_cast<std::uint64_t>(patch));

    // Hash number of levels
    const int nlevels = patchdata.amrcore->finestLevel() + 1;
    hash = hash_combine(hash, static_cast<std::uint64_t>(nlevels));

    // Hash each level's BoxArray and DistributionMapping
    for (int level = 0; level < nlevels; ++level) {
      const amrex::BoxArray &ba = patchdata.amrcore->boxArray(level);
      const amrex::DistributionMapping &dm =
          patchdata.amrcore->DistributionMap(level);

      // Hash BoxArray size
      hash = hash_combine(hash, static_cast<std::uint64_t>(ba.size()));

      // Hash each box
      for (int i = 0; i < ba.size(); ++i) {
        const amrex::Box &box = ba[i];
        const amrex::IntVect &lo = box.smallEnd();
        const amrex::IntVect &hi = box.bigEnd();

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          hash = hash_combine(hash, static_cast<std::uint64_t>(lo[d]));
          hash = hash_combine(hash, static_cast<std::uint64_t>(hi[d]));
        }
      }

      // Hash DistributionMapping
      hash = hash_combine(hash, static_cast<std::uint64_t>(dm.size()));
      const amrex::Vector<int> &proc_map = dm.ProcessorMap();
      for (int i = 0; i < proc_map.size(); ++i) {
        hash = hash_combine(hash, static_cast<std::uint64_t>(proc_map[i]));
      }
    }
  }

  return hash;
}

bool InterpTargetCache::coords_match(CCTK_INT npoints,
                                     const CCTK_REAL *globalsx,
                                     const CCTK_REAL *globalsy,
                                     const CCTK_REAL *globalsz) const {

  if (!cache_entry_.valid) {
    return false;
  }

  if (cache_entry_.npoints != npoints) {
    return false;
  }

  // Compare coordinates
  if (std::memcmp(cache_entry_.globalsx.data(), globalsx,
                  npoints * sizeof(CCTK_REAL)) != 0) {
    return false;
  }

  if (std::memcmp(cache_entry_.globalsy.data(), globalsy,
                  npoints * sizeof(CCTK_REAL)) != 0) {
    return false;
  }

  if (std::memcmp(cache_entry_.globalsz.data(), globalsz,
                  npoints * sizeof(CCTK_REAL)) != 0) {
    return false;
  }

  return true;
}

const InterpTargetCache::CacheEntry *
InterpTargetCache::get_cached_entry(CCTK_INT npoints, const CCTK_REAL *globalsx,
                                    const CCTK_REAL *globalsy,
                                    const CCTK_REAL *globalsz) {

  std::lock_guard<std::mutex> lock(mutex_);

  // Check if cache is valid
  if (!cache_entry_.valid) {
    ++stats_.misses;
    return nullptr;
  }

  // Compute current grid state hash
  const std::uint64_t current_hash = compute_grid_state_hash();

  // Check if grid has changed
  if (cache_entry_.grid_state_hash != current_hash) {
    ++stats_.invalidations;
    cache_entry_.valid = false;
    cache_entry_.containers.clear();
    return nullptr;
  }

  // Check if coordinates match
  if (!coords_match(npoints, globalsx, globalsy, globalsz)) {
    ++stats_.misses;
    return nullptr;
  }

  // Cache hit!
  ++stats_.hits;
  return &cache_entry_;
}

void InterpTargetCache::store_entry(
    CCTK_INT npoints, const CCTK_REAL *globalsx, const CCTK_REAL *globalsy,
    const CCTK_REAL *globalsz,
    std::vector<std::shared_ptr<Container> > &&containers) {

  std::lock_guard<std::mutex> lock(mutex_);

  // Store grid state hash
  cache_entry_.grid_state_hash = compute_grid_state_hash();

  // Store containers
  cache_entry_.containers = std::move(containers);

  // Store npoints
  cache_entry_.npoints = npoints;

  // Store coordinates for validation
  cache_entry_.globalsx.assign(globalsx, globalsx + npoints);
  cache_entry_.globalsy.assign(globalsy, globalsy + npoints);
  cache_entry_.globalsz.assign(globalsz, globalsz + npoints);

  // Mark as valid
  cache_entry_.valid = true;

  // Update statistics
  ++stats_.rebuilds;
}

void InterpTargetCache::invalidate() {
  std::lock_guard<std::mutex> lock(mutex_);

  if (cache_entry_.valid) {
    ++stats_.invalidations;
    cache_entry_.valid = false;
    cache_entry_.containers.clear();
  }
}

void InterpTargetCache::print_stats() const {
  DECLARE_CCTK_PARAMETERS;

  if (!verbose) {
    return;
  }

  std::lock_guard<std::mutex> lock(mutex_);

  const std::size_t total_accesses = stats_.hits + stats_.misses;
  if (total_accesses == 0) {
    return;
  }

  const double hit_rate = 100.0 * stats_.hits / total_accesses;

  CCTK_VINFO("InterpTargetCache statistics:");
  CCTK_VINFO("  Total accesses: %zu", total_accesses);
  CCTK_VINFO("  Hits:           %zu (%.2f)", stats_.hits, hit_rate);
  CCTK_VINFO("  Misses:         %zu", stats_.misses);
  CCTK_VINFO("  Rebuilds:       %zu", stats_.rebuilds);
  CCTK_VINFO("  Invalidations:  %zu", stats_.invalidations);
}

void InterpTargetCache::reset_stats() {
  std::lock_guard<std::mutex> lock(mutex_);
  stats_ = CacheStats{};
}

} // namespace CarpetX
