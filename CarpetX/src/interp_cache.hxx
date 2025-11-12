#ifndef CARPETX_CARPETX_INTERP_CACHE_HXX
#define CARPETX_CARPETX_INTERP_CACHE_HXX

#include <cctk.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <memory>
#include <vector>
#include <mutex>
#include <cstdint>

namespace CarpetX {

/**
 * Cache for interpolation particle containers to avoid recreation and
 * redistribution when the grid doesn't change via AMR.
 *
 * This is a Meyers singleton that caches:
 * - AMReX particle containers with interpolation target coordinates
 * - The grid state hash to detect when AMR has changed the grid
 *
 * When verbose output is enabled, it tracks and reports cache statistics.
 */
class InterpTargetCache {
public:
  using Container = amrex::AmrParticleContainer<3, 2>;
  using Particle = Container::ParticleType;

  struct CacheEntry {
    // Hash of the grid state (to detect AMR changes)
    std::uint64_t grid_state_hash{0};

    // Number of interpolation points
    CCTK_INT npoints{0};

    // Cached coordinates used for validation (global coordinates)
    std::vector<CCTK_REAL> globalsx;
    std::vector<CCTK_REAL> globalsy;
    std::vector<CCTK_REAL> globalsz;

    // Cached AMReX particle containers (one per patch) - THE EXPENSIVE PART
    // These are stored as shared_ptr to allow transfer without deep copy
    std::vector<std::shared_ptr<Container> > containers;

    // Whether cache entry is valid
    bool valid{false};
  };

  struct CacheStats {
    std::size_t hits{0};
    std::size_t misses{0};
    std::size_t rebuilds{0};
    std::size_t invalidations{0};
  };

  // Meyers singleton accessor
  static InterpTargetCache &instance() {
    static InterpTargetCache cache;
    return cache;
  }

  // Delete copy and move constructors/assignment
  InterpTargetCache(const InterpTargetCache &) = delete;
  InterpTargetCache &operator=(const InterpTargetCache &) = delete;
  InterpTargetCache(InterpTargetCache &&) = delete;
  InterpTargetCache &operator=(InterpTargetCache &&) = delete;

  /**
   * Compute a hash of the current grid state to detect AMR changes.
   * This includes:
   * - Number of patches
   * - Number of levels per patch
   * - BoxArray and DistributionMapping for each level
   */
  std::uint64_t compute_grid_state_hash() const;

  /**
   * Try to get cached particle containers for given interpolation coordinates.
   * Returns nullptr if cache miss or grid has changed.
   *
   * The containers can be used directly for interpolation, avoiding the
   * expensive particle creation and redistribution steps.
   */
  const CacheEntry *get_cached_entry(CCTK_INT npoints,
                                     const CCTK_REAL *globalsx,
                                     const CCTK_REAL *globalsy,
                                     const CCTK_REAL *globalsz);

  /**
   * Store particle containers in cache for future reuse.
   *
   * Takes shared_ptr to allow efficient sharing without copying the containers.
   */
  void store_entry(CCTK_INT npoints, const CCTK_REAL *globalsx,
                   const CCTK_REAL *globalsy, const CCTK_REAL *globalsz,
                   std::vector<std::shared_ptr<Container> > &&containers);

  /**
   * Invalidate the cache (called when grid changes via AMR).
   */
  void invalidate();

  /**
   * Get cache statistics.
   */
  const CacheStats &get_stats() const { return stats_; }

  /**
   * Print cache statistics (called when verbose=yes).
   */
  void print_stats() const;

  /**
   * Reset cache statistics.
   */
  void reset_stats();

private:
  InterpTargetCache() = default;
  ~InterpTargetCache() = default;

  // Check if coordinates match cached coordinates
  bool coords_match(CCTK_INT npoints, const CCTK_REAL *coordsx,
                    const CCTK_REAL *coordsy, const CCTK_REAL *coordsz) const;

  // Mutex for thread-safe access
  mutable std::mutex mutex_;

  // Cached entry
  CacheEntry cache_entry_;

  // Statistics
  CacheStats stats_;
};

} // namespace CarpetX

#endif // CARPETX_CARPETX_INTERP_CACHE_HXX
