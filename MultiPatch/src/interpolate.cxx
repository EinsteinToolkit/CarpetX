#include "multipatch.hxx"

// TODO: Don't include files from other thorn; create a proper interface
#include "../../CarpetX/src/driver.hxx"
#include "../../CarpetX/src/schedule.hxx"

#include <loop.hxx>

#include <cctk.h>

#include <array>
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace {
// <https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x>
constexpr inline std::size_t hash_combine(std::size_t h1, std::size_t h2) {
  return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}
} // namespace

namespace MultiPatch {

struct Location {
  int patch;
  int level;
  int index;
  int block;
};
} // namespace MultiPatch

namespace std {
template <> struct equal_to<MultiPatch::Location> {
  bool operator()(const MultiPatch::Location &x,
                  const MultiPatch::Location &y) const {
    return std::equal_to<std::array<int, 4> >()(
        std::array<int, 4>{x.patch, x.level, x.index, x.block},
        std::array<int, 4>{y.patch, y.level, y.index, y.block});
  }
};
template <> struct less<MultiPatch::Location> {
  bool operator()(const MultiPatch::Location &x,
                  const MultiPatch::Location &y) const {
    return std::less<std::array<int, 4> >()(
        std::array<int, 4>{x.patch, x.level, x.index, x.block},
        std::array<int, 4>{y.patch, y.level, y.index, y.block});
  }
};
template <> struct hash<MultiPatch::Location> {
  std::size_t operator()(const MultiPatch::Location &x) const {
    return hash_combine(hash_combine(hash_combine(std::hash<int>()(x.patch),
                                                  std::hash<int>()(x.level)),
                                     std::hash<int>()(x.index)),
                        std::hash<int>()(x.block));
  }
};
} // namespace std

namespace MultiPatch {
using SourcePoints = std::array<std::vector<CCTK_REAL>, dim>;

extern "C" void
MultiPatch1_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                        const CCTK_INT nvars_,
                        const CCTK_INT *restrict const varinds_) {
  assert(cctkGH_);
  assert(nvars_ >= 0);
  assert(varinds_);
  const cGH *const cctkGH = static_cast<const cGH *>(cctkGH_);
  const std::vector<int> varinds(varinds_, varinds_ + nvars_);

  // Step 0: Check input

  for (const int varind : varinds) {
    assert(varind >= 0);
    const int gi = CCTK_GroupIndexFromVarI(varind);
    assert(gi >= 0);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = varind - v0;
    assert(vi >= 0);
    cGroup group_data;
    const int ierr = CCTK_GroupData(gi, &group_data);
    assert(!ierr);
    assert(group_data.grouptype == CCTK_GF);
    assert(group_data.vartype == CCTK_VARIABLE_REAL);
    assert(group_data.disttype == CCTK_DISTRIB_DEFAULT);
    assert(group_data.dim == dim);
    // TODO: Check centering
  }

  // Step 1: Find coordinates where we need to interpolate

  std::map<Location, SourcePoints> source_mapping;
  loop_over_blocks(active_levels_t(), [&](int patch, int level, int index,
                                          int block, const cGH *cctkGH) {
    const Loop::GridDescBase grid(cctkGH);
    const std::array<int, dim> centering{0, 0, 0};
    const Loop::GF3D2layout layout(cctkGH, centering);

    const Location location{patch, level, index, block};

    const std::array<Loop::GF3D2<const CCTK_REAL>, dim> vcoords{
        Loop::GF3D2<const CCTK_REAL>(
            layout, static_cast<const CCTK_REAL *>(
                        CCTK_VarDataPtr(cctkGH, 0, "Coordinates::vcoordx"))),
        Loop::GF3D2<const CCTK_REAL>(
            layout, static_cast<const CCTK_REAL *>(
                        CCTK_VarDataPtr(cctkGH, 0, "Coordinates::vcoordy"))),
        Loop::GF3D2<const CCTK_REAL>(
            layout, static_cast<const CCTK_REAL *>(
                        CCTK_VarDataPtr(cctkGH, 0, "Coordinates::vcoordz")))};

    SourcePoints source_points;
    // Note: This includes symmetry points
    grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const auto &normal{p.NI};
      const auto &current_patch{the_patch_system->patches.at(patch)};
      const auto &patch_faces{current_patch.faces};

      // TODO: Do not continue, do a loop, clen up.
      // -x face
      if (normal[0] < 0 && patch_faces[0][0].is_outer_boundary) {
        return;
      }

      // +x face
      if (normal[0] > 0 && patch_faces[1][0].is_outer_boundary) {
        return;
      }

      // -y face
      if (normal[1] < 0 && patch_faces[0][1].is_outer_boundary) {
        return;
      }

      // +y face
      if (normal[1] > 0 && patch_faces[1][1].is_outer_boundary) {
        return;
      }

      // -z face
      if (normal[2] < 0 && patch_faces[0][2].is_outer_boundary) {
        return;
      }

      // +z face
      if (normal[2] > 0 && patch_faces[1][2].is_outer_boundary) {
        return;
      }

      for (int d = 0; d < dim; ++d)
        source_points[d].push_back(vcoords[d](p.I));
    });
#pragma omp critical
    source_mapping[location] = std::move(source_points);
  });

  // Step 2: Interpolate to these coordinates

  // Gather all coordinates
  std::array<std::vector<CCTK_REAL>, dim> coords;
  for (const auto &[location, source_points] : source_mapping)
    for (int d = 0; d < dim; ++d)
      coords[d].insert(coords[d].end(), source_points[d].begin(),
                       source_points[d].end());

  const std::size_t nvars = varinds.size();
  const std::size_t npoints = coords[0].size();

  std::vector<CCTK_INT> operations(nvars, 0);

  // Allocate memory for values
  std::vector<std::vector<CCTK_REAL> > results(nvars);
  std::vector<CCTK_REAL *> resultptrs(nvars);
  for (size_t n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
  }

  // Interpolate
  Interpolate(cctkGH, coords[0].size(), coords[0].data(), coords[1].data(),
              coords[2].data(), varinds.size(), varinds.data(),
              operations.data(), resultptrs.data());

  // Scatter interpolated values
  std::map<Location, std::vector<std::vector<CCTK_REAL> > > result_mapping;
  std::size_t pos = 0;
  for (const auto &[location, source_points] : source_mapping) {
    const std::size_t length = source_points[0].size();
    std::vector<std::vector<CCTK_REAL> > result_values(nvars);
    for (std::size_t n = 0; n < nvars; ++n)
      result_values.at(nvars).insert(result_values.at(nvars).begin(),
                                     &results.at(n).at(pos),
                                     &results.at(n).at(pos + length));
    result_mapping[location] = std::move(result_values);
  }

  // Step 3: Write back results

  loop_over_blocks(active_levels_t(), [&](int patch, int level, int index,
                                          int block, const cGH *cctkGH) {
    const Loop::GridDescBase grid(cctkGH);
    const std::array<int, dim> centering{0, 0, 0};
    const Loop::GF3D2layout layout(cctkGH, centering);

    const Location location{patch, level, index, block};
    const std::vector<std::vector<CCTK_REAL> > &result_values =
        result_mapping.at(location);

    for (std::size_t n = 0; n < nvars; ++n) {
      const std::vector<CCTK_REAL> &result_values_n = result_values.at(n);

      const Loop::GF3D2<CCTK_REAL> var(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, varinds.at(n))));
      std::size_t pos = 0;
      // Note: This includes symmetry points
      grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
        const auto &normal{p.NI};
        const auto &current_patch{the_patch_system->patches.at(patch)};
        const auto &patch_faces{current_patch.faces};

        // TODO: See line 123
        // -x face
        if (normal[0] < 0 && patch_faces[0][0].is_outer_boundary) {
          return;
        }

        // +x face
        if (normal[0] > 0 && patch_faces[1][0].is_outer_boundary) {
          return;
        }

        // -y face
        if (normal[1] < 0 && patch_faces[0][1].is_outer_boundary) {
          return;
        }

        // +y face
        if (normal[1] > 0 && patch_faces[1][1].is_outer_boundary) {
          return;
        }

        // -z face
        if (normal[2] < 0 && patch_faces[0][2].is_outer_boundary) {
          return;
        }

        // +z face
        if (normal[2] > 0 && patch_faces[1][2].is_outer_boundary) {
          return;
        }

        var(p.I) = result_values_n[pos++];
      });
      assert(pos == result_values.at(0).size());
    }
  });
}

} // namespace MultiPatch
