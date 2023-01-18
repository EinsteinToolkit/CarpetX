#ifndef MULTIPATCH_HXX
#define MULTIPATCH_HXX

#include <fixmath.hxx>
#include <cctk.h>

#include <loop.hxx>
#include <mat.hxx>
#include <tuple.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace MultiPatch {
using namespace Arith;

constexpr int dim = 3;

struct PatchFace {
  bool is_outer_boundary;
  int other_patch; // -1 if outer boundary
};

struct Patch {
  std::string name;
  vect<int, dim> ncells;
  vect<CCTK_REAL, dim> xmin, xmax; // cell boundaries
  bool is_cartesian;               // Jacobian is trivial
  vect<vect<PatchFace, dim>, 2> faces;
};

struct PatchTransformations {
  // Cartesian
  CCTK_REAL cartesian_xmax;
  CCTK_REAL cartesian_xmin;
  CCTK_REAL cartesian_ymax;
  CCTK_REAL cartesian_ymin;
  CCTK_REAL cartesian_zmax;
  CCTK_REAL cartesian_zmin;
  int cartesian_ncells_i;
  int cartesian_ncells_j;
  int cartesian_ncells_k;

  // Cubed sphere
  CCTK_REAL cubed_sphere_rmin;
  CCTK_REAL cubed_sphere_rmax;

  // Swirl
  int swirl_ncells_i;
  int swirl_ncells_j;
  int swirl_ncells_k;

  // Cake
  // The radius of the outer boundary
  CCTK_REAL cake_outer_boundary_radius;
  // Half the coordinate length of the central cartesian cube's face
  CCTK_REAL cake_inner_boundary_radius;
  // The number of cells in the x direction of the central cartesian cube
  int cake_cartesian_ncells_i;
  // The number of cells in the y direction of the central cartesian cube
  int cake_cartesian_ncells_j;
  // The number of cells in the z direction of the central cartesian cube
  int cake_cartesian_ncells_k;
  // The number of cells in the angular direction of spherical patches
  int cake_angular_cells;
  // The number of cells in the radial direction of spherical patches
  int cake_radial_cells;

  PatchTransformations();

  PatchTransformations(const PatchTransformations &) = default;
  PatchTransformations(PatchTransformations &&) = default;
  PatchTransformations &operator=(const PatchTransformations &) = default;
  PatchTransformations &operator=(PatchTransformations &&) = default;

  std_tuple<int, vec<CCTK_REAL, dim> > (*global2local)(
      const PatchTransformations &pt, const vec<CCTK_REAL, dim> &x) = 0;

  vec<CCTK_REAL, dim> (*local2global)(const PatchTransformations &pt, int patch,
                                      const vec<CCTK_REAL, dim> &a) = 0;

  // Calculating global derivatives d/dx from local derivatives d/da
  // requires the Jacobian da/dx, and also its derivative d^2/dx^2 for
  // second derivatives

  // da/dx[i,j] = da[i] / dx[j]
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> > (
      *dlocal_dglobal)(const PatchTransformations &pt, int patch,
                       const vec<CCTK_REAL, dim> &a) = 0;

  // d^2a/dx^2[i,j,k] = d^2a[i] / dx^2[j,k]
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
            vec<smat<CCTK_REAL, dim>, dim> > (*d2local_dglobal2)(
      const PatchTransformations &pt, int patch,
      const vec<CCTK_REAL, dim> &a) = 0;

  // Device functions mirroring the function above
  std_tuple<int, vec<CCTK_REAL, dim> > (*global2local_device)(
      const PatchTransformations &pt, const vec<CCTK_REAL, dim> &x) = 0;
  vec<CCTK_REAL, dim> (*local2global_device)(const PatchTransformations &pt,
                                             int patch,
                                             const vec<CCTK_REAL, dim> &a) = 0;
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> > (
      *dlocal_dglobal_device)(const PatchTransformations &pt, int patch,
                              const vec<CCTK_REAL, dim> &a) = 0;
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
            vec<smat<CCTK_REAL, dim>, dim> > (*d2local_dglobal2_device)(
      const PatchTransformations &pt, int patch,
      const vec<CCTK_REAL, dim> &a) = 0;
};

struct PatchSystem {
  std::string name;
  std::vector<Patch> patches;
  int num_patches() const { return patches.size(); }

  PatchTransformations transformations;

  PatchSystem() {}
  PatchSystem(std::string name, std::vector<Patch> patches,
              PatchTransformations transformations)
      : name(std::move(name)), patches(std::move(patches)),
        transformations(std::move(transformations)) {}
};

////////////////////////////////////////////////////////////////////////////////

PatchSystem SetupCartesian();
PatchSystem SetupCubedSphere();
PatchSystem SetupSwirl();
PatchSystem SetupCake();

extern std::unique_ptr<PatchSystem> the_patch_system;

} // namespace MultiPatch

#endif // #ifndef MULTIPATCH_HXX
