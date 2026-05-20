#ifndef CARPETX_CARPETX_INTERP_HXX
#define CARPETX_CARPETX_INTERP_HXX

#include <cctk.h>

#include <defs.hxx>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <vector>
#include <vect.hxx>

// Scheduled functions
extern "C" int CarpetX_InterpGridArrays(
    cGH const *const cGH, int const N_dims, int const local_interp_handle,
    int const param_table_handle, int const coord_system_handle,
    int const N_interp_points, int const interp_coords_type_code,
    void const *const interp_coords[], int const N_input_arrays,
    CCTK_INT const input_array_variable_indices[], int const N_output_arrays,
    CCTK_INT const output_array_type_codes[], void *const output_arrays[]);

extern "C" CCTK_INT CarpetX_DriverInterpolate(
    CCTK_POINTER_TO_CONST const cctkGH, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type_code,
    CCTK_POINTER_TO_CONST const coords_list[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_variable_indices[],
    CCTK_INT const N_output_arrays, CCTK_INT const output_array_type_codes[],
    CCTK_POINTER const output_arrays[]);

namespace CarpetX {

// Create a cache of data required for interpolation
class InterpolationSetup {
private:
  using Container = amrex::AmrParticleContainer<3, 2>;
  using Particle = Container::ParticleType;
  using PinnedParticleTile = amrex::ParticleContainer_impl<
      Particle, 0, 0, amrex::PinnedArenaAllocator>::ParticleTileType;
  using ParticleTile = Container::ParticleTileType;

  const int npoints{0};
  std::vector<bool> symmetry_reflected_z;
  std::vector<Container> containers{}; // [patch]

public:

  /*
   * Constructs an InterpolationSetup by distributing interpolation points
   * across MPI processes via AMReX particle containers.
   *
   * Given npoints global coordinates (globalsx, globalsy, globalsz), this
   * constructor performs all the setup work required before the actual field
   * interpolation:
   *
   * 1. Coordinate conversion: global coordinates are converted to patch-local
   *    coordinates via MultiPatch_GlobalToLocal2 when available, or assigned to
   *    patch 0 otherwise.
   *
   * 2. Symmetry handling: if the reflection_z parameter is active, points with
   *    z < domain_zmin are reflected about the lower-z boundary. Which points
   *    were reflected is recorded in symmetry_reflected_z for later use by
   *    Interpolate().
   *
   * 3. Domain projection: particle positions are clamped to lie at least half a
   *    coarse-grid cell spacing inside the domain so that AMReX does not
   *    silently discard out-of-domain particles during redistribution. The
   *    actual interpolation still uses the original (pre-clamped) local
   *    coordinates stored in the particle's real-data slots.
   *
   * 4. Particle distribution: particles are inserted into per-patch AMReX
   *    particle containers and redistributed with Redistribute() so that each
   *    MPI rank owns exactly the points that overlap with its AMReX boxes.
   *
   * @param cctkGH   Cactus grid hierarchy handle (must be in global mode).
   * @param npoints  Number of interpolation points.
   * @param globalsx x-coordinates of the interpolation points in the global
   *                 frame.
   * @param globalsy y-coordinates of the interpolation points in the global
   *                 frame.
   * @param globalsz z-coordinates of the interpolation points in the global
   *                 frame.
   */
  InterpolationSetup(CCTK_ATTRIBUTE_UNUSED const cGH *restrict const cctkGH,
                     const CCTK_INT npoints,
                     const CCTK_REAL *restrict const globalsx,
                     const CCTK_REAL *restrict const globalsy,
                     const CCTK_REAL *restrict const globalsz);

  /*
   * Interpolate performs the actual grid interpolation given a
   * pre-built InterpolationSetup.
   *
   * allowed_boundaries[patch][f][d] controls, for each patch face
   * (face f=0/1, direction d), whether an interpolation stencil is permitted to
   * anchor in the boundary/ghost-zone region on that face. The convention
   * matches the one used by the interpolator struct (see its comment near
   * i0_allowed):
   *
   * true = stencil may anchor right up to that face (ghost zone data there is
   * assumed valid).
   *
   * false = stencil is pushed inward by nghostzones on that face (ghost zone
   * data is considered unavailable).
   *
   * The per-box `allowed_boundaries` passed to each interpolator is derived
   * from this policy via grid.bbox[f][d], which AMReX sets to true when the
   * box face touches the outer boundary of its patch's AMReX domain:
   *
   * When bbox[f][d] is true, the box face is at the patch outer boundary. In
   * this case we read allowed_boundaries for that face. If allowed_boundaries
   * is true, that means some routine (i.e., BCs) has filled the ghost zones.
   * When it is false ghost zones are filled by inter-patch interpolation and we
   * push the stencil anchor inward.
   *
   * When bbox[f][d] is false, the box face is interior to the patch, between
   * AMReX boxes. In this case AMReX fill-patch operations guarantee these ghost
   * zones are valid.
   */
  void Interpolate(CCTK_ATTRIBUTE_UNUSED const cGH *restrict const cctkGH,
                   const CCTK_INT nvars, const CCTK_INT *restrict const varinds,
                   const CCTK_INT *restrict const operations,
                   const std::vector<Arith::vect<Arith::vect<bool, 3>, 2> >
                       &allowed_boundaries, //  [patch][face][direction]
                   const CCTK_POINTER resultptrs_) const;
};

// a dummy routine for now
// TODO: implement this for actual local interpolation
int InterpLocalUniform(int N_dims, int param_table_handle,
                       /***** coordinate system *****/
                       const CCTK_REAL coord_origin[],
                       const CCTK_REAL coord_delta[],
                       /***** interpolation points *****/
                       int N_interp_points, int interp_coords_type_code,
                       const void *const interp_coords[],
                       /***** input arrays *****/
                       int N_input_arrays, const CCTK_INT input_array_dims[],
                       const CCTK_INT input_array_type_codes[],
                       const void *const input_arrays[],
                       /***** output arrays *****/
                       int N_output_arrays,
                       const CCTK_INT output_array_type_codes[],
                       void *const output_arrays[]);
} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_INTERP_HXX
