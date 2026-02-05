#ifndef CARPETX_CARPETX_INTERP_HXX
#define CARPETX_CARPETX_INTERP_HXX

#include <cctk.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <vector>

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

struct InterpSetup {
  using Container = amrex::AmrParticleContainer<3, 2>;

  int npoints;
  bool allow_boundaries;

  // Redistributed particle containers (one per patch)
  std::vector<Container> containers;

  // Symmetry state from Phase B (empty if no z-reflection)
  std::vector<bool> symmetry_reflected_z;
};

// Internal setup: performs coordinate transformation, symmetry handling,
// particle creation, and redistribution (Phases A-E).
InterpSetup interp_setup(const cGH *cctkGH, int npoints,
                         const CCTK_REAL *globalsx,
                         const CCTK_REAL *globalsy,
                         const CCTK_REAL *globalsz,
                         bool allow_boundaries);

// Internal apply: performs per-variable interpolation, MPI result gathering,
// and symmetry application (Phases F-H) using a pre-built setup.
// Takes non-const InterpSetup& because AMReX's ParIter mutates internal state.
void interp_apply(InterpSetup &setup, const cGH *cctkGH,
                  int nvars, const CCTK_INT *varinds,
                  const CCTK_INT *operations,
                  CCTK_POINTER resultptrs);

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
