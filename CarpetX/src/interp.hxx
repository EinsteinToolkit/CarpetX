#ifndef CARPETX_CARPETX_INTERP_HXX
#define CARPETX_CARPETX_INTERP_HXX

#include <cctk.h>

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
