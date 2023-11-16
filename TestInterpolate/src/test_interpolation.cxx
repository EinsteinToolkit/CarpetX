#include <cctk.h>
#include <cctk_Arguments.h>
#include <util_Table.h>

namespace TestInterpolate {

#define DIM(v) ((int)(sizeof(v)/sizeof((v)[0])))

extern "C" void TestInterpolate_test_interpolation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestInterpolate_test_interpolation;

  constexpr CCTK_INT N_dims = 3;

  const CCTK_INT all_operations[1 + N_dims + N_dims*(N_dims-1)] = {0, 1, 2, 3, 11, 12, 13, 22, 23, 33};

  const CCTK_INT all_varinds[N_dims] = {
      CCTK_VarIndex("CoordinatesX::vcoordx"),
      CCTK_VarIndex("CoordinatesX::vcoordy"),
      CCTK_VarIndex("CoordinatesX::vcoordz"),
  };

  constexpr int nvars = DIM(all_varinds) * DIM(all_operations);

  CCTK_INT varinds[DIM(all_operations)][DIM(all_varinds)];
  CCTK_INT operations[DIM(all_operations)][DIM(all_varinds)];
  for (int op = 0 ; op < DIM(all_operations) ; op++) {
    for (int var = 0 ; var < DIM(all_varinds) ; var++) {
      varinds[op][var] = all_varinds[var];
      operations[op][var] = all_operations[op];
    }
  }

  constexpr int npoints = 4;
  const CCTK_REAL coords[N_dims][npoints] = {
      {0.1, 0.2, 0.3, 0.1},
      {0.1, 0.2, 0.3, 0.2},
      {0.1, 0.2, 0.3, 0.3},
  };

  CCTK_REAL *resultptrs[DIM(all_operations)][DIM(all_varinds)] = {
    interpolate_x_o1_op0,
    interpolate_y_o1_op0,
    interpolate_z_o1_op0,

    interpolate_x_o1_op1, interpolate_y_o1_op1, interpolate_z_o1_op1,
    interpolate_x_o1_op2, interpolate_y_o1_op2, interpolate_z_o1_op2,
    interpolate_x_o1_op3, interpolate_y_o1_op3, interpolate_z_o1_op3,

    interpolate_x_o1_op11, interpolate_y_o1_op11, interpolate_z_o1_op11,
    interpolate_x_o1_op12, interpolate_y_o1_op12, interpolate_z_o1_op12,
    interpolate_x_o1_op13, interpolate_y_o1_op13, interpolate_z_o1_op13,
    interpolate_x_o1_op22, interpolate_y_o1_op22, interpolate_z_o1_op22,
    interpolate_x_o1_op23, interpolate_y_o1_op23, interpolate_z_o1_op23,
    interpolate_x_o1_op33, interpolate_y_o1_op33, interpolate_z_o1_op33,
  };
  CCTK_REAL *resultptrs2[DIM(all_operations)][DIM(all_varinds)] = {
    driver_interpolate_x_o1_op0,
    driver_interpolate_y_o1_op0,
    driver_interpolate_z_o1_op0,

    driver_interpolate_x_o1_op1, driver_interpolate_y_o1_op1, driver_interpolate_z_o1_op1,
    driver_interpolate_x_o1_op2, driver_interpolate_y_o1_op2, driver_interpolate_z_o1_op2,
    driver_interpolate_x_o1_op3, driver_interpolate_y_o1_op3, driver_interpolate_z_o1_op3,

    driver_interpolate_x_o1_op11, driver_interpolate_y_o1_op11, driver_interpolate_z_o1_op11,
    driver_interpolate_x_o1_op12, driver_interpolate_y_o1_op12, driver_interpolate_z_o1_op12,
    driver_interpolate_x_o1_op13, driver_interpolate_y_o1_op13, driver_interpolate_z_o1_op13,
    driver_interpolate_x_o1_op22, driver_interpolate_y_o1_op22, driver_interpolate_z_o1_op22,
    driver_interpolate_x_o1_op23, driver_interpolate_y_o1_op23, driver_interpolate_z_o1_op23,
    driver_interpolate_x_o1_op33, driver_interpolate_y_o1_op33, driver_interpolate_z_o1_op33,
  };
  CCTK_REAL *referenceptrs[DIM(all_operations)][DIM(all_varinds)] = {
    reference_x_o1_op0,
    reference_y_o1_op0,
    reference_z_o1_op0,

    reference_x_o1_op1, reference_y_o1_op1, reference_z_o1_op1,
    reference_x_o1_op2, reference_y_o1_op2, reference_z_o1_op2,
    reference_x_o1_op3, reference_y_o1_op3, reference_z_o1_op3,

    reference_x_o1_op11, reference_y_o1_op11, reference_z_o1_op11,
    reference_x_o1_op12, reference_y_o1_op12, reference_z_o1_op12,
    reference_x_o1_op13, reference_y_o1_op13, reference_z_o1_op13,
    reference_x_o1_op22, reference_y_o1_op22, reference_z_o1_op22,
    reference_x_o1_op23, reference_y_o1_op23, reference_z_o1_op23,
    reference_x_o1_op33, reference_y_o1_op33, reference_z_o1_op33,
  };

  Interpolate(cctkGH, npoints, coords[0], coords[1], coords[2], nvars,
      (CCTK_INT const * const)varinds, (CCTK_INT const * const)operations, 1,
      (CCTK_REAL **)resultptrs);

  const void* interp_coords[N_dims] = {
    coords[0], coords[1], coords[2]
  };
  void *const *output_array = (void *const *)resultptrs2;

  /* DriverInterpolate arguments that aren't currently used */
  CCTK_INT const coord_system_handle = 0;
  CCTK_INT const interp_coords_type_code = 0;
  CCTK_INT const output_array_type_codes[1] = {0};
  int interp_handle = CCTK_InterpHandle("CarpetX");

  if(interp_handle < 0)
    CCTK_VERROR("Could not obtain inteprolator handle for built-in 'CarpetX' interpolator: %d", interp_handle);

  /* Table generation */
  int operands[DIM(all_operations)][DIM(all_varinds)];
  for (int op = 0 ; op < DIM(all_operations) ; op++) {
    for (int var = 0 ; var < DIM(all_varinds) ; var++) {
      operands[op][var] = var;
    }
  }

  int ierr;
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (param_table_handle < 0)
    CCTK_VERROR("Can't create parameter table: %d", param_table_handle);
  if ((ierr = Util_TableSetInt(param_table_handle, 1, "order")) < 0)
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, nvars, (int const*const)operands,
                            "operand_indices")) < 0)
    CCTK_VERROR("Can't set operand_indices array in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, nvars, (int const*const)operations,
                            "operation_codes")) < 0)
    CCTK_VERROR("Can't set operation_codes array in parameter table: %d", ierr);

  DriverInterpolate(
      cctkGH, N_dims, interp_handle, param_table_handle, coord_system_handle,
      npoints, interp_coords_type_code, interp_coords, nvars, (int const*const)varinds,
      nvars, output_array_type_codes, output_array);

  for (int op = 0 ; op < DIM(all_operations) ; op++) {
    for (int var = 0 ; var < DIM(all_varinds) ; var++) {
      for (int n = 0; n < npoints; ++n) {
        const int operand = operations[op][var];
        switch (operand) {
        case 0:
          referenceptrs[op][var][n] = coords[var][n];
          break;
        case 1:
        case 2:
        case 3:
          referenceptrs[op][var][n] = (operand - 1 == var ? 1. : 0.);
          break;
        case 11:
        case 12:
        case 13:
        case 22:
        case 23:
        case 33:
          referenceptrs[op][var][n] = 0.;
          break;
        default:
          CCTK_VERROR("Unsupported derivative operator %d", operand);
          break;
        }
      }
    }
  }
}

} // namespace TestInterpolate
