#ifndef MULTIPATCH_CAKE_JACOBIANS_HXX
#define MULTIPATCH_CAKE_JACOBIANS_HXX

namespace MultiPatch {
namespace Cake {

CCTK_DEVICE CCTK_HOST inline std_tuple<jac_t, djac_t>
cake_jacs(const PatchTransformations &pt, int patch,
          const svec &global_vars) {
  using std::pow;
  using std::sqrt;

  jac_t J{};
  djac_t dJ{};

  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  switch (patch) {

  case static_cast<int>(patch_piece::cartesian):
    J(0)(0) = 1;
    J(0)(1) = 0;
    J(0)(2) = 0;
    J(1)(0) = 0;
    J(1)(1) = 1;
    J(1)(2) = 0;
    J(2)(0) = 0;
    J(2)(1) = 0;
    J(2)(2) = 1;

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = 0;
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = 0;
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) = 0;
    dJ(2)(0, 1) = 0;
    dJ(2)(0, 2) = 0;
    dJ(2)(1, 0) = 0;
    dJ(2)(1, 1) = 0;
    dJ(2)(1, 2) = 0;
    dJ(2)(2, 0) = 0;
    dJ(2)(2, 1) = 0;
    dJ(2)(2, 2) = 0;
    break;

  case static_cast<int>(patch_piece::plus_x):
    J(0)(0) = -(z / pow(x, 2));
    J(0)(1) = 0;
    J(0)(2) = 1 / x;
    J(1)(0) = -(y / pow(x, 2));
    J(1)(1) = 1 / x;
    J(1)(2) = 0;
    J(2)
    (0) = d_global_to_local_cake_core_dvar(pt, z / x, y / x, x) -
          (y * d_global_to_local_cake_core_db(pt, z / x, y / x, x) +
           z * d_global_to_local_cake_core_da(pt, z / x, y / x, x)) /
              pow(x, 2);
    J(2)(1) = d_global_to_local_cake_core_db(pt, z / x, y / x, x) / x;
    J(2)(2) = d_global_to_local_cake_core_da(pt, z / x, y / x, x) / x;

    dJ(0)(0, 0) = (2 * z) / pow(x, 3);
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = -pow(x, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = -pow(x, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = (2 * y) / pow(x, 3);
    dJ(1)(0, 1) = -pow(x, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = -pow(x, -2);
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        (pow(x, 4) * d2_global_to_local_cake_core_dvar2(pt, z / x, y / x, x) +
         2 * x * y * d_global_to_local_cake_core_db(pt, z / x, y / x, x) -
         2 * pow(x, 2) * y *
             d2_global_to_local_cake_core_dbdvar(pt, z / x, y / x, x) +
         pow(y, 2) * d2_global_to_local_cake_core_db2(pt, z / x, y / x, x) +
         2 * x * z * d_global_to_local_cake_core_da(pt, z / x, y / x, x) -
         2 * pow(x, 2) * z *
             d2_global_to_local_cake_core_dadvar(pt, z / x, y / x, x) +
         2 * y * z * d2_global_to_local_cake_core_dadb(pt, z / x, y / x, x) +
         pow(z, 2) * d2_global_to_local_cake_core_da2(pt, z / x, y / x, x)) /
        pow(x, 4);
    dJ(2)(0, 1) =
        -((x * (d_global_to_local_cake_core_db(pt, z / x, y / x, x) -
                x * d2_global_to_local_cake_core_dbdvar(pt, z / x, y / x, x)) +
           y * d2_global_to_local_cake_core_db2(pt, z / x, y / x, x) +
           z * d2_global_to_local_cake_core_dadb(pt, z / x, y / x, x)) /
          pow(x, 3));
    dJ(2)(0, 2) =
        -((x * (d_global_to_local_cake_core_da(pt, z / x, y / x, x) -
                x * d2_global_to_local_cake_core_dadvar(pt, z / x, y / x, x)) +
           y * d2_global_to_local_cake_core_dadb(pt, z / x, y / x, x) +
           z * d2_global_to_local_cake_core_da2(pt, z / x, y / x, x)) /
          pow(x, 3));
    dJ(2)(1, 0) =
        -((x * (d_global_to_local_cake_core_db(pt, z / x, y / x, x) -
                x * d2_global_to_local_cake_core_dbdvar(pt, z / x, y / x, x)) +
           y * d2_global_to_local_cake_core_db2(pt, z / x, y / x, x) +
           z * d2_global_to_local_cake_core_dadb(pt, z / x, y / x, x)) /
          pow(x, 3));
    dJ(2)(1, 1) =
        d2_global_to_local_cake_core_db2(pt, z / x, y / x, x) / pow(x, 2);
    dJ(2)(1, 2) =
        d2_global_to_local_cake_core_dadb(pt, z / x, y / x, x) / pow(x, 2);
    dJ(2)(2, 0) =
        -((x * (d_global_to_local_cake_core_da(pt, z / x, y / x, x) -
                x * d2_global_to_local_cake_core_dadvar(pt, z / x, y / x, x)) +
           y * d2_global_to_local_cake_core_dadb(pt, z / x, y / x, x) +
           z * d2_global_to_local_cake_core_da2(pt, z / x, y / x, x)) /
          pow(x, 3));
    dJ(2)(2, 1) =
        d2_global_to_local_cake_core_dadb(pt, z / x, y / x, x) / pow(x, 2);
    dJ(2)(2, 2) =
        d2_global_to_local_cake_core_da2(pt, z / x, y / x, x) / pow(x, 2);
    break;

  case static_cast<int>(patch_piece::minus_x):
    J(0)(0) = z / pow(x, 2);
    J(0)(1) = 0;
    J(0)(2) = -(1 / x);
    J(1)(0) = -(y / pow(x, 2));
    J(1)(1) = 1 / x;
    J(1)(2) = 0;
    J(2)
    (0) = d_global_to_local_cake_core_dvar(pt, -(z / x), y / x, x) +
          (-(y * d_global_to_local_cake_core_db(pt, -(z / x), y / x, x)) +
           z * d_global_to_local_cake_core_da(pt, -(z / x), y / x, x)) /
              pow(x, 2);
    J(2)(1) = d_global_to_local_cake_core_db(pt, -(z / x), y / x, x) / x;
    J(2)(2) = -(d_global_to_local_cake_core_da(pt, -(z / x), y / x, x) / x);

    dJ(0)(0, 0) = (-2 * z) / pow(x, 3);
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = pow(x, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = pow(x, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = (2 * y) / pow(x, 3);
    dJ(1)(0, 1) = -pow(x, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = -pow(x, -2);
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        (pow(x, 4) *
             d2_global_to_local_cake_core_dvar2(pt, -(z / x), y / x, x) +
         2 * x * y * d_global_to_local_cake_core_db(pt, -(z / x), y / x, x) -
         2 * pow(x, 2) * y *
             d2_global_to_local_cake_core_dbdvar(pt, -(z / x), y / x, x) +
         pow(y, 2) * d2_global_to_local_cake_core_db2(pt, -(z / x), y / x, x) -
         2 * x * z * d_global_to_local_cake_core_da(pt, -(z / x), y / x, x) +
         2 * pow(x, 2) * z *
             d2_global_to_local_cake_core_dadvar(pt, -(z / x), y / x, x) -
         2 * y * z * d2_global_to_local_cake_core_dadb(pt, -(z / x), y / x, x) +
         pow(z, 2) * d2_global_to_local_cake_core_da2(pt, -(z / x), y / x, x)) /
        pow(x, 4);
    dJ(2)(0, 1) =
        (-(x * d_global_to_local_cake_core_db(pt, -(z / x), y / x, x)) +
         pow(x, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, -(z / x), y / x, x) -
         y * d2_global_to_local_cake_core_db2(pt, -(z / x), y / x, x) +
         z * d2_global_to_local_cake_core_dadb(pt, -(z / x), y / x, x)) /
        pow(x, 3);
    dJ(2)(0, 2) =
        (x * (d_global_to_local_cake_core_da(pt, -(z / x), y / x, x) -
              x * d2_global_to_local_cake_core_dadvar(pt, -(z / x), y / x, x)) +
         y * d2_global_to_local_cake_core_dadb(pt, -(z / x), y / x, x) -
         z * d2_global_to_local_cake_core_da2(pt, -(z / x), y / x, x)) /
        pow(x, 3);
    dJ(2)(1, 0) =
        (-(x * d_global_to_local_cake_core_db(pt, -(z / x), y / x, x)) +
         pow(x, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, -(z / x), y / x, x) -
         y * d2_global_to_local_cake_core_db2(pt, -(z / x), y / x, x) +
         z * d2_global_to_local_cake_core_dadb(pt, -(z / x), y / x, x)) /
        pow(x, 3);
    dJ(2)(1, 1) =
        d2_global_to_local_cake_core_db2(pt, -(z / x), y / x, x) / pow(x, 2);
    dJ(2)(1, 2) = -(d2_global_to_local_cake_core_dadb(pt, -(z / x), y / x, x) /
                    pow(x, 2));
    dJ(2)(2, 0) =
        (x * (d_global_to_local_cake_core_da(pt, -(z / x), y / x, x) -
              x * d2_global_to_local_cake_core_dadvar(pt, -(z / x), y / x, x)) +
         y * d2_global_to_local_cake_core_dadb(pt, -(z / x), y / x, x) -
         z * d2_global_to_local_cake_core_da2(pt, -(z / x), y / x, x)) /
        pow(x, 3);
    dJ(2)(2, 1) = -(d2_global_to_local_cake_core_dadb(pt, -(z / x), y / x, x) /
                    pow(x, 2));
    dJ(2)(2, 2) =
        d2_global_to_local_cake_core_da2(pt, -(z / x), y / x, x) / pow(x, 2);
    break;

  case static_cast<int>(patch_piece::plus_y):
    J(0)(0) = 0;
    J(0)(1) = -(z / pow(y, 2));
    J(0)(2) = 1 / y;
    J(1)(0) = -(1 / y);
    J(1)(1) = x / pow(y, 2);
    J(1)(2) = 0;
    J(2)(0) = -(d_global_to_local_cake_core_db(pt, z / y, -(x / y), y) / y);
    J(2)
    (1) = d_global_to_local_cake_core_dvar(pt, z / y, -(x / y), y) +
          (x * d_global_to_local_cake_core_db(pt, z / y, -(x / y), y) -
           z * d_global_to_local_cake_core_da(pt, z / y, -(x / y), y)) /
              pow(y, 2);
    J(2)(2) = d_global_to_local_cake_core_da(pt, z / y, -(x / y), y) / y;

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = (2 * z) / pow(y, 3);
    dJ(0)(1, 2) = -pow(y, -2);
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = -pow(y, -2);
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = pow(y, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = pow(y, -2);
    dJ(1)(1, 1) = (-2 * x) / pow(y, 3);
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        d2_global_to_local_cake_core_db2(pt, z / y, -(x / y), y) / pow(y, 2);
    dJ(2)(0, 1) =
        (y * (d_global_to_local_cake_core_db(pt, z / y, -(x / y), y) -
              y * d2_global_to_local_cake_core_dbdvar(pt, z / y, -(x / y), y)) -
         x * d2_global_to_local_cake_core_db2(pt, z / y, -(x / y), y) +
         z * d2_global_to_local_cake_core_dadb(pt, z / y, -(x / y), y)) /
        pow(y, 3);
    dJ(2)(0, 2) = -(d2_global_to_local_cake_core_dadb(pt, z / y, -(x / y), y) /
                    pow(y, 2));
    dJ(2)(1, 0) =
        (y * (d_global_to_local_cake_core_db(pt, z / y, -(x / y), y) -
              y * d2_global_to_local_cake_core_dbdvar(pt, z / y, -(x / y), y)) -
         x * d2_global_to_local_cake_core_db2(pt, z / y, -(x / y), y) +
         z * d2_global_to_local_cake_core_dadb(pt, z / y, -(x / y), y)) /
        pow(y, 3);
    dJ(2)(1, 1) =
        (pow(y, 4) *
             d2_global_to_local_cake_core_dvar2(pt, z / y, -(x / y), y) -
         2 * x * y * d_global_to_local_cake_core_db(pt, z / y, -(x / y), y) +
         2 * x * pow(y, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, z / y, -(x / y), y) +
         pow(x, 2) * d2_global_to_local_cake_core_db2(pt, z / y, -(x / y), y) +
         2 * y * z * d_global_to_local_cake_core_da(pt, z / y, -(x / y), y) -
         2 * pow(y, 2) * z *
             d2_global_to_local_cake_core_dadvar(pt, z / y, -(x / y), y) -
         2 * x * z * d2_global_to_local_cake_core_dadb(pt, z / y, -(x / y), y) +
         pow(z, 2) * d2_global_to_local_cake_core_da2(pt, z / y, -(x / y), y)) /
        pow(y, 4);
    dJ(2)(1, 2) =
        (-(y * d_global_to_local_cake_core_da(pt, z / y, -(x / y), y)) +
         pow(y, 2) *
             d2_global_to_local_cake_core_dadvar(pt, z / y, -(x / y), y) +
         x * d2_global_to_local_cake_core_dadb(pt, z / y, -(x / y), y) -
         z * d2_global_to_local_cake_core_da2(pt, z / y, -(x / y), y)) /
        pow(y, 3);
    dJ(2)(2, 0) = -(d2_global_to_local_cake_core_dadb(pt, z / y, -(x / y), y) /
                    pow(y, 2));
    dJ(2)(2, 1) =
        (-(y * d_global_to_local_cake_core_da(pt, z / y, -(x / y), y)) +
         pow(y, 2) *
             d2_global_to_local_cake_core_dadvar(pt, z / y, -(x / y), y) +
         x * d2_global_to_local_cake_core_dadb(pt, z / y, -(x / y), y) -
         z * d2_global_to_local_cake_core_da2(pt, z / y, -(x / y), y)) /
        pow(y, 3);
    dJ(2)(2, 2) =
        d2_global_to_local_cake_core_da2(pt, z / y, -(x / y), y) / pow(y, 2);
    break;

  case static_cast<int>(patch_piece::minus_y):
    J(0)(0) = 0;
    J(0)(1) = z / pow(y, 2);
    J(0)(2) = -(1 / y);
    J(1)(0) = -(1 / y);
    J(1)(1) = x / pow(y, 2);
    J(1)(2) = 0;
    J(2)(0) = -(d_global_to_local_cake_core_db(pt, -(z / y), -(x / y), y) / y);
    J(2)
    (1) = d_global_to_local_cake_core_dvar(pt, -(z / y), -(x / y), y) +
          (x * d_global_to_local_cake_core_db(pt, -(z / y), -(x / y), y) +
           z * d_global_to_local_cake_core_da(pt, -(z / y), -(x / y), y)) /
              pow(y, 2);
    J(2)(2) = -(d_global_to_local_cake_core_da(pt, -(z / y), -(x / y), y) / y);

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = (-2 * z) / pow(y, 3);
    dJ(0)(1, 2) = pow(y, -2);
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = pow(y, -2);
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = pow(y, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = pow(y, -2);
    dJ(1)(1, 1) = (-2 * x) / pow(y, 3);
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        d2_global_to_local_cake_core_db2(pt, -(z / y), -(x / y), y) / pow(y, 2);
    dJ(2)(0, 1) =
        -((-(y * d_global_to_local_cake_core_db(pt, -(z / y), -(x / y), y)) +
           pow(y, 2) *
               d2_global_to_local_cake_core_dbdvar(pt, -(z / y), -(x / y), y) +
           x * d2_global_to_local_cake_core_db2(pt, -(z / y), -(x / y), y) +
           z * d2_global_to_local_cake_core_dadb(pt, -(z / y), -(x / y), y)) /
          pow(y, 3));
    dJ(2)(0, 2) = d2_global_to_local_cake_core_dadb(pt, -(z / y), -(x / y), y) /
                  pow(y, 2);
    dJ(2)(1, 0) =
        -((-(y * d_global_to_local_cake_core_db(pt, -(z / y), -(x / y), y)) +
           pow(y, 2) *
               d2_global_to_local_cake_core_dbdvar(pt, -(z / y), -(x / y), y) +
           x * d2_global_to_local_cake_core_db2(pt, -(z / y), -(x / y), y) +
           z * d2_global_to_local_cake_core_dadb(pt, -(z / y), -(x / y), y)) /
          pow(y, 3));
    dJ(2)(1, 1) =
        (pow(y, 4) *
             d2_global_to_local_cake_core_dvar2(pt, -(z / y), -(x / y), y) -
         2 * x * y * d_global_to_local_cake_core_db(pt, -(z / y), -(x / y), y) +
         2 * x * pow(y, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, -(z / y), -(x / y), y) +
         pow(x, 2) *
             d2_global_to_local_cake_core_db2(pt, -(z / y), -(x / y), y) -
         2 * y * z * d_global_to_local_cake_core_da(pt, -(z / y), -(x / y), y) +
         2 * pow(y, 2) * z *
             d2_global_to_local_cake_core_dadvar(pt, -(z / y), -(x / y), y) +
         2 * x * z *
             d2_global_to_local_cake_core_dadb(pt, -(z / y), -(x / y), y) +
         pow(z, 2) *
             d2_global_to_local_cake_core_da2(pt, -(z / y), -(x / y), y)) /
        pow(y, 4);
    dJ(2)(1, 2) =
        -((-(y * d_global_to_local_cake_core_da(pt, -(z / y), -(x / y), y)) +
           pow(y, 2) *
               d2_global_to_local_cake_core_dadvar(pt, -(z / y), -(x / y), y) +
           x * d2_global_to_local_cake_core_dadb(pt, -(z / y), -(x / y), y) +
           z * d2_global_to_local_cake_core_da2(pt, -(z / y), -(x / y), y)) /
          pow(y, 3));
    dJ(2)(2, 0) = d2_global_to_local_cake_core_dadb(pt, -(z / y), -(x / y), y) /
                  pow(y, 2);
    dJ(2)(2, 1) =
        -((-(y * d_global_to_local_cake_core_da(pt, -(z / y), -(x / y), y)) +
           pow(y, 2) *
               d2_global_to_local_cake_core_dadvar(pt, -(z / y), -(x / y), y) +
           x * d2_global_to_local_cake_core_dadb(pt, -(z / y), -(x / y), y) +
           z * d2_global_to_local_cake_core_da2(pt, -(z / y), -(x / y), y)) /
          pow(y, 3));
    dJ(2)(2, 2) =
        d2_global_to_local_cake_core_da2(pt, -(z / y), -(x / y), y) / pow(y, 2);
    break;

  case static_cast<int>(patch_piece::plus_z):
    J(0)(0) = -(1 / z);
    J(0)(1) = 0;
    J(0)(2) = x / pow(z, 2);
    J(1)(0) = 0;
    J(1)(1) = 1 / z;
    J(1)(2) = -(y / pow(z, 2));
    J(2)(0) = -(d_global_to_local_cake_core_da(pt, -(x / z), y / z, z) / z);
    J(2)(1) = d_global_to_local_cake_core_db(pt, -(x / z), y / z, z) / z;
    J(2)
    (2) = d_global_to_local_cake_core_dvar(pt, -(x / z), y / z, z) +
          (-(y * d_global_to_local_cake_core_db(pt, -(x / z), y / z, z)) +
           x * d_global_to_local_cake_core_da(pt, -(x / z), y / z, z)) /
              pow(z, 2);

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = pow(z, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = pow(z, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = (-2 * x) / pow(z, 3);
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = 0;
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = 0;
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = -pow(z, -2);
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = -pow(z, -2);
    dJ(1)(2, 2) = (2 * y) / pow(z, 3);
    dJ(2)(0, 0) =
        d2_global_to_local_cake_core_da2(pt, -(x / z), y / z, z) / pow(z, 2);
    dJ(2)(0, 1) = -(d2_global_to_local_cake_core_dadb(pt, -(x / z), y / z, z) /
                    pow(z, 2));
    dJ(2)(0, 2) =
        (z * (d_global_to_local_cake_core_da(pt, -(x / z), y / z, z) -
              z * d2_global_to_local_cake_core_dadvar(pt, -(x / z), y / z, z)) +
         y * d2_global_to_local_cake_core_dadb(pt, -(x / z), y / z, z) -
         x * d2_global_to_local_cake_core_da2(pt, -(x / z), y / z, z)) /
        pow(z, 3);
    dJ(2)(1, 0) = -(d2_global_to_local_cake_core_dadb(pt, -(x / z), y / z, z) /
                    pow(z, 2));
    dJ(2)(1, 1) =
        d2_global_to_local_cake_core_db2(pt, -(x / z), y / z, z) / pow(z, 2);
    dJ(2)(1, 2) =
        (-(z * d_global_to_local_cake_core_db(pt, -(x / z), y / z, z)) +
         pow(z, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, -(x / z), y / z, z) -
         y * d2_global_to_local_cake_core_db2(pt, -(x / z), y / z, z) +
         x * d2_global_to_local_cake_core_dadb(pt, -(x / z), y / z, z)) /
        pow(z, 3);
    dJ(2)(2, 0) =
        (z * (d_global_to_local_cake_core_da(pt, -(x / z), y / z, z) -
              z * d2_global_to_local_cake_core_dadvar(pt, -(x / z), y / z, z)) +
         y * d2_global_to_local_cake_core_dadb(pt, -(x / z), y / z, z) -
         x * d2_global_to_local_cake_core_da2(pt, -(x / z), y / z, z)) /
        pow(z, 3);
    dJ(2)(2, 1) =
        (-(z * d_global_to_local_cake_core_db(pt, -(x / z), y / z, z)) +
         pow(z, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, -(x / z), y / z, z) -
         y * d2_global_to_local_cake_core_db2(pt, -(x / z), y / z, z) +
         x * d2_global_to_local_cake_core_dadb(pt, -(x / z), y / z, z)) /
        pow(z, 3);
    dJ(2)(2, 2) =
        (pow(z, 4) *
             d2_global_to_local_cake_core_dvar2(pt, -(x / z), y / z, z) +
         2 * y * z * d_global_to_local_cake_core_db(pt, -(x / z), y / z, z) -
         2 * y * pow(z, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, -(x / z), y / z, z) +
         pow(y, 2) * d2_global_to_local_cake_core_db2(pt, -(x / z), y / z, z) -
         2 * x * z * d_global_to_local_cake_core_da(pt, -(x / z), y / z, z) +
         2 * x * pow(z, 2) *
             d2_global_to_local_cake_core_dadvar(pt, -(x / z), y / z, z) -
         2 * x * y * d2_global_to_local_cake_core_dadb(pt, -(x / z), y / z, z) +
         pow(x, 2) * d2_global_to_local_cake_core_da2(pt, -(x / z), y / z, z)) /
        pow(z, 4);
    break;

  case static_cast<int>(patch_piece::minus_z):
    J(0)(0) = -(1 / z);
    J(0)(1) = 0;
    J(0)(2) = x / pow(z, 2);
    J(1)(0) = 0;
    J(1)(1) = -(1 / z);
    J(1)(2) = y / pow(z, 2);
    J(2)(0) = -(d_global_to_local_cake_core_da(pt, -(x / z), -(y / z), z) / z);
    J(2)(1) = -(d_global_to_local_cake_core_db(pt, -(x / z), -(y / z), z) / z);
    J(2)
    (2) = d_global_to_local_cake_core_dvar(pt, -(x / z), -(y / z), z) +
          (y * d_global_to_local_cake_core_db(pt, -(x / z), -(y / z), z) +
           x * d_global_to_local_cake_core_da(pt, -(x / z), -(y / z), z)) /
              pow(z, 2);

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = pow(z, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = pow(z, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = (-2 * x) / pow(z, 3);
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = 0;
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = 0;
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = pow(z, -2);
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = pow(z, -2);
    dJ(1)(2, 2) = (-2 * y) / pow(z, 3);
    dJ(2)(0, 0) =
        d2_global_to_local_cake_core_da2(pt, -(x / z), -(y / z), z) / pow(z, 2);
    dJ(2)(0, 1) = d2_global_to_local_cake_core_dadb(pt, -(x / z), -(y / z), z) /
                  pow(z, 2);
    dJ(2)(0, 2) =
        -((-(z * d_global_to_local_cake_core_da(pt, -(x / z), -(y / z), z)) +
           pow(z, 2) *
               d2_global_to_local_cake_core_dadvar(pt, -(x / z), -(y / z), z) +
           y * d2_global_to_local_cake_core_dadb(pt, -(x / z), -(y / z), z) +
           x * d2_global_to_local_cake_core_da2(pt, -(x / z), -(y / z), z)) /
          pow(z, 3));
    dJ(2)(1, 0) = d2_global_to_local_cake_core_dadb(pt, -(x / z), -(y / z), z) /
                  pow(z, 2);
    dJ(2)(1, 1) =
        d2_global_to_local_cake_core_db2(pt, -(x / z), -(y / z), z) / pow(z, 2);
    dJ(2)(1, 2) =
        -((-(z * d_global_to_local_cake_core_db(pt, -(x / z), -(y / z), z)) +
           pow(z, 2) *
               d2_global_to_local_cake_core_dbdvar(pt, -(x / z), -(y / z), z) +
           y * d2_global_to_local_cake_core_db2(pt, -(x / z), -(y / z), z) +
           x * d2_global_to_local_cake_core_dadb(pt, -(x / z), -(y / z), z)) /
          pow(z, 3));
    dJ(2)(2, 0) =
        -((-(z * d_global_to_local_cake_core_da(pt, -(x / z), -(y / z), z)) +
           pow(z, 2) *
               d2_global_to_local_cake_core_dadvar(pt, -(x / z), -(y / z), z) +
           y * d2_global_to_local_cake_core_dadb(pt, -(x / z), -(y / z), z) +
           x * d2_global_to_local_cake_core_da2(pt, -(x / z), -(y / z), z)) /
          pow(z, 3));
    dJ(2)(2, 1) =
        -((-(z * d_global_to_local_cake_core_db(pt, -(x / z), -(y / z), z)) +
           pow(z, 2) *
               d2_global_to_local_cake_core_dbdvar(pt, -(x / z), -(y / z), z) +
           y * d2_global_to_local_cake_core_db2(pt, -(x / z), -(y / z), z) +
           x * d2_global_to_local_cake_core_dadb(pt, -(x / z), -(y / z), z)) /
          pow(z, 3));
    dJ(2)(2, 2) =
        (pow(z, 4) *
             d2_global_to_local_cake_core_dvar2(pt, -(x / z), -(y / z), z) -
         2 * y * z * d_global_to_local_cake_core_db(pt, -(x / z), -(y / z), z) +
         2 * y * pow(z, 2) *
             d2_global_to_local_cake_core_dbdvar(pt, -(x / z), -(y / z), z) +
         pow(y, 2) *
             d2_global_to_local_cake_core_db2(pt, -(x / z), -(y / z), z) -
         2 * x * z * d_global_to_local_cake_core_da(pt, -(x / z), -(y / z), z) +
         2 * x * pow(z, 2) *
             d2_global_to_local_cake_core_dadvar(pt, -(x / z), -(y / z), z) +
         2 * x * y *
             d2_global_to_local_cake_core_dadb(pt, -(x / z), -(y / z), z) +
         pow(x, 2) *
             d2_global_to_local_cake_core_da2(pt, -(x / z), -(y / z), z)) /
        pow(z, 4);
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("No jacobians available for patch %s",
                piece_name(static_cast<patch_piece>(patch)).c_str());
#else
    assert(0);
#endif
    break;
  }

  return std_make_tuple(J, dJ);
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_JACOBIANS_HXX
