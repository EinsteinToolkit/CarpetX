#include <loop.hxx>
#include <loop_device.hxx>

#include <vec.hxx>
#include <arith.hxx>

#include <fixmath.hxx>
#include <cctk.h>

#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

namespace TestMultiPatch {

/**
 * Stores a 3-vector
 */
template <typename real_t> using svec = Arith::vec<real_t, Loop::dim>;

/**
 * Stores a 4-vector
 */
template <typename real_t> using qvec = Arith::vec<real_t, 4>;

template <typename real_t>
CCTK_DEVICE CCTK_HOST static inline real_t
plane_wave_omega(const svec<real_t> &wave_number) {
  using std::sqrt;
  return sqrt(wave_number(0) * wave_number(0) +
              wave_number(1) * wave_number(1) +
              wave_number(2) * wave_number(2));
}

template <typename real_t>
CCTK_DEVICE CCTK_HOST static inline qvec<real_t>
plane_wave_n(const svec<real_t> &wave_number, const qvec<real_t> &offsets,
             const qvec<real_t> coords) {
  return qvec<real_t>{plane_wave_omega(wave_number) * (coords(0) - offsets(0)),
                      wave_number(0) * (coords(1) - offsets(1)),
                      wave_number(1) * (coords(2) - offsets(2)),
                      wave_number(2) * (coords(3) - offsets(3))};
}

template <typename real_t>
CCTK_DEVICE CCTK_HOST static inline real_t
plane_wave_w(const svec<real_t> &wave_number, const qvec<real_t> &offsets,
             const qvec<real_t> coords) {
  const auto n = plane_wave_n(wave_number, offsets, coords);
  return 2 * M_PI * (n(0) + n(1) + n(2) + n(3));
}

extern "C" void TestMultiPatch_initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestMultiPatch_initialize;
  DECLARE_CCTK_PARAMETERS;

  using std::cos;

  const Loop::GridDescBaseDevice grid(cctkGH);

  const std::array<int, Loop::dim> indextype_vc = {0, 0, 0};
  const Loop::GF3D2layout layout_vc(cctkGH, indextype_vc);

  const Loop::GF3D2<const CCTK_REAL> gf_vcoordx(layout_vc, vcoordx);
  const Loop::GF3D2<const CCTK_REAL> gf_vcoordy(layout_vc, vcoordy);
  const Loop::GF3D2<const CCTK_REAL> gf_vcoordz(layout_vc, vcoordz);

  const Loop::GF3D2<CCTK_REAL> gf_test_gf(layout_vc, test_gf);

  const svec<CCTK_REAL> wave_numbers{wave_number[0], wave_number[1],
                                     wave_number[2]};
  const qvec<CCTK_REAL> offsets{time_offset, space_offset[0], space_offset[1],
                                space_offset[2]};

  const auto loop_lambda =
      [=] ARITH_DEVICE(const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index(layout_vc, p.I);

        const qvec<CCTK_REAL> coords{
            CCTK_REAL{0},
            gf_vcoordx(index),
            gf_vcoordy(index),
            gf_vcoordz(index),
        };

        gf_test_gf(index) = cos(plane_wave_w(wave_numbers, offsets, coords));
      };

  grid.loop_all_device<0, 0, 0>(grid.nghostzones, loop_lambda);
}

} // namespace TestMultiPatch
