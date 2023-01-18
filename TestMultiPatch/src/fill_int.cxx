#include "test_multipatch.hxx"

extern "C" void TestMultiPatch::TestMultiPatch_fill_int(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestMultiPatch_fill_int;
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

  grid.loop_int_device<0, 0, 0>(grid.nghostzones, loop_lambda);
}
