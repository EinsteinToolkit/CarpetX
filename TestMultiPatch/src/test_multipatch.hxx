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
 * Stores a 3-vector with all indices up
 */
template <typename real_t>
using svec = Arith::vec<real_t, Loop::dim>;

/**
 * Stores a 4-vector with all indices up
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
  return qvec<real_t>{plane_wave_omega(wave_number) *
                            (coords(0) - offsets(0)),
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

extern "C" void TestMultiPatch_fill_all(CCTK_ARGUMENTS);
extern "C" void TestMultiPatch_fill_int(CCTK_ARGUMENTS);

} // namespace TestMultiPatch
