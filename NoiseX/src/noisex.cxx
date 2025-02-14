#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "noisex.hxx"

#include <cmath>
#include <random>

namespace NoiseX {

struct NoiseData {
  std::default_random_engine default_engine{};
  std::mt19937 mt_engine{};
  std::uniform_real_distribution<CCTK_REAL> distrib{};
};

static NoiseData noise_data{};

extern "C" int NoiseX_setup_globals() {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(seed_type, "hardware random device")) {
    std::random_device device{};
    noise_data.default_engine = std::default_random_engine(device());
    noise_data.mt_engine = std::mt19937(device());
  } else if (CCTK_EQUALS(seed_type, "fixed")) {
    noise_data.default_engine = std::default_random_engine(fixed_seed_value);
    noise_data.mt_engine = std::mt19937(fixed_seed_value);
  } else {
    CCTK_VERROR("Unrecognized seed type \"%s\"", seed_type);
  }

  noise_data.distrib = std::uniform_real_distribution<CCTK_REAL>(
      -noise_amplitude, noise_amplitude);

  return 0;
}

void add(const cGH *restrict const cctkGH, Loop::GF3D2<CCTK_REAL> var) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GridDescBaseDevice grid(cctkGH);

  if (CCTK_EQUALS(random_engine, "default")) {
    grid.loop_int<0, 0, 0>(grid.nghostzones,
                           [&] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 var(p.I) += noise_data.distrib(
                                     noise_data.default_engine);
                               });

  } else if (CCTK_EQUALS(random_engine, "mersenne twister")) {
    grid.loop_int<0, 0, 0>(grid.nghostzones,
                           [&] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 var(p.I) +=
                                     noise_data.distrib(noise_data.mt_engine);
                               });

  } else if (CCTK_EQUALS(random_engine, "sine wave")) {
    grid.loop_int<0, 0, 0>(grid.nghostzones,
                           [=] CCTK_DEVICE CCTK_HOST(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 using std::sin;
                                 var(p.I) += noise_amplitude *
                                             sin(noise_frequency * p.x / p.dx) *
                                             sin(noise_frequency * p.y / p.dy) *
                                             sin(noise_frequency * p.z / p.dz);
                               });

  } else {
    CCTK_VERROR("Unrecognized random engine \"%s\"", random_engine);
  }
}

} // namespace NoiseX