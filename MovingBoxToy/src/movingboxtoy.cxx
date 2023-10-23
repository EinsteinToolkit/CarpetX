#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

template <typename T>
constexpr auto pow2(T val) noexcept -> T {
    return val * val;
}

extern "C" void MovingBoxToy_MoveBoxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MovingBoxToy_MoveBoxes;
  DECLARE_CCTK_PARAMETERS;

  using std::cos;

  const CCTK_REAL omega{M_PI/4};

  // Initial positions of box 1
  const auto x0_1{position_x_1};

  // Initial positions of box 2
  const auto x0_2{position_x_2};

  // Trajectory of box 1
  position_x[0] = x0_1 * cos(omega * cctk_time);
  position_y[0] = x0_1 * sin(omega * cctk_time);

  // Trajectory of box 2
  position_x[1] = x0_2 * cos(omega * cctk_time);
  position_y[1] = x0_2 * sin(omega * cctk_time);
}