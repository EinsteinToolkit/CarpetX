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
  using std::sqrt;

  const CCTK_REAL omega{M_PI/4};

  // Radius of each box
  const auto R1 = sqrt(pow2(position_x[0]) + pow2(position_y[0]) + pow2(position_z[0]));
  const auto R2 = sqrt(pow2(position_x[1]) + pow2(position_y[1]) + pow2(position_z[1])); 

  // Trajectory of box 1
  position_x[0] = R1 * cos(omega * cctk_time);
  position_y[0] = R1 * sin(omega * cctk_time);

  // Trajectory of box 2
  position_x[1] = -R2 * cos(omega * cctk_time);
  position_y[1] = -R2 * sin(omega * cctk_time);
}