#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace TestBoxInBox {

extern "C" void TestBoxInBox_Update(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestBoxInBox_Update;

  position_x[0] = -0.5 + 0.1 * cctk_iteration;
}

} // namespace BoxInBox
