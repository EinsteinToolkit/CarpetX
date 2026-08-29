#include "checked.hxx"

namespace Arith {

void TestChecked() {
  static_assert(checked_pos(+1) == +(+1));
  static_assert(checked_pos(-1) == +(-1));

  static_assert(checked_neg(+1) == -(+1));
  static_assert(checked_neg(-1) == -(-1));

  static_assert(checked_add(+2, +3) == (+2) + (+3));
  static_assert(checked_add(+2, -3) == (+2) + (-3));
  static_assert(checked_add(-2, +3) == (-2) + (+3));
  static_assert(checked_add(-2, -3) == (-2) + (-3));

  static_assert(checked_sub(+2, +3) == (+2) - (+3));
  static_assert(checked_sub(+2, -3) == (+2) - (-3));
  static_assert(checked_sub(-2, +3) == (-2) - (+3));
  static_assert(checked_sub(-2, -3) == (-2) - (-3));

  static_assert(checked_mul(+2, +3) == (+2) * (+3));
  static_assert(checked_mul(+2, -3) == (+2) * (-3));
  static_assert(checked_mul(-2, +3) == (-2) * (+3));
  static_assert(checked_mul(-2, -3) == (-2) * (-3));

  static_assert(checked_div(+2, +3) == (+2) / (+3));
  static_assert(checked_div(+2, -3) == (+2) / (-3));
  static_assert(checked_div(-2, +3) == (-2) / (+3));
  static_assert(checked_div(-2, -3) == (-2) / (-3));

  static_assert(checked_mod(+2, +3) == (+2) % (+3));
  static_assert(checked_mod(+2, -3) == (+2) % (-3));
  static_assert(checked_mod(-2, +3) == (-2) % (+3));
  static_assert(checked_mod(-2, -3) == (-2) % (-3));
}

} // namespace Arith
