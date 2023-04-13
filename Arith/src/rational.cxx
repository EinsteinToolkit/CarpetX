#include "rational.hxx"

#include <numeric>

namespace Arith {
using namespace std;

template <typename I> constexpr bool isgood(const rational<I> &x) {
  return x.den >= 0 && gcd(x.num, x.den) == 1;
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestRational() {
  using rat64 = rational<int64_t>;
  static_assert(rat64().num == 0);
  static_assert(rat64(1).num == 1);
  static_assert(rat64(1, 2).num == 1);
  static_assert(rat64(2, 4).num == 1);
  static_assert(rat64().den == 1);
  static_assert(rat64(1).den == 1);
  static_assert(rat64(1, 2).den == 2);
  static_assert(rat64(2, 4).den == 2);

  constexpr auto n = rat64();
  constexpr auto e = rat64(1);
  constexpr auto x = rat64(1, 2);
  constexpr auto y = rat64(2, 3);
  constexpr auto z = rat64(3, 4);
  constexpr auto a = 4;

  +x;
  -x;

  x + y;
  x - y;
  x *y;
  x / y;

  a + y;
  a - y;
  a *y;
  a / y;

  x + a;
  x - a;
  x *a;
  x / a;

  auto r = x;
  r = x;
  r += y;
  r = x;
  r -= y;
  r = x;
  r *= y;
  r = x;
  r /= y;

  r = x;
  r += a;
  r = x;
  r -= a;
  r = x;
  r *= a;
  r = x;
  r /= a;

  x == y;
  x != y;
  x < y;
  x > y;
  x <= y;
  x >= y;

  a == y;
  a != y;
  a < y;
  a > y;
  a <= y;
  a >= y;

  x == a;
  x != a;
  x < a;
  x > a;
  x <= a;
  x >= a;

  using std::abs;
  static_assert(isgood(abs(x)));

  using std::max, std::min;
  static_assert(isgood(max(x, y)));
  static_assert(isgood(min(x, y)));

  using std::pow;
  static_assert(isgood(pown(x, a)));

  static_assert(+x == n + x);
  static_assert(x + n == x);
  static_assert(n + x == x);
  static_assert(x + (y + z) == (x + y) + z);
  static_assert(x + y == y + x);

  static_assert(-x == n - x);
  static_assert(-(-x) == x);
  static_assert(x + (-y) == x - y);
  static_assert(x - x == n);

  static_assert(x * e == x);
  static_assert(e * x == x);
  static_assert(x * (y * z) == (x * y) * z);
  static_assert(x * y == y * x);
  static_assert(n * x == n);

  static_assert(1 / (1 / x) == x);
  static_assert(x / x == e);

  static_assert(x * (y + z) == x * y + x * z);
}

} // namespace Arith
