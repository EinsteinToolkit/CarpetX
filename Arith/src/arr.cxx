#include "arr.hxx"

#include <array>
#include <cassert>

namespace Arith {

void TestArr() {
  constexpr arr<char, 0> c0{};
  constexpr arr<int, 1> i1{1};
  constexpr arr<float, 2> f2{2, 3};
  constexpr arr<double, 4> d4{4, 5, 6, 7};

  static_assert(i1.at(0) == 1);
  static_assert(f2.at(0) == 2);
  static_assert(f2.at(1) == 3);
  static_assert(d4.at(0) == 4);
  static_assert(d4.at(1) == 5);
  static_assert(d4.at(2) == 6);
  static_assert(d4.at(3) == 7);

  static_assert(i1[0] == 1);
  static_assert(f2[0] == 2);
  static_assert(f2[1] == 3);
  static_assert(d4[0] == 4);
  static_assert(d4[1] == 5);
  static_assert(d4[2] == 6);
  static_assert(d4[3] == 7);

  static_assert(*d4.begin() == 4);
  static_assert(*d4.rbegin() == 7);
  static_assert(*(d4.end() - 1) == 7);
  static_assert(*(d4.rend() - 1) == 4);

  static_assert(&*d4.begin() == d4.data() + 0);
  static_assert(&*d4.rbegin() == d4.data() + 3);
  static_assert(&*d4.end() == d4.data() + 4);
  // static_assert(&*d4.rend() == d4.data() - 1);
  static_assert(&*(d4.begin() + 1) == d4.data() + 1);
  static_assert(&*(d4.rbegin() + 1) == d4.data() + 2);
  static_assert(&*(d4.end() - 1) == d4.data() + 3);
  static_assert(&*(d4.rend() - 1) == d4.data() + 0);
  static_assert(d4.begin()[1] == d4[1]);
  static_assert(d4.rbegin()[1] == d4[2]);
  static_assert(d4.end()[-1] == d4[3]);
  static_assert(d4.rend()[-1] == d4[0]);

  static_assert(d4.front() == 4);
  static_assert(d4.back() == 7);

  static_assert(c0.empty());
  static_assert(!i1.empty());
  static_assert(!f2.empty());
  static_assert(!d4.empty());

  static_assert(c0.size() == 0);
  static_assert(i1.size() == 1);
  static_assert(f2.size() == 2);
  static_assert(d4.size() == 4);

  static_assert(c0.max_size() == 0);
  static_assert(i1.max_size() == 1);
  static_assert(f2.max_size() == 2);
  static_assert(d4.max_size() == 4);
}

template class arr<char, 0>;
template class arr<int, 1>;
template class arr<float, 2>;
template class arr<double, 4>;

} // namespace Arith
