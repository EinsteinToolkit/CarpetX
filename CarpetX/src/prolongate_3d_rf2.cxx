#include "prolongate_3d_rf2.hxx"

#include <cassert>
#include <iostream>

namespace CarpetX {

std::ostream &operator<<(std::ostream &os, const centering_t cent) {
  switch (cent) {
  case centering_t::vertex:
    return os << "vertex";
  case centering_t::cell:
    return os << "cell";
  }
  assert(false);
}

std::ostream &operator<<(std::ostream &os, const interpolation_t intp) {
  switch (intp) {
  case interpolation_t::poly:
    return os << "poly";
  case interpolation_t::hermite:
    return os << "hermite";
  case interpolation_t::cons:
    return os << "cons";
  case interpolation_t::eno:
    return os << "eno";
  case interpolation_t::eno_star:
    return os << "eno_star";
  case interpolation_t::minmod:
    return os << "minmod";
  }
  assert(false);
}

std::ostream &operator<<(std::ostream &os, const fallback_t fb) {
  switch (fb) {
  case fallback_t::none:
    return os << "none";
  case fallback_t::constant:
    return os << "constant";
  case fallback_t::linear:
    return os << "linear";
  }
  assert(false);
}

} // namespace CarpetX
