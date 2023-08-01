#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, POS, NEG>() const;

} // namespace CarpetX
