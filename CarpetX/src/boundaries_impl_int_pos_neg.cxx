#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<INT, POS, NEG>() const;

} // namespace CarpetX
