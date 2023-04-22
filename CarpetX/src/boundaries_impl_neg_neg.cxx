#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, NEG, NEG>() const;
template void BoundaryCondition::apply_on_face<INT, NEG, NEG>() const;
template void BoundaryCondition::apply_on_face<POS, NEG, NEG>() const;

} // namespace CarpetX
