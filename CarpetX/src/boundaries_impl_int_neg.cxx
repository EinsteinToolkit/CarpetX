#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, INT, NEG>() const;
template void BoundaryCondition::apply_on_face<INT, INT, NEG>() const;
template void BoundaryCondition::apply_on_face<POS, INT, NEG>() const;

} // namespace CarpetX
