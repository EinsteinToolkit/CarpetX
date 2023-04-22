#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, NEG, POS>() const;
template void BoundaryCondition::apply_on_face<INT, NEG, POS>() const;
template void BoundaryCondition::apply_on_face<POS, NEG, POS>() const;

} // namespace CarpetX
