#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, POS, POS>() const;
template void BoundaryCondition::apply_on_face<INT, POS, POS>() const;
template void BoundaryCondition::apply_on_face<POS, POS, POS>() const;

} // namespace CarpetX
