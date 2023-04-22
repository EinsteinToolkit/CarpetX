#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, INT, POS>() const;
template void BoundaryCondition::apply_on_face<INT, INT, POS>() const;
template void BoundaryCondition::apply_on_face<POS, INT, POS>() const;

} // namespace CarpetX
