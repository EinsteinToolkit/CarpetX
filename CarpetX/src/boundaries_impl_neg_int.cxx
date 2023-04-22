#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, NEG, INT>() const;
template void BoundaryCondition::apply_on_face<INT, NEG, INT>() const;
template void BoundaryCondition::apply_on_face<POS, NEG, INT>() const;

} // namespace CarpetX
