#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, POS, INT>() const;
template void BoundaryCondition::apply_on_face<INT, POS, INT>() const;
template void BoundaryCondition::apply_on_face<POS, POS, INT>() const;

} // namespace CarpetX
