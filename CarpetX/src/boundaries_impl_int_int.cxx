#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<NEG, INT, INT>() const;
//  template void BoundaryCondition::apply_on_face<INT, INT, INT>() const;
template void BoundaryCondition::apply_on_face<POS, INT, INT>() const;

} // namespace CarpetX
