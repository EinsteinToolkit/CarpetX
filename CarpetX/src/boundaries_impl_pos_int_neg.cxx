#include "boundaries_impl.hxx"

namespace CarpetX {

template void BoundaryCondition::apply_on_face<POS, INT, NEG>() const;

} // namespace CarpetX
