#ifndef CARPETX_SUBCYCLING_DENSE_CONTROL_BUILDER_HXX
#define CARPETX_SUBCYCLING_DENSE_CONTROL_BUILDER_HXX

#include "subcycling_dense_mfab_state.hxx"
#include "subcycling_dense_stencil.hxx"

#include <cstddef>
#include <memory>
#include <vector>

namespace CarpetX {

struct DenseMFabRawSampleView {
  DenseSampleConstraint constraint;
  const OwnedMultiFabDenseState *state;
};

class DenseMFabControlSet {
public:
  SubcyclingODEMethod method() const noexcept;
  double parent_dt() const noexcept;
  std::size_t control_count() const noexcept;
  const OwnedMultiFabDenseState &control(std::size_t index) const;

  std::vector<std::unique_ptr<OwnedMultiFabDenseState>>
  release_controls() &&;

  DenseMFabControlSet(const DenseMFabControlSet &) = delete;
  DenseMFabControlSet &operator=(const DenseMFabControlSet &) = delete;
  DenseMFabControlSet(DenseMFabControlSet &&) noexcept = default;
  DenseMFabControlSet &operator=(DenseMFabControlSet &&) noexcept = default;

private:
  friend DenseMFabControlSet build_reference_dense_controls(
      SubcyclingODEMethod, double,
      const std::vector<DenseMFabRawSampleView> &);

  DenseMFabControlSet(
      SubcyclingODEMethod method, double parent_dt,
      std::vector<std::unique_ptr<OwnedMultiFabDenseState>> controls);

  SubcyclingODEMethod method_;
  double parent_dt_;
  std::vector<std::unique_ptr<OwnedMultiFabDenseState>> controls_;
};

DenseMFabControlSet build_reference_dense_controls(
    SubcyclingODEMethod method, double parent_dt,
    const std::vector<DenseMFabRawSampleView> &samples);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_DENSE_CONTROL_BUILDER_HXX
