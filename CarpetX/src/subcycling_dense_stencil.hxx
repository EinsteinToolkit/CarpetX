#ifndef CARPETX_SUBCYCLING_DENSE_STENCIL_HXX
#define CARPETX_SUBCYCLING_DENSE_STENCIL_HXX

#include "subcycling_step_context.hxx"

#include <cstddef>
#include <vector>

namespace CarpetX {

enum class DenseSampleKind { value, scaled_derivative };

struct DenseSampleConstraint {
  double theta;
  DenseSampleKind kind;
};

struct DenseStencilSpecification {
  SubcyclingODEMethod method;
  int endpoint_order;
  int dense_uniform_order;
  int stage_count;
  int extra_rhs_evaluations;
  std::vector<DenseSampleConstraint> constraints;
};

class DenseStencil {
public:
  explicit DenseStencil(DenseStencilSpecification specification);

  const DenseStencilSpecification &specification() const noexcept;
  std::size_t control_count() const noexcept;
  std::size_t sample_count() const noexcept;
  const std::vector<double> &weights() const noexcept;
  double weight(std::size_t control, std::size_t sample) const;
  std::vector<double>
  make_controls(const std::vector<double> &samples) const;

private:
  static std::vector<double>
  build_weights(const DenseStencilSpecification &specification);

  const DenseStencilSpecification specification_;
  const std::vector<double> weights_;
};

const DenseStencil &reference_dense_stencil(SubcyclingODEMethod method);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_DENSE_STENCIL_HXX
