#ifndef CARPETX_SUBCYCLING_DENSE_OUTPUT_HXX
#define CARPETX_SUBCYCLING_DENSE_OUTPUT_HXX

#include "subcycling_step_context.hxx"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace CarpetX {

using TableauFingerprint = std::array<std::uint8_t, 32>;

struct DenseCapability {
  SubcyclingODEMethod method;
  TableauFingerprint tableau_fingerprint;
  int endpoint_order;
  int dense_uniform_order;
  int stage_count;
  int extra_rhs_evaluations;
  int persistent_vector_count;
  bool arbitrary_theta;
  bool verified;
};

struct DenseIntervalId {
  int level;
  step_clock_t begin_clock;
  step_clock_t end_clock;
  double begin_time;
  double end_time;
  SubcyclingODEMethod method;
  TableauFingerprint tableau_fingerprint;
};

class DenseStateVector {
public:
  virtual ~DenseStateVector() = default;
  virtual bool compatible(const DenseStateVector &other) const noexcept = 0;
  virtual void copy_from(const DenseStateVector &other) = 0;
  virtual void linear_combination(
      const std::vector<double> &weights,
      const std::vector<const DenseStateVector *> &sources) = 0;
};

class DenseInterval {
public:
  const DenseIntervalId &id() const noexcept;
  const DenseCapability &capability() const noexcept;
  std::size_t control_count() const noexcept;
  void evaluate(double theta, DenseStateVector &destination) const;
  void evaluate(step_clock_t theta, DenseStateVector &destination) const;

private:
  friend class DenseIntervalBuilder;
  DenseInterval(DenseCapability capability, DenseIntervalId id,
                std::vector<std::unique_ptr<DenseStateVector>> controls);

  DenseCapability capability_;
  DenseIntervalId id_;
  std::vector<std::unique_ptr<DenseStateVector>> controls_;
};

class DenseIntervalBuilder {
public:
  DenseIntervalBuilder(DenseCapability capability, DenseIntervalId id);
  ~DenseIntervalBuilder() = default;

  DenseIntervalBuilder(const DenseIntervalBuilder &) = delete;
  DenseIntervalBuilder &operator=(const DenseIntervalBuilder &) = delete;
  DenseIntervalBuilder(DenseIntervalBuilder &&other) noexcept;
  DenseIntervalBuilder &operator=(DenseIntervalBuilder &&other);

  void add_control(std::unique_ptr<DenseStateVector> control);
  std::shared_ptr<const DenseInterval> seal();

private:
  enum class State { active, sealed, moved_from };
  void require_active() const;

  DenseCapability capability_{};
  DenseIntervalId id_{};
  std::vector<std::unique_ptr<DenseStateVector>> controls_;
  State state_{State::active};
};

class DenseOutputProvider final {
public:
  explicit DenseOutputProvider(DenseCapability capability);
  const DenseCapability &capability() const noexcept;
  std::unique_ptr<DenseIntervalBuilder>
  begin_interval(const DenseIntervalId &id) const;

private:
  const DenseCapability capability_;
};

class DenseOutputRegistry {
public:
  void register_provider(std::shared_ptr<const DenseOutputProvider> provider);
  std::shared_ptr<const DenseOutputProvider>
  require(SubcyclingODEMethod method,
          const TableauFingerprint &tableau_fingerprint) const;

private:
  struct Entry {
    DenseCapability capability;
    std::shared_ptr<const DenseOutputProvider> provider;
  };
  std::vector<Entry> entries_;
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_DENSE_OUTPUT_HXX
