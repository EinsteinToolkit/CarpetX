#ifndef CARPETX_SUBCYCLING_NATIVE_GATE_HXX
#define CARPETX_SUBCYCLING_NATIVE_GATE_HXX

#include "subcycling_step_context.hxx"

#include <cctk.h>

#include <cstddef>
#include <memory>

namespace CarpetX {

struct NativeGateMethodContract {
  SubcyclingODEMethod method;
  const char *ode_parameter_value;
  std::size_t extra_rhs_evaluations;
  std::size_t control_count;
};

namespace native_gate_detail {

// Kept independent of driver types so the authoritative patch clock envelope
// can be regression-tested without constructing a Cactus hierarchy.
struct PatchContextInput {
  int active_min_level;
  int active_max_level;
  int active_min_patch;
  int active_max_patch;
  int iteration;
  int time_refinement_factor;
  step_clock_t level_iteration;
  step_clock_t delta_iteration;
  double time;
  double delta_time;
};

[[nodiscard]] StepContext
make_patch_step_context(const PatchContextInput &input,
                        SubcyclingODEMethod method);

} // namespace native_gate_detail

[[nodiscard]] NativeGateMethodContract
native_gate_method_contract(int evolution_iteration);

struct NativeGateReceipt;

class NativeGateSession final {
public:
  NativeGateSession() noexcept;
  ~NativeGateSession();
  NativeGateSession(const NativeGateSession &) = delete;
  NativeGateSession &operator=(const NativeGateSession &) = delete;
  NativeGateSession(NativeGateSession &&) noexcept;
  NativeGateSession &operator=(NativeGateSession &&) noexcept;
  bool active() const noexcept;

private:
  struct Storage;
  explicit NativeGateSession(std::unique_ptr<Storage> storage) noexcept;
  std::unique_ptr<Storage> storage_;
  friend NativeGateSession begin_native_gate(cGH *, int);
  friend struct NativeGateReceipt;
  friend NativeGateReceipt end_native_gate(NativeGateSession &&);
};

struct NativeGateReceipt {
  SubcyclingODEMethod method;
  std::size_t extra_rhs_evaluations;
  std::size_t control_count;
};

// Inventory mode is deliberately non-promoting. The path and build provenance
// are supplied through CARPETX_NATIVE_GATE_EXPECTATION.
void write_native_gate_inventory(cGH *cctkGH, int cctk_itlast);

// Native mode consumes the frozen inventory, creates the certified transaction
// and installs the StepContext for the immediately following ODESolvers call.
[[nodiscard]] NativeGateSession begin_native_gate(cGH *cctkGH,
                                                  int cctk_itlast);
[[nodiscard]] NativeGateReceipt end_native_gate(NativeGateSession &&session);

} // namespace CarpetX

#endif
