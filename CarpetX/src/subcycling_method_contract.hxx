#ifndef CARPETX_SUBCYCLING_METHOD_CONTRACT_HXX
#define CARPETX_SUBCYCLING_METHOD_CONTRACT_HXX

#include "subcycling_dense_output.hxx"
#include "subcycling_scratch_state_transaction.hxx"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace CarpetX {

struct SubcyclingMethodContract {
  std::string provider_id;
  std::string parameter_value;
  DenseCapability dense;
};

struct SubcyclingGroupSchema {
  std::vector<ScratchGroupPair> ordered_group_pairs;
  std::vector<int> dependent_groups;
};

struct SubcyclingMethodContractSnapshot {
  SubcyclingMethodContract contract;
  std::optional<SubcyclingGroupSchema> group_schema;
  std::uint64_t registration_epoch{0};
};

class SubcyclingMethodContractRegistry final {
public:
  bool empty() const noexcept;

  std::uint64_t
  register_contract(const SubcyclingMethodContract &contract);
  std::uint64_t
  register_group_schema(std::string_view provider_id,
                        const SubcyclingGroupSchema &group_schema);

  SubcyclingMethodContractSnapshot require_snapshot() const;
  SubcyclingMethodContractSnapshot require_complete_snapshot() const;
  void
  validate_snapshot(const SubcyclingMethodContractSnapshot &snapshot) const;

private:
  std::optional<SubcyclingMethodContractSnapshot> snapshot_;
};

// These functions bind the pure registry above to the active CarpetX GHExt.
// Registration happens after GHExt construction at CCTK_PARAMCHECK. Returned
// snapshots are immutable copies suitable for a future HierarchyStepper
// session; validate them immediately before touching primary state.
std::uint64_t
register_subcycling_method_contract(const SubcyclingMethodContract &contract);
std::uint64_t
register_subcycling_group_schema(std::string_view provider_id,
                                 const SubcyclingGroupSchema &group_schema);

SubcyclingMethodContractSnapshot
require_registered_subcycling_method_contract();
SubcyclingMethodContractSnapshot
require_complete_registered_subcycling_method_contract();
void validate_registered_subcycling_method_contract(
    const SubcyclingMethodContractSnapshot &snapshot);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_METHOD_CONTRACT_HXX
