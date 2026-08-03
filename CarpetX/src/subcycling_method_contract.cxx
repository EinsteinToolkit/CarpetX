#include "subcycling_method_contract.hxx"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace CarpetX {
namespace {

bool known_method(const SubcyclingODEMethod method) noexcept {
  switch (method) {
  case SubcyclingODEMethod::rk4:
  case SubcyclingODEMethod::rkf78_order7:
  case SubcyclingODEMethod::dp87_order8:
    return true;
  }
  return false;
}

bool zero_fingerprint(const TableauFingerprint &fingerprint) noexcept {
  return std::all_of(fingerprint.begin(), fingerprint.end(),
                     [](const std::uint8_t value) { return value == 0; });
}

bool same_dense(const DenseCapability &left,
                const DenseCapability &right) noexcept {
  return left.method == right.method &&
         left.tableau_fingerprint == right.tableau_fingerprint &&
         left.endpoint_order == right.endpoint_order &&
         left.dense_uniform_order == right.dense_uniform_order &&
         left.stage_count == right.stage_count &&
         left.extra_rhs_evaluations == right.extra_rhs_evaluations &&
         left.persistent_vector_count == right.persistent_vector_count &&
         left.arbitrary_theta == right.arbitrary_theta &&
         left.verified == right.verified;
}

bool same_contract(const SubcyclingMethodContract &left,
                   const SubcyclingMethodContract &right) noexcept {
  return left.provider_id == right.provider_id &&
         left.parameter_value == right.parameter_value &&
         same_dense(left.dense, right.dense);
}

bool same_schema(const SubcyclingGroupSchema &left,
                 const SubcyclingGroupSchema &right) noexcept {
  if (left.ordered_group_pairs.size() !=
          right.ordered_group_pairs.size() ||
      left.dependent_groups != right.dependent_groups)
    return false;
  for (std::size_t index = 0; index < left.ordered_group_pairs.size();
       ++index) {
    const auto &left_pair = left.ordered_group_pairs[index];
    const auto &right_pair = right.ordered_group_pairs[index];
    if (left_pair.evolved_group != right_pair.evolved_group ||
        left_pair.rhs_group != right_pair.rhs_group)
      return false;
  }
  return true;
}

void validate_contract(const SubcyclingMethodContract &contract) {
  const auto &dense = contract.dense;
  if (contract.provider_id.empty() || contract.parameter_value.empty())
    throw std::invalid_argument(
        "subcycling method contract identity must not be empty");
  if (!known_method(dense.method) ||
      zero_fingerprint(dense.tableau_fingerprint) ||
      dense.endpoint_order <= 0 ||
      dense.dense_uniform_order < dense.endpoint_order ||
      dense.stage_count <= 0 || dense.extra_rhs_evaluations < 0 ||
      dense.persistent_vector_count <= 0 || !dense.arbitrary_theta ||
      !dense.verified)
    throw std::invalid_argument(
        "subcycling method contract has an invalid dense capability");
}

} // namespace

bool SubcyclingMethodContractRegistry::empty() const noexcept {
  return !snapshot_.has_value();
}

std::uint64_t SubcyclingMethodContractRegistry::register_contract(
    const SubcyclingMethodContract &contract) {
  validate_contract(contract);
  if (snapshot_) {
    if (snapshot_->contract.provider_id != contract.provider_id)
      throw std::logic_error(
          "another provider already owns the subcycling method contract");
    if (!same_contract(snapshot_->contract, contract))
      throw std::logic_error(
          "the subcycling method contract is frozen for this GHExt");
    return snapshot_->registration_epoch;
  }
  snapshot_.emplace(SubcyclingMethodContractSnapshot{
      contract, std::nullopt, std::uint64_t{1}});
  return snapshot_->registration_epoch;
}

std::uint64_t SubcyclingMethodContractRegistry::register_group_schema(
    const std::string_view provider_id,
    const SubcyclingGroupSchema &group_schema) {
  if (!snapshot_)
    throw std::logic_error(
        "cannot register a group schema before a method contract");
  if (provider_id != snapshot_->contract.provider_id)
    throw std::logic_error(
        "group schema provider does not own the method contract");
  if (snapshot_->group_schema) {
    if (!same_schema(*snapshot_->group_schema, group_schema))
      throw std::logic_error(
          "the subcycling group schema is frozen for this GHExt");
    return snapshot_->registration_epoch;
  }
  if (snapshot_->registration_epoch ==
      std::numeric_limits<std::uint64_t>::max())
    throw std::overflow_error(
        "subcycling method-contract registration epoch exhausted");
  snapshot_->group_schema = group_schema;
  ++snapshot_->registration_epoch;
  return snapshot_->registration_epoch;
}

SubcyclingMethodContractSnapshot
SubcyclingMethodContractRegistry::require_snapshot() const {
  if (!snapshot_)
    throw std::out_of_range("no subcycling method contract is registered");
  return *snapshot_;
}

SubcyclingMethodContractSnapshot
SubcyclingMethodContractRegistry::require_complete_snapshot() const {
  const auto snapshot = require_snapshot();
  if (!snapshot.group_schema)
    throw std::out_of_range(
        "the subcycling method contract has no immutable group schema");
  return snapshot;
}

void SubcyclingMethodContractRegistry::validate_snapshot(
    const SubcyclingMethodContractSnapshot &snapshot) const {
  if (!snapshot_ ||
      snapshot.registration_epoch != snapshot_->registration_epoch ||
      !same_contract(snapshot.contract, snapshot_->contract) ||
      snapshot.group_schema.has_value() !=
          snapshot_->group_schema.has_value() ||
      (snapshot.group_schema &&
       !same_schema(*snapshot.group_schema, *snapshot_->group_schema)))
    throw std::logic_error(
        "subcycling method-contract snapshot is stale or modified");
}

} // namespace CarpetX
