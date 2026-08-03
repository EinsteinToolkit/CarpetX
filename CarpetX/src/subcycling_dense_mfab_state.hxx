#ifndef CARPETX_SUBCYCLING_DENSE_MFAB_STATE_HXX
#define CARPETX_SUBCYCLING_DENSE_MFAB_STATE_HXX

#include "subcycling_dense_output.hxx"

#include <AMReX_MultiFab.H>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace CarpetX {

struct DenseMFabKey {
  std::int64_t hierarchy_epoch;
  int patch;
  int level;
  int group_index;
};

struct DenseMFabView {
  DenseMFabKey key;
  const amrex::MultiFab *multifab;
};

class OwnedMultiFabDenseState final : public DenseStateVector {
public:
  ~OwnedMultiFabDenseState() override = default;

  OwnedMultiFabDenseState(const OwnedMultiFabDenseState &) = delete;
  OwnedMultiFabDenseState &operator=(const OwnedMultiFabDenseState &) = delete;
  OwnedMultiFabDenseState(OwnedMultiFabDenseState &&) = delete;
  OwnedMultiFabDenseState &operator=(OwnedMultiFabDenseState &&) = delete;

  static std::unique_ptr<OwnedMultiFabDenseState>
  copy_of(const std::vector<DenseMFabView> &views);

  static std::unique_ptr<OwnedMultiFabDenseState>
  allocate_like(const std::vector<DenseMFabView> &views);

  static std::unique_ptr<OwnedMultiFabDenseState>
  empty_like(const OwnedMultiFabDenseState &source);

  static std::unique_ptr<OwnedMultiFabDenseState> linear_combination_of(
      const std::vector<double> &weights,
      const std::vector<const OwnedMultiFabDenseState *> &sources);

  std::size_t entry_count() const noexcept;
  const DenseMFabKey &key(std::size_t index) const;
  const amrex::MultiFab &multifab(std::size_t index) const;

  bool compatible(const DenseStateVector &other) const noexcept override;
  void copy_from(const DenseStateVector &other) override;
  void linear_combination(
      const std::vector<double> &weights,
      const std::vector<const DenseStateVector *> &sources) override;

private:
  struct Entry {
    DenseMFabKey key;
    std::unique_ptr<amrex::MultiFab> multifab;
  };

  explicit OwnedMultiFabDenseState(std::vector<Entry> entries);

  std::vector<Entry> entries_;
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_DENSE_MFAB_STATE_HXX
