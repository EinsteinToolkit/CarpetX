#include "subcycling_scratch_adapter_internal.hxx"

#include <AMReX_MultiFab.H>
#include <mpi.h>

#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace CarpetX;
using namespace CarpetX::subcycling_detail;

struct Source {
  amrex::MultiFab multifab;
  std::vector<ScratchValidity> validity;
  int token;
  Source(const int rank, const bool ok = true)
      : multifab(2, true, ok, 100 + rank),
        validity{{true, false, true}, {false, true, false}},
        token(100 + rank) {}
};

struct RecordingMpiCollective final : CollectiveOps {
  MpiCollectiveOps mpi;
  std::vector<std::int64_t> events;
  explicit RecordingMpiCollective(MPI_Comm &comm) : mpi(&comm) {}
  std::int64_t reduce_min(const ScratchCollectivePhase phase,
                          const std::size_t ordinal,
                          const std::int64_t local) override {
    events.push_back(100000 + static_cast<std::int64_t>(phase) * 1000 +
                     static_cast<std::int64_t>(ordinal));
    return mpi.reduce_min(phase, ordinal, local);
  }
  std::int64_t reduce_max(const ScratchCollectivePhase phase,
                          const std::size_t ordinal,
                          const std::int64_t local) override {
    events.push_back(200000 + static_cast<std::int64_t>(phase) * 1000 +
                     static_cast<std::int64_t>(ordinal));
    return mpi.reduce_max(phase, ordinal, local);
  }
  bool reduce_and(const ScratchCollectivePhase phase,
                  const std::size_t ordinal, const bool local) override {
    events.push_back(300000 + static_cast<std::int64_t>(phase) * 1000 +
                     static_cast<std::int64_t>(ordinal));
    return mpi.reduce_and(phase, ordinal, local);
  }
};

LocalScratchEntry entry(Source &source, const int group,
                        std::vector<std::int64_t> fields = {}) {
  if (fields.empty())
    fields = {83, 0, 1, group, 2, 2, 1, 0, 0, 4, 4, 4,
              7, 7, 7, 2, 0, 1, 1, 8, 8, 8, 0, 1, 1, 0, 1};
  const auto field_count = static_cast<std::int64_t>(fields.size());
  return LocalScratchEntry{ScratchGFKey{83, 1, 0, group}, &source.token,
                           &source.multifab, source.validity,
                           std::move(fields), field_count};
}

LocalScratchBatch batch(Source &source, const int group = 1) {
  return LocalScratchBatch{ScratchAdapterFailure::none, 83, 0, 1,
                           {entry(source, group)}, 1};
}

bool transact(LocalScratchBatch &local, RecordingMpiCollective &collective,
              const std::int64_t epoch, const bool fail_copy) {
  try {
    static_cast<void>(run_certification_transaction(
        local, collective, [epoch] { return epoch; },
        [fail_copy](const UncertifiedScratchLevelManifest &manifest,
                    const std::vector<ScratchGFView> &views) {
          if (fail_copy)
            throw std::runtime_error("rank-local copy failure");
          return ScratchLevelFrame::copy_tl0(manifest, views);
        }));
    return true;
  } catch (const ScratchAdapterCollectiveError &) {
    return false;
  }
}

void require_all(MPI_Comm comm, const bool local, const bool expected) {
  int value = local ? 1 : 0;
  int minimum = 0;
  int maximum = 0;
  MPI_Allreduce(&value, &minimum, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce(&value, &maximum, 1, MPI_INT, MPI_MAX, comm);
  assert(minimum == maximum);
  assert((minimum != 0) == expected);
}

void require_same_event_trace(
    MPI_Comm comm, const std::vector<std::int64_t> &local_events) {
  const int local_count = static_cast<int>(local_events.size());
  int minimum_count = 0;
  int maximum_count = 0;
  MPI_Allreduce(&local_count, &minimum_count, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce(&local_count, &maximum_count, 1, MPI_INT, MPI_MAX, comm);
  assert(minimum_count == maximum_count);
  std::vector<std::int64_t> peer(local_events.size());
  const int rank = [&] {
    int value = -1;
    MPI_Comm_rank(comm, &value);
    return value;
  }();
  MPI_Sendrecv(local_events.data(), local_count, MPI_INT64_T, 1 - rank, 9,
               peer.data(), local_count, MPI_INT64_T, 1 - rank, 9, comm,
               MPI_STATUS_IGNORE);
  assert(peer == local_events);
}

void require_outcome_and_trace(MPI_Comm comm, const bool local_result,
                               const bool expected,
                               const RecordingMpiCollective &collective) {
  require_all(comm, local_result, expected);
  require_same_event_trace(comm, collective.events);
}

void test_agreed_success(MPI_Comm comm, const int rank) {
  Source source(rank);
  auto local = batch(source);
  RecordingMpiCollective collective(comm);
  require_outcome_and_trace(comm, transact(local, collective, 83, false),
                            true, collective);
}

void test_one_rank_preflight_failure(MPI_Comm comm, const int rank) {
  Source source(rank);
  auto local = batch(source);
  if (rank == 1)
    local.preflight_status = ScratchAdapterFailure::invalid_tl0;
  RecordingMpiCollective collective(comm);
  require_outcome_and_trace(comm, transact(local, collective, 83, false),
                            false, collective);

  auto invalid_size = batch(source);
  if (rank == 1)
    invalid_size.entries[0].prevalidated_schema_length = -1;
  RecordingMpiCollective size_collective(comm);
  require_outcome_and_trace(
      comm, transact(invalid_size, size_collective, 83, false), false,
      size_collective);
}

void test_one_rank_count_schema_ok_epoch_copy_failures(MPI_Comm comm,
                                                       const int rank) {
  {
    Source source(rank), extra(rank + 10);
    auto local = batch(source);
    if (rank == 1) {
      local.entries.push_back(entry(extra, 2));
      local.prevalidated_entry_count = 2;
    }
    RecordingMpiCollective collective(comm);
    require_outcome_and_trace(comm, transact(local, collective, 83, false),
                              false, collective);
  }
  {
    Source source(rank);
    auto local = batch(source);
    if (rank == 1) {
      local.entries[0].schema_fields.push_back(101);
      ++local.entries[0].prevalidated_schema_length;
    }
    RecordingMpiCollective collective(comm);
    require_outcome_and_trace(comm, transact(local, collective, 83, false),
                              false, collective);
  }
  {
    Source source(rank);
    auto local = batch(source);
    if (rank == 1)
      local.entries[0].schema_fields[4] = 99;
    RecordingMpiCollective collective(comm);
    require_outcome_and_trace(comm, transact(local, collective, 83, false),
                              false, collective);
  }
  {
    Source source(rank, rank == 0);
    auto local = batch(source);
    RecordingMpiCollective collective(comm);
    require_outcome_and_trace(comm, transact(local, collective, 83, false),
                              false, collective);
  }
  {
    Source source(rank);
    auto local = batch(source);
    RecordingMpiCollective collective(comm);
    require_outcome_and_trace(
        comm, transact(local, collective, rank == 0 ? 83 : 84, false), false,
        collective);
  }
  {
    Source source(rank);
    auto local = batch(source);
    RecordingMpiCollective collective(comm);
    require_outcome_and_trace(
        comm, transact(local, collective, 83, rank == 1), false, collective);
  }
}

void test_one_rank_box_processor_validity_mismatch(MPI_Comm comm,
                                                   const int rank) {
  for (const std::size_t ordinal : {9U, 20U, 25U}) {
    Source source(rank);
    auto local = batch(source);
    if (rank == 1)
      local.entries[0].schema_fields[ordinal] += 1;
    RecordingMpiCollective collective(comm);
    require_outcome_and_trace(comm, transact(local, collective, 83, false),
                              false, collective);
  }
}

void test_collective_event_order(MPI_Comm comm, const int rank) {
  Source source(rank);
  auto local = batch(source);
  RecordingMpiCollective collective(comm);
  require_outcome_and_trace(comm, transact(local, collective, 83, false),
                            true, collective);
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = -1;
  int size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  assert(size == 2);
  MPI_Comm comm = MPI_COMM_WORLD;
  test_agreed_success(comm, rank);
  test_one_rank_preflight_failure(comm, rank);
  test_one_rank_count_schema_ok_epoch_copy_failures(comm, rank);
  test_one_rank_box_processor_validity_mismatch(comm, rank);
  test_collective_event_order(comm, rank);
  if (rank == 0)
    std::cout << "Phase 8B1 two-rank MPI protocol tests passed\n";
  MPI_Finalize();
}
