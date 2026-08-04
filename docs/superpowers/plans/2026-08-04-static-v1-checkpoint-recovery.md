# Static-v1 P2 checkpoint/recovery plan

## Goal

Prove one bounded, strict checkpoint/recovery path for the existing CarpetX
static two-level subcycling adapter without widening the geometry, schedule, or
science scope validated by P1.

## Global constraints

- Work only in the existing detached worktree
  `/home/cloud/ET/Cactus/scratch-subcycling-integration2-20260803`.
- Keep one patch, two levels, factor-two spatial/time refinement,
  `regrid_every=0`, CPU MPI 1 / OpenMP 1, and `TestODESolvers2` only.
- A checkpoint is legal only after a completed coarse/fine synchronization.
  Never serialize an in-flight dense interval or scratch transaction.
- Reconstruct the integration state from the strict checkpoint endpoint:
  root iteration `N`, root physical time `T`, exact clocks `{N,N}`, and
  accepted-step counts `{N,2N}`. The AMR topology epoch remains process-local.
- The recovery build and both halves of the smoke must use the same executable
  SHA-256, ODE method, tableau fingerprint, group schema, parameter physics,
  and grid layout.
- Reuse the existing local Silo/HDF5 installation for the checkpoint backend;
  do not rent a GPU or rebuild unrelated configurations.
- Follow TDD: every production change begins with a focused failing unit test.
- Preserve the P1 executable and artifacts unchanged.

## Task 1: Recovery seed contract and adapter integration

1. Add failing unit tests for a recovered synchronized seed at epoch 2 and
   for rejection of desynchronized, malformed, or non-finite recovery state.
2. Add a failing trace-equivalence unit test proving that an uninterrupted
   `N -> N+1` epoch and a newly reconstructed stepper produce the same clocks,
   accepted-step counts, physical time, and observer epoch.
3. Implement the smallest pure recovery-seed helper and use it in static-v1
   preflight/adapter construction. New starts retain the existing zero seed.
4. On strict recovery, canonicalize the freshly rebuilt CarpetX level clocks
   to the recovered root iteration only after all recoverable invariants pass.
5. Keep dynamic regrid, general thorns, mid-step recovery, and changed method
   contracts fail-closed.

## Task 2: Incremental recovery-capable build

1. Reuse the current p8d-native configuration and exact integration/flesh
   source links used by P1.
2. Add only the HDF5/Silo capabilities required by CarpetX checkpoint I/O,
   reusing the already installed local Silo tree.
3. Build incrementally, install a new immutable executable, record its hash,
   and restore canonical source links on every exit path.

## Task 3: One bounded native smoke

Use RKF78 and the P1 static-v1 smoke geometry with final synchronized time
`T=0.05`.

1. Run an uninterrupted no-checkpoint reference.
2. Run a seed case with Silo checkpoints at coarse epochs 2 and 4.
3. Isolate the epoch-2 checkpoint, then perform strict recovery from epoch 2
   to epoch 4 using the same executable.
4. Require recovery receipt `(iteration=2, time=0.025)`, resumed epochs 3 and
   4 exactly once, final `(iteration=4, time=0.05)`, `Done.`, and no hard fatal
   or stale dependent-state error.
5. Compare both levels' final `TestODESolvers2::state` against the
   uninterrupted reference. Require identical topology and values within a
   scale-aware roundoff tolerance; retain full raw TSVs and logs.

## Task 4: Evidence and handoff

1. Freeze pars, build/run logs, checkpoint inventory, comparison output,
   executable/source hashes, patch, and `SHA256SUMS` under the canonical
   Windows analysis tree.
2. Record what is validated and explicitly exclude dynamic regrid,
   Cottonmouth/BBH, tracker/AH/QLM, reflection-z, MPI scaling, and GPU.
3. Obtain an independent spec/code/evidence review before declaring P2 done.

