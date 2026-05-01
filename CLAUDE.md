# CarpetX

CarpetX is a Cactus driver based on AMReX for block-structured adaptive mesh refinement (AMR), targeting the Einstein Toolkit for numerical relativity simulations. It is written in C++ and supports both CPU and GPU accelerators.

## Project Structure & Module Organization

This repository is organized as Cactus thorns. Each top-level directory such as `CarpetX/`, `ODESolvers/`, `TestSubcyclingMC/`, or `TestSubcyclingMC2/` is a thorn with Cactus metadata files like `interface.ccl`, `param.ccl`, `schedule.ccl`, and implementation under `src/`. Tests usually live beside the thorn in `test/` and use parameter files plus checked-in reference output, for example `TestSubcyclingMC2/test/gaussian-interprocessonly.par` and `TestSubcyclingMC2/test/gaussian-interprocessonly/*.tsv`.

Do **not** read reference test data files (e.g. `*/test/*/*.tsv`, `*/test/*/*.h5`, or other large numerical output) unless the user explicitly asks for it. They are large, binary-like, and rarely informative for code-level reasoning. Inspect the `.par` parameter file and the thorn source instead.

## Build & Test

```bash
./agent_scripts/build.sh
./agent_scripts/test.sh
```

Both scripts run the build/tests inside a Docker container — the host
checkout is bind-mounted, but build artifacts (object files, configs,
testsuite logs) live **inside the container**, not on the host. To
inspect them (e.g. `nm` on a `.o`, tailing a testsuite log), exec into
the same container the scripts use:

```bash
docker exec "$CONTAINERLOCAL" zsh -c '
  find "$CONTAINERLOCALCACTUS/configs/sim-carpetx" -name "<file>"
'
```

The relevant env vars (`CONTAINERLOCAL`, `CONTAINERLOCALCACTUS`,
`CONTAINERLOCALMACHINE`) are already set in the host shell that runs
`agent_scripts/*.sh`.

## Commit & Pull Request Guidelines

Recent history favors short subject lines with a thorn prefix, for example `ODESolvers: update comments` or `CarpetX: add par poison_undefined_arrays`. Keep commits focused on one thorn or one behavior change. PRs should explain the scientific or runtime impact, list the tests run, and link the relevant issue or discussion. Include plots or output snippets when a change affects diagnostics, IO, or visible simulation results.

## Further Reading

- `.github/workflows/ci.yml` — CI matrix and environment variable combinations
