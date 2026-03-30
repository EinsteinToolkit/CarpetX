# CarpetX

CarpetX is a Cactus driver based on AMReX for block-structured adaptive mesh refinement (AMR), targeting the Einstein Toolkit for numerical relativity simulations. It is written in C++ and supports both CPU and GPU accelerators.

## Project Structure & Module Organization

This repository is organized as Cactus thorns. Each top-level directory such as `CarpetX/`, `ODESolvers/`, `WaveToyX/`, or `TestSubcycling/` is a thorn with Cactus metadata files like `interface.ccl`, `param.ccl`, `schedule.ccl`, and implementation under `src/`. Tests usually live beside the thorn in `test/` and use parameter files plus checked-in reference output, for example `WaveToyX/test/presync.par` and `WaveToyX/test/presync/*.tsv`.

## Build & Test

```bash
./agent_scripts/build.sh
./agent_scripts/test.sh
```

## Commit & Pull Request Guidelines

Recent history favors short subject lines with a thorn prefix, for example `ODESolvers: update comments` or `CarpetX: add par poison_undefined_arrays`. Keep commits focused on one thorn or one behavior change. PRs should explain the scientific or runtime impact, list the tests run, and link the relevant issue or discussion. Include plots or output snippets when a change affects diagnostics, IO, or visible simulation results.

## Further Reading

- `README.md` — project overview, links to talks and getting-started wiki
- `CarpetX/doc/` — architecture and design documentation
- `.github/workflows/ci.yml` — CI matrix and environment variable combinations
