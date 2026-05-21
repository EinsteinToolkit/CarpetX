# CarpetX

CarpetX is a Cactus driver based on AMReX for block-structured adaptive mesh refinement (AMR), targeting the Einstein Toolkit for numerical relativity simulations. It is written in C++ and supports both CPU and GPU accelerators.

## Project Structure & Module Organization

This repository is organized as Cactus thorns. Each top-level directory such as `CarpetX/`, `ODESolvers/`, `TestSubcyclingMC/`, or `TestSubcyclingMC2/` is a thorn with Cactus metadata files like `interface.ccl`, `param.ccl`, `schedule.ccl`, and implementation under `src/`. Tests usually live beside the thorn in `test/` and use parameter files plus checked-in reference output, for example `TestSubcyclingMC2/test/gaussian.par` and `TestSubcyclingMC2/test/gaussian/*.tsv`.

Do **not** read reference test data files (e.g. `*/test/*/*.tsv`, `*/test/*/*.h5`, or other large numerical output) unless the user explicitly asks for it. They are large, binary-like, and rarely informative for code-level reasoning. Inspect the `.par` parameter file and the thorn source instead.

## Build & Test

```bash
./agent_scripts/build.sh
./agent_scripts/test.sh
```

Both scripts run inside the running container named `$CONTAINERLOCAL`; build artifacts and testsuite logs live **inside the container**, not on the host. `$CONTAINERLOCAL` and `$CONTAINERLOCALCACTUS` are pre-set in the host shell — pass `CONTAINERLOCALCACTUS` through with `-e` (it's not in the container's default env, and a bare `"$CONTAINERLOCALCACTUS"` inside the `zsh -c` string expands to empty):

```bash
docker exec -e CONTAINERLOCALCACTUS="$CONTAINERLOCALCACTUS" \
  "$CONTAINERLOCAL" zsh -c 'tail -80 "$CONTAINERLOCALCACTUS/<path>"'
```

Per-test logs: `$CONTAINERLOCALCACTUS/TEST/sim-carpetx/<Thorn>/<test>.log`.
Build artifacts: `$CONTAINERLOCALCACTUS/configs/sim-carpetx/`.

## Commit & Pull Request Guidelines

Recent history favors short subject lines with a thorn prefix, for example `ODESolvers: update comments` or `CarpetX: add par poison_undefined_arrays`. Keep commits focused on one thorn or one behavior change. PRs should explain the scientific or runtime impact, list the tests run, and link the relevant issue or discussion. Include plots or output snippets when a change affects diagnostics, IO, or visible simulation results.

## Further Reading

- `.github/workflows/ci.yml` — CI matrix and environment variable combinations
