# Provenance and reproducibility

## Research reference

- Repository: https://github.com/rdm375/transient-pipe-state-finder
- Branch: `feature/shared-trajectory-matrix-free`
- Commit: `31b04d7c3c7446b579767a97997afef7a1f18b36`
- Referenced documents: `docs/M27_PHASE3_CLOSURE.md`, `docs/M28_OPERATIONAL_FORECAST_REPORT.md`, `docs/CHARACTERISTIC_TIMESCALES.md`, `docs/ENGINEERING_CASE_STUDY.md`.
- This archive's `src/forward.cpp` is a new, minimal adapter based on the research `src/cpp/forward_model.cpp` interface, not a verbatim copy.

## Forward solver

- Expected sibling checkout: `../isothermal-pipe-sim`.
- CMake target: `pipe_sim_cpp_core`.
- **Revision: UNRESOLVED.** The forward solver's remote could not be retrieved in this environment.
- Record an exact commit and review its license before production or distribution.
- No forward solver source code is bundled.

## Verification status

- Python config unit tests: runnable without external dependencies.
- C++ smoke test: provided but **not executed here** because the solver checkout is missing.
- No gradient, optimizer, uncertainty, or forecast acceptance tests have been run.
- Configuration schema is E1-only and does not yet specify measured time series.

## Reproducible run checklist

1. Checkout a known-good forward solver revision.
2. Record `git -C /path/to/isothermal-pipe-sim rev-parse HEAD`.
3. Configure and build this project against that checkout.
4. Run `ctest --test-dir build --output-on-failure`.
5. Run Python configuration tests.
6. Record compiler, CMake, platform, and numerical test outputs.
