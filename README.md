# TransientPipeStateFinder

Research implementation for single-pipe isothermal transient state estimation, uncertainty propagation, and engineering decision reliability.

> The purpose of our research is not simply to reconstruct a pipeline state. It is to determine how much measurement information is needed, how rapidly uncertainty loses its influence on future behavior, and whether the resulting predictions become sufficiently reliable in time to support an engineering decision.

**Status:** E5.7b research checkpoint completed and validated.
The project implements constrained full-space reconstruction,
local identifiability analysis, exact discrete tangent-linear
forecast sensitivities, uncertainty propagation, and sampled-time
engineering decision classification.

The E5.7b experiment passed 46/46 regression tests at commit
`ae83939` on branch `feature/e4-scalable-sqp`.

These are research results for synthetic single-pipe cases.
They are not nonlinear uncertainty certificates or claims of
production scalability or operational certification.

The engineering report generator remains planned.

See [Research roadmap](docs/RESEARCH_ROADMAP.md) and
[E5.7c experimental specification](docs/E57C_DECISION_SUFFICIENCY.md).

## Quick start

Place `isothermal-pipe-sim` alongside this directory, or provide its absolute path:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPIPE_SIM_SOURCE_DIR="$HOME/Sync/isothermal-pipe-sim"
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

If the solver is absent, configuration fails clearly; no mock physics is substituted.

Standalone configuration validation (does not require C++):

```bash
python3 -m unittest discover -s tests/python -v
python3 scripts/validate_config.py examples/single_pipe.json
```

## Sources and boundaries

- Forward solver: external `isothermal-pipe-sim` C++ target `pipe_sim_cpp_core`. **Exact solver commit not yet resolved** (repository was unavailable through the GitHub connection). Pin before declaring a reproducible physics build.
- Research: [rdm375/transient-pipe-state-finder](https://github.com/rdm375/transient-pipe-state-finder), branch `feature/shared-trajectory-matrix-free`, commit `31b04d7c3c7446b579767a97997afef7a1f18b36`.
- The research repository is **not** a runtime dependency.
- See `docs/TRANSFER_MANIFEST.md` and `docs/ENGINEERING_REQUIREMENTS.md`.

## License

Project intent: `AGPL-3.0-only`. Review forward-solver license and any future transferred files before distribution. See `LICENSE-SPDX` and `docs/LICENSING.md`.
