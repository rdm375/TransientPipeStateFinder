# Engineering requirements

> The purpose of our research is not simply to reconstruct a pipeline state. It is to determine how much measurement information is needed, how rapidly uncertainty loses its influence on future behavior, and whether the resulting predictions become sufficiently reliable in time to support an engineering decision.

## Requirements

1. Use the validated C++ forward solver directly; never implement an alternate physics kernel.
2. Reconstruct from measured boundary histories with explicit units, noise assumptions, and input roles.
3. Verify adjoint and tangent-linear gradients against independent finite differences.
4. Preserve discrete global mass conservation; distinguish physical inlet–outlet imbalance from numerical conservation defect.
5. Report optimizer convergence and initialization sensitivity, not just objective reduction.
6. Distinguish structural rank, practical identifiability, and predictive uncertainty.
7. Forecast operational pressure, inventory, and threshold crossings.
8. Report characteristic times separately, including forecast settling and operational warning intervals.
9. Separate retrospective validation from prospective uncertainty guarantees.
10. Record solver/research commit hashes, configuration, numerical tolerances, and data provenance.

## E1 acceptance

- C++20 CMake project with direct forward solver linkage.
- Validated, versioned JSON configuration schema and sample.
- Smoke test exercises the real solver (not a stub).
- Transfer manifest and provenance/limitations.
- Build and tests must run in an environment with the forward solver installed.

## Current status

The E1 foundation has progressed through the E5.7b research
checkpoint. The repository now contains constrained reconstruction,
measurement-information analysis, exact tangent-linear forecast
sensitivities, local uncertainty propagation, and engineering
decision classification.

At commit `ae83939`, the E5.7b experiment and the complete
46-test regression suite passed.

The current evidence is based on synthetic single-pipe experiments.
Large-scale solver performance, nonlinear uncertainty guarantees,
measurement-design optimization, and the engineering report
generator remain open work.

The forward solver is linked directly from the external
`isothermal-pipe-sim` source tree. Its exact revision and
licensing must be recorded for reproducible distribution.

See `docs/RESEARCH_ROADMAP.md`.
