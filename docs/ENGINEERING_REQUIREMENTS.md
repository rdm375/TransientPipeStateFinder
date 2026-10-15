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

Foundation files and schema unit tests are included. C++ integration has not been compiled in this packaging environment because the forward solver checkout is unavailable. No estimator or report generator is claimed to exist.
