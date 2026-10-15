# TransientPipeStateFinder

**Research software for reconstructing transient pipeline states and determining when those states are sufficiently accurate to initialize subsequent engineering simulations.**

> The purpose of our research is not simply to reconstruct a pipeline state. It is to determine how much measurement information is needed, how rapidly uncertainty loses its influence on future behavior, and whether the resulting predictions become sufficiently reliable in time to support an engineering decision.

## Primary Research Objective: State Qualification Time

The central research question is:

**At what point during a statefinding analysis does the reconstructed state become sufficiently reliable to serve as the initial condition for a subsequent transient simulation?**

We call this the **State Qualification Time**, denoted \(t_Q\).

At each candidate analysis time, the statefinder uses only measurements available through that time. The reconstructed state then initializes an independent nonlinear forward simulation.

The resulting pressure and mass-flow errors are compared against a synthetic reference over a specified restart horizon.

A state qualifies when both errors satisfy their prescribed tolerances.

The qualification time is the earliest tested analysis time for which every subsequently tested candidate also qualifies.

## Current Research Status

**Current milestone: R5.1 — Causal State Qualification and Nonlinear Restart Validation.**

Implemented capabilities include:

- Direct reuse of the validated isothermal transient-pipe forward solver.
- Constrained full-space state reconstruction.
- Four-channel boundary measurement assimilation.
- Dense nonlinear SQP optimization.
- Exact discrete tangent-linear and adjoint operations.
- Local identifiability and measurement-information analysis.
- Forecast continuation.
- Initial-state sensitivity and attenuation analysis.
- Engineering decision classification.
- Engineering report generation.
- Causal state qualification with independent nonlinear restarts.

### R5.1 Results

| Quantity | Result |
|---|---:|
| Candidate analysis times | 1–8 seconds |
| Restart horizon | 4 seconds |
| Pressure tolerance | 10,000 Pa |
| Mass-flow tolerance | 0.1 kg/s |
| Reconstruction convergence | 8/8 |
| Qualified candidates | 8/8 |
| Earliest tested qualifying time | 1 second |
| Maximum pressure restart error | Approximately 9.31e-10 Pa |
| Maximum flow restart error | Approximately 4.03e-13 kg/s |

**Important:** The one-second result is the earliest tested qualifying time, not an established physical statefinder spin-up time.

The experiment is noise-free, uses identical forward and inverse models, and fixes most initial-state coordinates.

It validates the numerical restart-qualification infrastructure but does not yet demonstrate the dissipation of meaningful initial-state deficiencies.

## Next Research Milestone: R5.2

**Initial-State Deficiency and Qualification-Time Study**

The immediate objective is to determine how long meaningful deficiencies in reconstructed pipeline states persist before those states become suitable for subsequent simulation.

Research tasks:

1. Release additional physically admissible initial-state unknowns.
2. Introduce controlled pressure and mass-flow errors.
3. Introduce reproducible measurement noise.
4. Perform causal reconstructions at successive analysis times.
5. Restart the nonlinear forward solver from each reconstructed state.
6. Compare restarted trajectories against synthetic truth.
7. Determine persistent qualification times.
8. Investigate scenarios where qualification is never attained.

The principal diagnostic is the normalized restart error \(Q(t,H)\), with qualification threshold \(Q=1\).

## Build and Test

Requirements:

- CMake 3.20 or newer.
- C++20 compiler.
- Eigen3.
- Compatible checkout of `isothermal-pipe-sim`.

Configure and build:

    cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DPIPE_SIM_SOURCE_DIR="$HOME/Sync/isothermal-pipe-sim-statefinder"

    cmake --build build -j"$(nproc)"

Run tests:

    ctest --test-dir build --output-on-failure

Run the engineering demonstration:

    ./build/tpsf_engineering_demo reports/engineering-demo

## Engineering Reports

The engineering report must ultimately lead with the state qualification assessment.

Current generated reports include:

- `state_qualification.json`
- `state_qualification.md`
- `engineering_assessment.json`
- `engineering_assessment.md`
- `decision_boundary_sweep.json`
- `decision_boundary_sweep.md`

State qualification is currently reported separately. Integration into the unified engineering report remains outstanding.

## Scientific Limitations

Current results are based on small synthetic single-pipe experiments.

The project has not established:

- A universal statefinder spin-up time.
- Qualification under arbitrary future boundary conditions.
- Operational accuracy using field measurements.
- Nonlinear uncertainty guarantees.
- Production scalability.
- Qualification without a trustworthy reference state.

Numerical optimizer convergence is not equivalent to state qualification.

Sensitivity attenuation is not equivalent to state qualification.

The principal research objective is to determine when reconstructed states produce sufficiently accurate subsequent simulations.

## License

Intended license: **AGPL-3.0-only**.

External forward-solver licensing and redistribution requirements must be verified before distribution.

See [docs/LICENSING.md](docs/LICENSING.md).

---

**Central research question:** At what time during an analysis does an initially deficient reconstructed pipeline state become sufficiently reliable to initialize subsequent engineering simulations, and what evidence supports that qualification?
