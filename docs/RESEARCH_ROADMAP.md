# TransientPipeStateFinder research roadmap

> The purpose of our research is not simply to reconstruct a pipeline state. It is to determine how much measurement information is needed, how rapidly uncertainty loses its influence on future behavior, and whether the resulting predictions become sufficiently reliable in time to support an engineering decision.

## Research boundaries

- Reuse the existing `isothermal-pipe-sim` C++ forward solver.
- Do not create a replacement physics solver.
- Keep reconstruction, uncertainty analysis, optimization,
  and reporting in TransientPipeStateFinder.
- Use AGPL-3.0-only for the statefinder project, subject to
  dependency and redistribution review.
- Preserve discrete mass conservation.
- Distinguish physical imbalance from numerical mass defect.
- Reserve theta for the time-integration weighting parameter.
- Reserve p/P for pressure, not generic parameter vectors.
- Use bold xi for generic scenario parameters.
- Keep future boundary conditions fixed in current
  forecast-sensitivity experiments.
- Distinguish local linearized results from nonlinear guarantees.

## Completed research checkpoints

| Milestone | Scope | Status |
|---|---|---|
| E1-E3 | Forward integration, reconstruction, continuation | Completed research foundation |
| E4 | Scalable SQP research infrastructure | Implemented; scalability remains open |
| E5.1 | Reconstruction quality and identifiability | Completed |
| E5.2 | Engineering-output sensitivity | Completed |
| E5.3 | Measurement uncertainty | Completed |
| E5.4 | Exact tangent-linear forecast attenuation | Completed |
| E5.5 | Decision reliability and information geometry | Completed |
| E5.6a | Persistent sensitivity attenuation | Completed |
| E5.6b | Operational decision timescales | Completed |
| E5.7a | Observation timing and information attenuation | Completed |
| E5.7b | Measurement sufficiency and critical limits | Completed |
| E5.7c | Decision sufficiency without complete identifiability | Next |

## E5.7b validated results

Checkpoint: `ae83939`, branch `feature/e4-scalable-sqp`.
Validation: 46/46 tests passed.

| Observation times [s] | Leading information singular value | Critical pressure limit [Pa] |
|---|---:|---:|
| None | 0 | No finite limit under unbounded unresolved uncertainty |
| 0 | 34.97744442462 | 4999519.425071 |
| 1 | 11.18497294641 | 4999218.957275 |
| 2 | 4.940619264993 | 4998660.676650 |
| 0, 1, 2 | 37.05313693888 | 4999527.337893 |

The critical pressure is the largest constant minimum-pressure
requirement supported over the sampled forecast by the local
linearized uncertainty model.

Observed configurations passed independent classifier checks
at one pascal below and above the calculated boundary.

No-measurement reconstruction converged but retained an
approximately 20 kPa initial-pressure error and an unresolved
forecast-sensitive direction.

These are synthetic results, not universal pipeline properties.

## E5.7c — Incomplete identifiability

See `E57C_DECISION_SUFFICIENCY.md`.

Goal: establish whether a decision can be reliable when the
state remains incompletely identified.

## E5.8 — Measurement design

Proposed research:
- Measurement placement and timing.
- Measurement precision and cost.
- Minimum measurement sets supporting specified decisions.
- Decision-oriented versus state-oriented information criteria.

## E5.9 — Nonlinear robustness

Proposed research:
- Finite perturbation validation.
- Measurement noise and model mismatch.
- Sensitivity to rank thresholds and scaling.
- Empirical uncertainty-envelope coverage.
- Uncertain future boundary conditions.
- Time discretization and forecast sampling sensitivity.

## E6 — Engineering reporting

Planned deliverable:
- Reproducible model and measurement provenance.
- Reconstruction and convergence diagnostics.
- Identifiability and unresolved-direction analysis.
- Forecast and conservation diagnostics.
- Decision classification and characteristic times.
- Uncertainty assumptions and validation evidence.
- Explicit limitations and engineering interpretation.

## Remaining engineering gates

1. Pin the external forward-solver revision.
2. Verify dependency licensing and distribution requirements.
3. Establish scalability on representative larger problems.
4. Establish nonlinear robustness and uncertainty coverage.
5. Implement and validate engineering report generation.
