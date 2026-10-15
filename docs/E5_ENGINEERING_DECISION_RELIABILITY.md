# E5 — Engineering Decision Reliability

## Research objective

The purpose of our research is not simply to reconstruct a pipeline
state. It is to determine how much measurement information is needed,
how rapidly uncertainty loses its influence on future behavior, and
whether the resulting predictions become sufficiently reliable in
time to support an engineering decision.

## Architecture decision

Use the existing full-space SQP implementation as the initial
statefinding engine.

Preserve the full-space trajectory representation, exact discrete
constraints, tangent/adjoint operations, dense reference SQP,
matrix-free KKT infrastructure, and forecast continuation.

Scalable matrix-free nonlinear SQP remains valuable research,
but it is not a prerequisite for E5/E6 engineering reports.

## E5.1 — Statefinding quality

At the reconstructed scaled coordinates u, assemble

    A = dc/du
    B = dr/du

where c is the scaled discrete constraint residual and r is the
normalized available-measurement residual.

Find an orthonormal numerical nullspace basis Z:

    A Z = 0

The locally feasible measurement sensitivity is

    S = B Z

Report:

- SQP termination status
- Nonlinear constraint residual norm
- Normalized measurement residual norm and objective
- Projected gradient norm
- Constraint rank
- Dimension of the feasible tangent space
- Rank and singular values of S
- Number and basis of unresolved feasible directions

The singular-value threshold is a numerical convention and must
be reported. Rank depends on coordinate scaling and threshold.

The assessment is LOCAL and LINEARIZED. It is not proof of global
uniqueness, convergence from arbitrary initialization, or accuracy
of the reconstructed state.

An unobservable direction does not automatically imply uncertainty
in every engineering output. Output-specific sensitivity is required.

A fully observed tangent space does not automatically imply small
uncertainty. Weak singular directions and measurement-noise
assumptions remain important.

No statistical confidence bounds are inferred from rank alone.

## E5.2 — Engineering-output sensitivity

Evaluate linepack, minimum pressure, and relevant boundary
quantities along feasible and weakly observed directions.

Use physical-unit gradients with the scaled-coordinate mapping.

Distinguish physical inlet/outlet imbalance from numerical
mass-conservation defects.

## E5.3 — Uncertainty models

Introduce explicit assumptions for measurement noise, initial
uncertainty, model discrepancy, and uncertain boundary conditions.

Do not present inverse singular values as confidence intervals
without a defensible probabilistic or deterministic uncertainty
model.

## E5.4 — Forecast propagation

Propagate admissible reconstruction perturbations through the
existing validated forward solver.

Assess forecast output sensitivity as a function of horizon.

Do not assume that all state errors decay, or that a single
exponential attenuation rate describes the dynamics.

## E5.5 — Decision reliability

Compare predicted engineering quantities and defensible uncertainty
envelopes with operational thresholds.

Allow the report to conclude:

- supported under stated assumptions
- not supported
- insufficient information

## E6 — Engineering report

Report numerical convergence, observability, reconstruction
sensitivity, forecast sensitivity, conservation diagnostics,
operational thresholds, assumptions, and limitations.

A successful optimizer termination is not itself an engineering
decision certificate.

## Initial test

A one-time-level two-cell fixture has one feasible pressure
direction. An outlet-pressure measurement resolves it; removing
the measurement leaves one unresolved direction.

This tests local information structure, not nonlinear SQP recovery.
The existing E2 reconstruction tests remain the nonlinear baseline.
