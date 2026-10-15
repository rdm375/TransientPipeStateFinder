# Scored Forward Trajectories

## Research objective

The purpose of our research is not simply to reconstruct a
pipeline state. It is to determine how much measurement
information is needed, how rapidly uncertainty loses its
influence on future behavior, and whether the resulting
predictions become sufficiently reliable in time to support
an engineering decision.

## State recording

The recording-capable scheduled-forward integrator preserves
the complete physical state at every successfully converged
timestep, including the initial state.

The state ordering is:

    [p_0, m_in, m_0, p_1, ..., m_(n-1), p_n]

Every snapshot contains:

- Timestep index.
- Simulation time.
- Complete physical state.
- Prescribed inlet pressure and outlet mass flow.
- Four physical boundary channels.
- Numerical quality metrics.
- Optional measurement and uncertainty metrics.

## Numerical quality

### Dynamics residual

For k >= 1:

    residual_linf = ||F_(k-1)||_infinity

The residual is assembled using the validated pipe solver.

The initial state's residual is unavailable because there
is no preceding timestep.

### Mass-conservation defect

For k >= 1:

    defect = abs(
        (M_k - M_(k-1))/dt
        - theta*(m_in,k - m_out,k)
        - (1-theta)*(m_in,k-1 - m_out,k-1)
    )

Units: kg/s.

The implementation preserves the validated integrator's
arithmetic grouping.

The initial conservation defect is unavailable.

### Physical mass imbalance

    imbalance = m_in - m_out

Units: kg/s.

This is a physical operating diagnostic, not a numerical
conservation error.

Nonzero imbalance can indicate changing linepack.

### Linepack

Linepack is calculated using pipe_sim::linepack_cz.

Units: kg.

## Measurement scoring

Measurement scores require actual observations,
availability masks, and measurement uncertainty scales.

Normalized residual:

    r = (predicted - observed)/sigma

A weighted RMS score may be computed over the available
measurements.

Missing measurements must not be interpreted as zero error.

Measurement scoring is performed by the reconstruction layer,
not by the forward integrator.

## State uncertainty

State uncertainty requires an uncertainty model, ensemble,
covariance approximation, or admissible state set.

Pressure and mass-flow errors require appropriate scaling.

Uncomputed uncertainty scores remain unavailable.

## Forecast reliability

Prediction reliability depends on:

- Quantity of interest.
- Prediction horizon.
- Future boundary conditions.
- State uncertainty.
- Acceptance criteria.

Forecast-specific reliability assessments are separate from
basic numerical quality scores.

## Restart readiness

A state is suitable for restarting a prediction only when
the required numerical, measurement, and predictive
criteria are satisfied.

The planned readiness classifications are:

    not_evaluated
    not_ready
    ready

A numerically converged state is not automatically a
reliable reconstructed state.

Readiness is specific to an engineering prediction.

## Implementation

The APIs are:

    run_forward_scheduled(...)
    run_forward_scheduled_recorded(...)

Both share one integration implementation.

The recorded API returns:

    RecordedForwardResult {
        Result summary;
        ScoredTrajectory trajectory;
    }

The existing non-recording API retains its behavior.

Successful integration produces steps + 1 snapshots.

Early termination produces completed_steps + 1 snapshots.

Only successfully converged states are recorded.

The implementation does not modify the validated pipe
equations, Newton iteration, or timestep discretization.

## Required verification

- Existing scheduled-forward regression tests pass.
- Initial snapshot equals the supplied initial state.
- Every converged timestep is recorded.
- Recorded final state matches the original forward result.
- Boundary histories agree with the original forward result.
- Numerical scores agree with independent calculations.
- Initial residual and conservation scores are unavailable.
- Measurement and uncertainty scores remain unavailable.
- Failed timesteps are not recorded as converged.
- Recorded trajectories satisfy full-space constraints.
- Debug and Release test suites pass.

## Future development

- Measurement scoring.
- State uncertainty propagation.
- Error damping and persistence analysis.
- Forecast uncertainty estimation.
- Restart-readiness classification.
- Engineering report generation.
- Persistent trajectory storage and reproducibility metadata.
