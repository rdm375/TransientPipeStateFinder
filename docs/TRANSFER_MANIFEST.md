# Nine-component transfer manifest

Research source: `rdm375/transient-pipe-state-finder`, `feature/shared-trajectory-matrix-free`,
commit `31b04d7c3c7446b579767a97997afef7a1f18b36`.

**Classification:** `E1 integrated` means this archive includes a thin adapter.
`Planned` means only the research provenance and intended interface are documented;
the numerical algorithm has **not** been copied or implemented.

| # | Component | Research source examples | E1 status | Required gate |
|---|---|---|---|---|
| 1 | Forward-solver integration | `src/cpp/forward_model.cpp`, `include/state_finder/forward_model.hpp`, `tests/cpp/forward_smoke.cpp` | E1 integrated, unbuilt here | Pin solver SHA; compile and run CTest |
| 2 | Discrete adjoint | `src/cpp/adjoint.cpp`, `src/cpp/old_state_vjp.cpp`, `docs/M15A_NUMERICAL_VERIFICATION.md` | Planned | Gradient/FD parity; long trajectories |
| 3 | Scaled state | `include/state_finder/scaled_initial_state.hpp` | Planned | Dimensional and coordinate-gradient tests |
| 4 | Observations/residuals | `include/state_finder/observations.hpp`, `src/cpp/objective.cpp` | Planned | Channel/units/noise/covariance tests |
| 5 | Matrix-free operators | `src/cpp/tangent_linear.cpp`, `analysis/matrix_free_normal_operator.cpp` | Planned | Matrix-free vs explicit matrix parity |
| 6 | Regularized optimization | `src/cpp/optimizer.cpp`, `include/state_finder/regularized_objective.hpp`, `docs/M27_PHASE3_CLOSURE.md` | Planned | Stationarity, realistic initialization, solver diagnostics |
| 7 | Observability SVD | `analysis/m24_observability_svd.py`, `analysis/m27_identifiability_analysis.py` | Planned | Whitened scaling, rank threshold, nullspace tests |
| 8 | Predictive uncertainty | `analysis/predictive_observability.py`, `docs/M28_OPERATIONAL_FORECAST_REPORT.md` | Planned | Uncertainty assumptions and empirical coverage |
| 9 | Modal/information decay | `analysis/information_decay.py`, `analysis/matrix_free_predictive_modes.py` | Planned optional | Mode convergence and time-scale definition |

## Numerical reference results (historical; not rerun here)

- M15A: componentwise discrete-adjoint verification and high-precision friction derivatives.
- M27: 128 fixed-budget reconstructions; none met the configured stationarity criterion.
- M28: 32 paired realizations per scenario; documented forecast inventory-error attenuation under known future boundary histories.
- M28 reports a maximum per-step mass-conservation defect of about 6.53e-10 kg in its campaign.
- These are research results, **not** claims of production performance or operational reliability.

## Source handling

Do not copy experimental algorithms wholesale. Each future transfer should record:
original path, research commit, licensing check, code changes, independent tests,
and any changed tolerances or assumptions.

The research repository is not a runtime dependency of this application.
