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

---

## Mathematical Formulation

### 1. Full-Space Constrained State Estimation

The statefinder reconstructs a transient pipeline trajectory by minimizing
discrepancies between measured and reconstructed boundary quantities,
subject to the discrete governing equations of the pipeline.

The optimization problem is formulated in full space: pipeline states
at every analysis time, boundary inputs, and reconstructed boundary
observation channels are represented explicitly.

#### 1.1 Time Discretization

Consider an analysis interval discretized as

$$
t_k = k\Delta t,
\qquad k=0,\ldots,N.
$$

Let

$$
\mathbf{x}_k \in \mathbb{R}^{2n+2}
$$

denote the discrete physical pipeline state at time $t_k$.

The state contains the pressure and mass-flow degrees of freedom
associated with the spatial discretization.

The physical boundary-input vector is

$$
\mathbf{q}_k =
\begin{bmatrix}
p_{\mathrm{in},k}\\
\dot m_{\mathrm{out},k}
\end{bmatrix}.
$$

These are the two independently prescribed physical boundary quantities
in the forward model.

The four reconstructed boundary channels are

$$
\mathbf{b}_k =
\begin{bmatrix}
p_{\mathrm{in},k}^{\mathrm{rec}}\\
\dot m_{\mathrm{in},k}^{\mathrm{rec}}\\
p_{\mathrm{out},k}^{\mathrm{rec}}\\
\dot m_{\mathrm{out},k}^{\mathrm{rec}}
\end{bmatrix}.
$$

The superscript "rec" distinguishes reconstructed boundary channels
from measured data.

The complete physical optimization trajectory is

$$
\mathbf{z} =
\left\{
\mathbf{x}_k,
\mathbf{q}_k,
\mathbf{b}_k
\right\}_{k=0}^{N}.
$$

#### 1.2 Measurement Model

Let

$$
y_{k,c}
$$

denote the observed value of boundary channel $c$ at time $t_k$.

The four channels are:

1. Inlet pressure.
2. Inlet mass flow.
3. Outlet pressure.
4. Outlet mass flow.

Let

$$
\sigma_{k,c}>0
$$

be the measurement residual scale, and define the availability indicator

$$
a_{k,c} =
\begin{cases}
1, & \text{measurement available},\\
0, & \text{measurement unavailable}.
\end{cases}
$$

For available measurements, the normalized residual is

$$
r_{k,c} =
\frac{b_{k,c}-y_{k,c}}{\sigma_{k,c}}.
$$

The measurement objective is

$$
\boxed{
J(\mathbf{z}) =
\frac12
\sum_{k=0}^{N}
\sum_{c=1}^{4}
a_{k,c}
\left(
\frac{b_{k,c}-y_{k,c}}{\sigma_{k,c}}
\right)^2
}
$$

Unavailable measurements are excluded from the objective.

The implementation permits arbitrary positive residual scales.
Interpreting these scales as measurement standard deviations requires
an appropriate statistical noise model.

#### 1.3 Discrete Pipeline Dynamics

The existing isothermal-pipe-sim forward solver defines the discrete
nonlinear governing equations.

For each time interval, write

$$
\boxed{
\mathbf{F}_k
\left(
\mathbf{x}_{k+1},
\mathbf{x}_k,
\mathbf{q}_{k+1},
\mathbf{q}_k;
\boldsymbol{\xi}
\right)
=
\mathbf{0}
}
$$

for

$$
k=0,\ldots,N-1.
$$

Here,

$$
\boldsymbol{\xi}
$$

denotes the fixed physical and numerical parameters.

The statefinder uses the actual discrete forward-model residuals.
It does not introduce a separate approximation to the pipeline physics.

This is important for discrete consistency between state
reconstruction and subsequent forward simulation.

#### 1.4 Boundary Observation Constraints

The reconstructed boundary channels are related to the physical
pipeline state and boundary inputs through an observation operator

$$
\mathbf{H}(\mathbf{x}_k,\mathbf{q}_k).
$$

The corresponding equality constraints are

$$
\boxed{
\mathbf{b}_k -
\mathbf{H}(\mathbf{x}_k,\mathbf{q}_k)
=
\mathbf{0}
}
$$

for

$$
k=0,\ldots,N.
$$

The observation constraints ensure that reconstructed measurements
correspond to physically consistent boundary quantities.

#### 1.5 Complete Optimization Problem

The full-space constrained nonlinear least-squares problem is

$$
\boxed{
\begin{aligned}
\underset{\mathbf{z}}{\operatorname{minimize}}
\quad&
\frac12
\sum_{k=0}^{N}
\sum_{c=1}^{4}
a_{k,c}
\left(
\frac{b_{k,c}-y_{k,c}}{\sigma_{k,c}}
\right)^2
\\[6pt]
\text{subject to}\quad&
\mathbf{F}_k
\left(
\mathbf{x}_{k+1},
\mathbf{x}_k,
\mathbf{q}_{k+1},
\mathbf{q}_k;
\boldsymbol{\xi}
\right)
=
\mathbf{0},
\\
&\hspace{35mm} k=0,\ldots,N-1,
\\[4pt]
&
\mathbf{b}_k-
\mathbf{H}(\mathbf{x}_k,\mathbf{q}_k)
=
\mathbf{0},
\\
&\hspace{35mm} k=0,\ldots,N,
\\[4pt]
&
z_i=z_i^{\mathrm{ref}},
\qquad i\in\mathcal{I}_{\mathrm{fixed}}.
\end{aligned}
}
$$

The final constraint represents variables deliberately held fixed
in a particular reconstruction experiment.

The set of fixed variables is configurable.

Consequently, the mathematical framework supports estimation of
initial-state variables, even though some current experiments
restrict the number of unknown initial-state coordinates.

#### 1.6 Scaled Optimization Coordinates

The implementation uses dimensionless or otherwise scaled
optimization coordinates to improve numerical conditioning.

Let

$$
\mathcal{I}_{\mathrm{free}}
$$

denote the indices of free physical variables.

For each free variable,

$$
\boxed{
z_i =
z_i^{\mathrm{ref}} + s_i u_i,
\qquad
i\in\mathcal{I}_{\mathrm{free}},
}
$$

where

- $z_i^{\mathrm{ref}}$ is the reference value used by the coordinate map;
- $s_i>0$ is the variable scale;
- $u_i$ is the corresponding optimization coordinate.

For fixed variables,

$$
z_i=z_i^{\mathrm{ref}}.
$$

The reference values used for coordinate scaling are not necessarily
the true physical state. In synthetic experiments, care is required
to avoid inadvertently supplying unknown truth through fixed variables.

The equality constraints are also multiplied by positive
constraint-specific scaling factors.

These transformations do not change the physical equality-constrained
problem, provided the scaling factors remain positive and finite.

#### 1.7 Nonlinear Optimization

The current dense reference solver uses sequential quadratic
programming (SQP).

At each iteration, the nonlinear equality constraints are linearized
and a constrained least-squares subproblem is solved.

The implementation includes:

- Exact discrete directional derivatives.
- Exact transpose derivative operations.
- Dense constrained least-squares subproblems.
- Constraint feasibility diagnostics.
- Projected-gradient stationarity diagnostics.
- Merit-function line search and convergence monitoring.

The dense solver is intended as a verified research reference.
Production-scale optimization remains an open development objective.

#### 1.8 Relationship to State Qualification

The optimization problem reconstructs a state trajectory from
measurements and the discrete governing equations.

It does not, by itself, establish whether the reconstructed state
is sufficiently accurate for subsequent engineering simulation.

That question is addressed separately through state qualification.

At each candidate analysis time, the reconstructed state is used
to initialize an independent nonlinear forward simulation.

The resulting forecast is compared against a reference trajectory
over a specified restart horizon.

Thus, the research separates two questions:

**Reconstruction:**

> Which physically consistent state trajectory best explains the
> available measurements?

**Qualification:**

> At what analysis time does the reconstructed state become
> sufficiently reliable to initialize subsequent simulations?

The second question is the principal engineering objective
of TransientPipeStateFinder.
