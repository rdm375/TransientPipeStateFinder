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

```math
t_k = k\Delta t,
\qquad k=0,\ldots,N.
```
Let

```math
\mathbf{x}_k \in \mathbb{R}^{2n+2}
```
denote the discrete physical pipeline state at time $t_k$.

The state contains the pressure and mass-flow degrees of freedom
associated with the spatial discretization.

The physical boundary-input vector is

```math
\mathbf{q}_k =
\begin{bmatrix}
p_{\mathrm{in},k}\\
\dot m_{\mathrm{out},k}
\end{bmatrix}.
```
These are the two independently prescribed physical boundary quantities
in the forward model.

The four reconstructed boundary channels are

```math
\mathbf{b}_k =
\begin{bmatrix}
p_{\mathrm{in},k}^{\mathrm{rec}}\\
\dot m_{\mathrm{in},k}^{\mathrm{rec}}\\
p_{\mathrm{out},k}^{\mathrm{rec}}\\
\dot m_{\mathrm{out},k}^{\mathrm{rec}}
\end{bmatrix}.
```
The superscript "rec" distinguishes reconstructed boundary channels
from measured data.

The complete physical optimization trajectory is

```math
\mathbf{z} =
\left\{
\mathbf{x}_k,
\mathbf{q}_k,
\mathbf{b}_k
\right\}_{k=0}^{N}.
```
#### 1.2 Measurement Model

Let

```math
y_{k,c}
```
denote the observed value of boundary channel $c$ at time $t_k$.

The four channels are:

1. Inlet pressure.
2. Inlet mass flow.
3. Outlet pressure.
4. Outlet mass flow.

Let

```math
\sigma_{k,c}>0
```
be the measurement residual scale, and define the availability indicator

```math
a_{k,c} =
\begin{cases}
1, & \text{measurement available},\\
0, & \text{measurement unavailable}.
\end{cases}
```
For available measurements, the normalized residual is

```math
r_{k,c} =
\frac{b_{k,c}-y_{k,c}}{\sigma_{k,c}}.
```
The measurement objective is

```math
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
```
Unavailable measurements are excluded from the objective.

The implementation permits arbitrary positive residual scales.
Interpreting these scales as measurement standard deviations requires
an appropriate statistical noise model.

#### 1.3 Discrete Pipeline Dynamics

The existing isothermal-pipe-sim forward solver defines the discrete
nonlinear governing equations.

For each time interval, write

```math
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
```
for

```math
k=0,\ldots,N-1.
```
Here,

```math
\boldsymbol{\xi}
```
denotes the fixed physical and numerical parameters.

The statefinder uses the actual discrete forward-model residuals.
It does not introduce a separate approximation to the pipeline physics.

This is important for discrete consistency between state
reconstruction and subsequent forward simulation.

#### 1.4 Boundary Observation Constraints

The reconstructed boundary channels are related to the physical
pipeline state and boundary inputs through an observation operator

```math
\mathbf{H}(\mathbf{x}_k,\mathbf{q}_k).
```
The corresponding equality constraints are

```math
\boxed{
\mathbf{b}_k -
\mathbf{H}(\mathbf{x}_k,\mathbf{q}_k)
=
\mathbf{0}
}
```
for

```math
k=0,\ldots,N.
```
The observation constraints ensure that reconstructed measurements
correspond to physically consistent boundary quantities.

#### 1.5 Complete Optimization Problem

The full-space constrained nonlinear least-squares problem is

```math
\boxed{
\begin{aligned}
\min_{\mathbf{z}}
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
```
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

```math
\mathcal{I}_{\mathrm{free}}
```
denote the indices of free physical variables.

For each free variable,

```math
\boxed{
z_i =
z_i^{\mathrm{ref}} + s_i u_i,
\qquad
i\in\mathcal{I}_{\mathrm{free}},
}
```
where

- $z_i^{\mathrm{ref}}$ is the reference value used by the coordinate map;
- $s_i>0$ is the variable scale;
- $u_i$ is the corresponding optimization coordinate.

For fixed variables,

```math
z_i=z_i^{\mathrm{ref}}.
```
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


---

### 2. How Reconstruction Is Solved

#### 2.1 The Inverse Problem

The reconstruction problem is to determine a physically consistent
pipeline trajectory that best explains the available boundary
measurements.

At analysis time $t_k$, only measurements collected at or before
$t_k$ may be used.

Let the available measurement history be

```math
\mathcal{Y}_k =
\left\{
y_{j,c}:
0\le j\le k,\;
a_{j,c}=1
\right\}.
```

The corresponding reconstruction problem is

```math
\boxed{
\begin{aligned}
\min_{\mathbf{u}}\quad&
\frac12
\left\|
\mathbf{r}_k(\mathbf{u})
\right\|_2^2
\\
\text{subject to}\quad&
\mathbf{c}_k(\mathbf{u})=\mathbf{0}.
\end{aligned}
}
```

Here:

- $\mathbf{u}$ contains the scaled free optimization coordinates.
- $\mathbf{r}_k$ contains normalized residuals for available measurements.
- $\mathbf{c}_k$ contains the scaled discrete dynamics and observation constraints.

The physical trajectory is recovered through the affine coordinate
transformation defined in Section 1.6.

The reconstruction is a constrained nonlinear least-squares problem.

The optimization variables represent an entire discrete trajectory,
not merely the initial state.

#### 2.2 The Lagrangian

Introduce Lagrange multipliers $\boldsymbol{\lambda}$ associated
with the equality constraints.

The Lagrangian is

```math
\boxed{
\mathcal{L}
\left(
\mathbf{u},\boldsymbol{\lambda}
\right)
=
\frac12
\left\|
\mathbf{r}_k(\mathbf{u})
\right\|_2^2
+
\boldsymbol{\lambda}^{T}
\mathbf{c}_k(\mathbf{u}).
}
```

Define the Jacobians

```math
\mathbf{R}_k =
\frac{\partial\mathbf{r}_k}{\partial\mathbf{u}},
\qquad
\mathbf{A}_k =
\frac{\partial\mathbf{c}_k}{\partial\mathbf{u}}.
```

The first-order Karush-Kuhn-Tucker conditions are

```math
\boxed{
\begin{aligned}
\mathbf{R}_k^{T}\mathbf{r}_k
+
\mathbf{A}_k^{T}\boldsymbol{\lambda}
&=\mathbf{0},
\\
\mathbf{c}_k(\mathbf{u})
&=\mathbf{0}.
\end{aligned}
}
```

The first equation expresses stationarity of the constrained
measurement objective.

The second expresses consistency with the discrete governing equations
and the boundary observation relationships.

Satisfaction of these conditions is a necessary first-order
optimality condition under suitable constraint qualifications.
It does not by itself establish uniqueness or global optimality.

#### 2.3 Sequential Quadratic Programming

The current reference implementation uses sequential quadratic
programming (SQP).

At iteration $r$, let the current optimization coordinates be
$\mathbf{u}^{(r)}$.

The measurement residuals and equality constraints are linearized
about the current iterate.

The resulting constrained least-squares subproblem is

```math
\boxed{
\begin{aligned}
\min_{\mathbf{d}}\quad&
\frac12
\left\|
\mathbf{r}_k^{(r)}
+
\mathbf{R}_k^{(r)}\mathbf{d}
\right\|_2^2
\\
\text{subject to}\quad&
\mathbf{c}_k^{(r)}
+
\mathbf{A}_k^{(r)}\mathbf{d}
=
\mathbf{0}.
\end{aligned}
}
```

The vector $\mathbf{d}$ is the proposed optimization step.

Its associated first-order optimality system has the form

```math
\begin{bmatrix}
\mathbf{R}^{T}\mathbf{R} & \mathbf{A}^{T}\\
\mathbf{A} & \mathbf{0}
\end{bmatrix}
\begin{bmatrix}
\mathbf{d}\\
\boldsymbol{\mu}
\end{bmatrix}
=
-
\begin{bmatrix}
\mathbf{R}^{T}\mathbf{r}\\
\mathbf{c}
\end{bmatrix}.
```

Iteration superscripts have been omitted for clarity.

This is the Gauss-Newton equality-constrained subproblem.

The matrix may be singular or rank deficient when the
measurements do not sufficiently determine the free variables.

The implementation uses rank-aware dense constrained
least-squares methods rather than assuming that this
block matrix is always nonsingular.

The proposed step is incorporated through a globalization
procedure involving a merit function and line search:

```math
\mathbf{u}^{(r+1)}
=
\mathbf{u}^{(r)}
+
\alpha_r\mathbf{d}^{(r)}.
```

The step length satisfies

```math
0<\alpha_r\le 1.
```

The solver monitors:

- Measurement objective reduction.
- Equality-constraint feasibility.
- Projected-gradient stationarity.
- Linearized constraint rank.
- Measurement-information rank.
- Unresolved optimization directions.
- Step acceptance and termination status.

The current implementation is a dense reference solver.
Large-scale computational performance remains an open objective.

#### 2.4 Exact Discrete Derivatives

The reconstruction problem requires derivatives of the discrete
pipeline equations.

Rather than differentiating an alternative continuous model,
the project differentiates the actual discrete residual operators.

For a trajectory perturbation $\mathbf{v}$, the tangent-linear
operation computes

```math
\mathbf{A}_k\mathbf{v}.
```

For a constraint-space vector $\mathbf{w}$, the corresponding
transpose operation computes

```math
\mathbf{A}_k^T\mathbf{w}.
```

Analogous operations are available for the measurement residual
Jacobian.

These Jacobian-vector and transpose-Jacobian-vector products
support derivative verification, constrained optimization,
sensitivity analysis, and adjoint calculations.

They also provide a foundation for future matrix-free
large-scale optimization methods.

#### 2.5 Reconstruction Result

At the conclusion of the analysis ending at $t_k$,
the optimization produces an estimated trajectory

```math
\hat{\mathbf{z}}^{(k)}.
```

The reconstructed state at the end of the analysis interval is

```math
\boxed{
\hat{\mathbf{x}}_k
=
\mathcal{E}_k
\hat{\mathbf{z}}^{(k)}.
}
```

Here $\mathcal{E}_k$ extracts the physical state at time $t_k$.

The superscript $(k)$ is important: each analysis time may
correspond to a separate reconstruction using a different
measurement prefix.

The estimated state is the result of Reconstruction.

It is not yet a qualified state.

A small measurement residual, low constraint defect, or
successful SQP termination does not establish that the state
will produce sufficiently accurate future predictions.

---

### 3. How Qualification Is Determined

#### 3.1 Qualification Is a Predictive Question

Reconstruction determines a physically consistent state estimate.

Qualification determines whether that state is sufficiently
reliable to initialize a subsequent transient simulation.

The qualification question is therefore answered by studying
the future behavior of the reconstructed state.

For each candidate analysis time $t_k$, the statefinder
produces a reconstructed state

```math
\hat{\mathbf{x}}_k.
```

This state becomes the initial condition for a new,
independent nonlinear forward simulation.

#### 3.2 Nonlinear Forecast Continuation

Let

```math
\mathcal{S}_{k\rightarrow j}
```

denote the discrete nonlinear solution operator of the
existing forward simulator.

For $j\ge k$, the restarted trajectory is

```math
\boxed{
\hat{\mathbf{x}}_{j}^{F}
=
\mathcal{S}_{k\rightarrow j}
\left(
\hat{\mathbf{x}}_k;
\mathbf{q}_{k:j},
\boldsymbol{\xi}
\right).
}
```

The superscript $F$ identifies the independently
restarted forward trajectory.

The subsequent boundary conditions are prescribed
identically for the reference and reconstructed simulations.

The statefinder does not use future measurements when
computing $\hat{\mathbf{x}}_k$.

Future reference data are used only for retrospective
validation of the qualification result.

#### 3.3 Reference Trajectory

In synthetic validation experiments, a reference trajectory

```math
\mathbf{x}_j^*
```

is available from an independently specified forward
simulation.

The reference and restarted simulations use the same
discrete physics and future boundary conditions.

The comparison isolates the effect of errors in the
reconstructed initial state under those assumptions.

This is a controlled numerical experiment, not a test
of forward-model discrepancy.

#### 3.4 Forecast Error Measures

Let $H$ denote the restart horizon.

Define the maximum pressure error over the restarted
trajectory as

```math
\boxed{
E_p(t_k,H)
=
\max_{\substack{
j:\;t_k\le t_j\le t_k+H\\
i\in\mathcal{I}_p
}}
\left|
\hat p_{i,j}^{F}
-
p_{i,j}^{*}
\right|.
}
```

Here $\mathcal{I}_p$ denotes the pressure degrees of freedom.

Similarly, define the maximum mass-flow error:

```math
\boxed{
E_{\dot m}(t_k,H)
=
\max_{\substack{
j:\;t_k\le t_j\le t_k+H\\
i\in\mathcal{I}_{\dot m}
}}
\left|
\hat{\dot m}_{i,j}^{F}
-
\dot m_{i,j}^{*}
\right|.
}
```

The errors are measured over the complete tested
restart horizon, not merely at its endpoint.

#### 3.5 Normalized Qualification Metric

Let the engineering tolerances be

```math
\varepsilon_p>0,
\qquad
\varepsilon_{\dot m}>0.
```

Define

```math
\boxed{
Q(t_k,H)
=
\max\left(
\frac{E_p(t_k,H)}{\varepsilon_p},
\frac{E_{\dot m}(t_k,H)}{\varepsilon_{\dot m}}
\right).
}
```

The candidate state is acceptable if

```math
\boxed{
Q(t_k,H)\le 1
}
```

and both reconstruction and nonlinear restart succeed.

This is an engineering accuracy criterion, not
a statistical confidence statement.

The tolerances must be selected according to the
intended engineering application.

#### 3.6 Persistent State Qualification Time

A candidate may satisfy the qualification criterion
at one analysis time and fail at a later time.

To avoid declaring qualification based on an isolated
threshold crossing, the current implementation requires
persistent acceptance across all remaining tested candidates.

Let

```math
\mathcal{T}
=
\{t_1,t_2,\ldots,t_M\}
```

be the ordered set of candidate analysis times.

Define the acceptance indicator

```math
A_k =
\begin{cases}
1, & \text{reconstruction and restart succeed, and }
Q(t_k,H)\le 1,\\
0, & \text{otherwise}.
\end{cases}
```

The sampled persistent qualification time is

```math
\boxed{
t_Q
=
\min
\left\{
t_k\in\mathcal{T}:
A_j=1
\text{ for every tested }j\ge k
\right\}.
}
```

If no such candidate exists, the result is

**NOT QUALIFIED**.

This definition establishes persistence only across the
tested analysis times. It does not prove acceptance
at untested intermediate times or under other future
operating conditions.

#### 3.7 Interpretation of the R5.1 Experiment

The current R5.1 demonstration evaluates candidate
analysis times from 1 through 8 seconds.

The nonlinear restart horizon is 4 seconds.

The tolerances are:

```math
\varepsilon_p=10^4\;\mathrm{Pa},
\qquad
\varepsilon_{\dot m}=0.1\;\mathrm{kg/s}.
```

All tested candidates satisfy the qualification criterion.

The earliest tested qualifying time is therefore

```math
t_Q=1\;\mathrm{s}.
```

However, this does not establish a physical
statefinder spin-up time.

The experiment uses noise-free measurements,
identical forward and inverse models, and a restricted
set of unknown initial-state variables.

Its primary result is successful validation of the
causal reconstruction and nonlinear restart machinery.

#### 3.8 Qualification Without Known Truth

The current qualification metric uses a known reference
trajectory.

That is appropriate for synthetic experiments.

In an operational application, the true internal pipeline
state is generally unavailable.

Consequently, the actual future prediction errors cannot
be evaluated directly at the moment a qualification
decision must be made.

An operational qualification method will require a
defensible bound or estimate of future prediction error.

Conceptually, if an uncertainty analysis provides
conservative bounds

```math
B_p(t_k,H),
\qquad
B_{\dot m}(t_k,H),
```

then a sufficient qualification condition would be

```math
\boxed{
\max\left(
\frac{B_p(t_k,H)}{\varepsilon_p},
\frac{B_{\dot m}(t_k,H)}{\varepsilon_{\dot m}}
\right)
\le 1.
}
```

This implication is valid only if the bounds genuinely
cover the relevant prediction errors under their
stated assumptions.

Constructing and validating such bounds remains
future research.

The present implementation does not provide an
unrestricted operational qualification certificate.

#### 3.9 Relationship Between Reconstruction and Qualification

The two problems form a sequential computational procedure.

**Reconstruction**

```math
\mathcal{Y}_k
\longrightarrow
\hat{\mathbf{z}}^{(k)}
\longrightarrow
\hat{\mathbf{x}}_k.
```

**Nonlinear continuation**

```math
\hat{\mathbf{x}}_k
\longrightarrow
\left\{
\hat{\mathbf{x}}_j^F
\right\}_{j\ge k}.
```

**Qualification**

```math
\left\{
\hat{\mathbf{x}}_j^F
\right\}
\longrightarrow
Q(t_k,H)
\longrightarrow
t_Q.
```

The distinction between the two research questions is:

**Reconstruction asks:**

> Which physically consistent pipeline trajectory best explains
> the measurements available so far?

**Qualification asks:**

> When does the reconstructed state become sufficiently reliable
> to initialize subsequent engineering simulations?

The second question determines the principal engineering
result of TransientPipeStateFinder.

