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



---

### 4. Governing Equations and Conservation

#### 4.1 Isothermal Compressible Flow

The forward model describes transient compressible flow through a
single pipe under an isothermal approximation.

Let:

- $x$ be distance along the pipe;
- $t$ be time;
- $\rho(x,t)$ be fluid density;
- $v(x,t)$ be axial velocity;
- $p(x,t)$ be pressure;
- $\dot m(x,t)$ be mass flow rate;
- $A$ be pipe cross-sectional area;
- $D$ be pipe diameter.

For a constant-area pipe,

```math
\dot m=\rho A v.
```

Conservation of mass gives

```math
\boxed{
\frac{\partial \rho}{\partial t}
+
\frac{1}{A}
\frac{\partial \dot m}{\partial x}
=0.
}
```

A conventional one-dimensional momentum balance can be written

```math
\boxed{
\frac{\partial \dot m}{\partial t}
+
\frac{\partial}{\partial x}
\left(
\frac{\dot m^2}{\rho A}+Ap
\right)
=
-\frac{f\dot m|\dot m|}{2D\rho A}
+
\rho A g_x.
}
```

Here $f$ denotes the Darcy friction factor and $g_x$ is
the signed axial component of gravitational acceleration.

This is a continuum reference formulation. The precise
momentum terms, friction closure, and numerical approximations
used by the linked forward solver must be established by
inspection of that solver's implementation.

An isothermal equation of state closes the system:

```math
\boxed{
\rho=\rho(p,T,\boldsymbol{\xi}_{\mathrm{fluid}}),
\qquad T=\text{constant}.
}
```

For a simple constant-compressibility idealized model,

```math
\rho=\frac{p}{ZRT},
```

with consistent definitions of the specific gas constant $R$
and compressibility factor $Z$.

The repository's actual forward-model equation of state
takes precedence over this illustrative expression.

#### 4.2 Staggered Spatial Variables

The linked forward solver uses a staggered state representation.

For $n$ pipe cells, pressure is represented at $n+1$
locations and interior mass flow at $n$ locations.

The full physical state additionally contains the inlet
mass-flow degree of freedom.

Its dimension is

```math
\dim(\mathbf{x}_k)=2n+2.
```

The physical boundary inputs are inlet pressure and
outlet mass flow.

The statefinder preserves the forward solver's native
state ordering and calls its discrete residual assembly.

#### 4.3 Implicit Time Discretization

The discrete dynamics are represented abstractly by

```math
\mathbf F_k(
\mathbf x_{k+1},
\mathbf x_k,
\mathbf q_{k+1},
\mathbf q_k;
\boldsymbol{\xi}
)=\mathbf0.
```

The forward solver uses an implicit time discretization
with a time-weighting parameter $\theta$.

For a generic quantity $a$, define

```math
a^\theta
=
(1-\theta)a^k+\theta a^{k+1}.
```

The parameter $\theta$ is reserved for time integration
and is not used to denote the optimization parameter vector.

The exact discrete residual, including boundary half-cell
contributions, is defined by the external forward solver.

#### 4.4 Pipeline Inventory and Mass Conservation

Define the pipeline mass inventory

```math
\boxed{
M(t)=\int_0^L A\rho(x,t)\,dx.
}
```

Integrating continuity over the pipe gives

```math
\boxed{
\frac{dM}{dt}
=
\dot m_{\mathrm{in}}
-
\dot m_{\mathrm{out}}.
}
```

For the conservative discrete scheme, the corresponding
global mass balance is

```math
\boxed{
\frac{M^{k+1}-M^k}{\Delta t}
=
\dot m_{\mathrm{in}}^\theta
-
\dot m_{\mathrm{out}}^\theta.
}
```

This identity is expected for a converged conservative
forward step, up to numerical residual and roundoff effects.

It provides an important verification target for
reconstruction and nonlinear restart simulations.

The physical mass imbalance is

```math
\Delta\dot m
=
\dot m_{\mathrm{in}}
-
\dot m_{\mathrm{out}}.
```

The discrete conservation defect is

```math
\boxed{
\mathcal R_M^k
=
\frac{M^{k+1}-M^k}{\Delta t}
-
\left(
\dot m_{\mathrm{in}}^\theta
-
\dot m_{\mathrm{out}}^\theta
\right).
}
```

These quantities must not be confused.

A nonzero physical mass imbalance changes linepack.
It does not imply a numerical conservation error.

The conservation defect is a separate diagnostic from
measurement fitting and state qualification.

---

### 5. The Implemented Dense SQP Method

#### 5.1 Linearized Equality-Constrained Least Squares

At an SQP iterate, let

```math
\mathbf A=\frac{\partial\mathbf c}{\partial\mathbf u},
\qquad
\mathbf B=\frac{\partial\mathbf r}{\partial\mathbf u}.
```

The Gauss-Newton SQP subproblem is

```math
\boxed{
\begin{aligned}
\min_{\mathbf d}\quad&
\frac12\|\mathbf B\mathbf d+\mathbf r\|_2^2\\
\text{subject to}\quad&
\mathbf A\mathbf d=-\mathbf c.
\end{aligned}
}
```

The implemented dense linear solver handles this problem
using singular-value decompositions.

#### 5.2 Particular Feasible Step

A rank-aware particular solution is obtained from
the retained singular subspace of $\mathbf A$.

Symbolically,

```math
\boxed{
\mathbf d_p=-\mathbf A^\dagger\mathbf c,
}
```

where the dagger denotes a truncated SVD pseudoinverse.

The implementation checks whether

```math
\|\mathbf A\mathbf d_p+\mathbf c\|_2
```

satisfies its numerical compatibility tolerance.

An incompatible linearization is reported rather than
treated as an exactly feasible subproblem.

#### 5.3 Null-Space Measurement Correction

Let $\mathbf Z$ contain a basis for the numerical null
space of $\mathbf A$:

```math
\mathbf A\mathbf Z\approx\mathbf0.
```

Every compatible linearized step can be expressed as

```math
\mathbf d=\mathbf d_p+\mathbf Z\mathbf v.
```

Substituting into the measurement objective gives

```math
\boxed{
\min_{\mathbf v}
\frac12
\left\|
\mathbf B\mathbf Z\mathbf v+
\mathbf B\mathbf d_p+\mathbf r
\right\|_2^2.
}
```

The implemented solver uses another SVD to obtain
a rank-aware minimum-norm correction.

It reports:

- Numerical constraint rank.
- Reduced measurement rank.
- Number of unresolved directions.
- Retained singular-value diagnostics.
- Linearized constraint compatibility.
- Final linearized residuals.

The null-space formulation is central to understanding
measurement sufficiency.

#### 5.4 Merit Function and Globalization

The nonlinear solver uses an $L^1$ equality-constraint
penalty merit function:

```math
\boxed{
\Psi_\rho(\mathbf u)
=
\frac12\|\mathbf r(\mathbf u)\|_2^2
+
\rho\|\mathbf c(\mathbf u)\|_1.
}
```

The penalty $\rho>0$ may be increased when necessary
to obtain a sufficiently descending search direction.

For a proposed step $\mathbf d$, the line search seeks
a step length $\alpha$ satisfying an Armijo condition:

```math
\boxed{
\Psi_\rho(\mathbf u+\alpha\mathbf d)
\le
\Psi_\rho(\mathbf u)
+
\gamma\alpha
D\Psi_\rho(\mathbf u;\mathbf d),
}
```

where $0<\gamma<1$ and the directional derivative
includes the nonsmooth $L^1$ constraint contribution.

The implementation rejects unsuccessful physical trial
evaluations and backtracks.

It also detects non-descent directions and line-search
failure.

#### 5.5 Termination

The implemented convergence test combines:

1. A scaled equality-constraint feasibility tolerance.
2. A projected objective-gradient stationarity tolerance.

The projected gradient removes the component of the
objective gradient lying in the numerical range of
the transpose constraint Jacobian.

These tests establish numerical first-order consistency
to the specified tolerances.

They do not establish global optimality, uniqueness,
or qualification for future prediction.

---

### 6. Discrete Tangent-Linear and Adjoint Equations

#### 6.1 Tangent-Linear Dynamics

Consider the implicit discrete step

```math
\mathbf F_k(
\mathbf x_{k+1},\mathbf x_k,
\mathbf q_{k+1},\mathbf q_k
)=\mathbf0.
```

Define its partial derivatives:

```math
\mathbf F_{k,+}
=
\frac{\partial\mathbf F_k}{\partial\mathbf x_{k+1}},
\qquad
\mathbf F_{k,-}
=
\frac{\partial\mathbf F_k}{\partial\mathbf x_k}.
```

For perturbations in state and boundary inputs,

```math
\boxed{
\begin{aligned}
\mathbf F_{k,+}\delta\mathbf x_{k+1}
+
\mathbf F_{k,-}\delta\mathbf x_k
+
\mathbf F_{k,q+}\delta\mathbf q_{k+1}
+
\mathbf F_{k,q-}\delta\mathbf q_k
=
\mathbf0.
\end{aligned}
}
```

When $\mathbf F_{k,+}$ is nonsingular, this defines
a tangent-linear state update:

```math
\boxed{
\delta\mathbf x_{k+1}
=
-\mathbf F_{k,+}^{-1}
\left(
\mathbf F_{k,-}\delta\mathbf x_k+
\mathbf F_{k,q+}\delta\mathbf q_{k+1}+
\mathbf F_{k,q-}\delta\mathbf q_k
\right).
}
```

The inverse notation represents a linear solve.
An implementation need not form the inverse explicitly.

#### 6.2 State-Transition Sensitivity

With fixed future boundary inputs, define

```math
\mathbf T_k
=
-\mathbf F_{k,+}^{-1}\mathbf F_{k,-}.
```

The sensitivity of the future state to the initial
state is

```math
\boxed{
\mathbf\Phi_{k\rightarrow j}
=
\mathbf T_{j-1}\cdots\mathbf T_k.
}
```

Thus,

```math
\delta\mathbf x_j
\approx
\mathbf\Phi_{k\rightarrow j}
\delta\mathbf x_k.
```

This is the mathematical basis for studying
initial-state memory and perturbation attenuation.

A sensitivity attenuation threshold is not equivalent
to the nonlinear state qualification criterion.

#### 6.3 Discrete Adjoint

Suppose a scalar forecast quantity is

```math
G=G(\mathbf x_0,\ldots,\mathbf x_N).
```

For fixed boundary inputs, introduce multipliers
$\boldsymbol\lambda_k$ for the discrete dynamics:

```math
\mathcal L
=
G+
\sum_{k=0}^{N-1}
\boldsymbol\lambda_k^T\mathbf F_k.
```

The terminal adjoint equation is

```math
\boxed{
\mathbf F_{N-1,+}^T\boldsymbol\lambda_{N-1}
=
-\frac{\partial G}{\partial\mathbf x_N}.
}
```

For intermediate time levels, the backward recurrence is

```math
\boxed{
\mathbf F_{k-1,+}^T\boldsymbol\lambda_{k-1}
=
-
\frac{\partial G}{\partial\mathbf x_k}
-
\mathbf F_{k,-}^T\boldsymbol\lambda_k,
}
```

for $k=N-1,\ldots,1$.

The gradient with respect to the initial state is

```math
\boxed{
\frac{dG}{d\mathbf x_0}
=
\frac{\partial G}{\partial\mathbf x_0}
+
\mathbf F_{0,-}^T\boldsymbol\lambda_0.
}
```

This is the discrete implicit adjoint formulation
for a differentiable objective and nonsingular
linearized timestep operators.

It illustrates the backward structure of the
adjoint equations.

The repository's current full-space transpose
derivative operations provide the corresponding
building blocks. A separate production-scale
backward adjoint implementation should not be
inferred solely from these formulas.

---

### 7. Measurement Information and Identifiability

#### 7.1 Physically Feasible Perturbations

At a feasible reconstructed trajectory, let

```math
\mathbf A
=
\frac{\partial\mathbf c}{\partial\mathbf u}
```

be the scaled equality-constraint Jacobian.

A first-order physically feasible perturbation
must satisfy

```math
\mathbf A\delta\mathbf u=\mathbf0.
```

Let $\mathbf Z$ span its null space:

```math
\delta\mathbf u=\mathbf Z\mathbf v.
```

#### 7.2 Reduced Measurement Jacobian

Let

```math
\mathbf B
=
\frac{\partial\mathbf r}{\partial\mathbf u}
```

be the normalized measurement-residual Jacobian.

The measurement sensitivity restricted to
feasible perturbations is

```math
\boxed{
\mathbf B_{\mathrm{red}}=\mathbf B\mathbf Z.
}
```

Its singular values describe the local
distinguishability of feasible perturbations.

A zero singular value corresponds to a locally
unresolved feasible direction.

A small singular value corresponds to a
poorly constrained direction relative to the
chosen variable and measurement scaling.

#### 7.3 Local Information Matrix

For independent Gaussian measurement errors with
correctly specified residual scales, the
Gauss-Newton local information matrix is

```math
\boxed{
\mathcal I
=
\mathbf Z^T\mathbf B^T\mathbf B\mathbf Z.
}
```

The information matrix is positive semidefinite.

It is positive definite only when the reduced
measurement Jacobian has full column rank.

Without that rank condition, the data do not
locally determine every feasible perturbation.

These statements concern local linearized
information and depend on the measurement model.

They do not establish global identifiability.

#### 7.4 Reconstruction Uncertainty

If the reduced information matrix is positive
definite and the assumed statistical model is
appropriate, a local covariance approximation is

```math
\boxed{
\mathbf C_v\approx\mathcal I^{-1}.
}
```

The corresponding covariance in physical
optimization coordinates is

```math
\mathbf C_u
\approx
\mathbf Z\mathbf C_v\mathbf Z^T.
```

This approximation is not valid as an unrestricted
finite covariance when the information matrix
is singular.

A pseudoinverse describes uncertainty only on
the retained identifiable subspace and does not
resolve uncertainty in the null directions.

Prior information, physical bounds, or additional
measurements would be required to constrain
those directions.

---

### 8. From Reconstruction Uncertainty to Prediction Uncertainty

#### 8.1 Linearized Forecast Propagation

Let $\mathbf C_k$ represent a defensible covariance
estimate for the reconstructed physical state.

For fixed future boundary inputs, the
local-linearized forecast covariance is

```math
\boxed{
\mathbf C_j
\approx
\mathbf\Phi_{k\rightarrow j}
\mathbf C_k
\mathbf\Phi_{k\rightarrow j}^T.
}
```

If future boundary conditions are uncertain,
their contribution must also be included.

This equation does not account automatically
for model discrepancy, non-Gaussian uncertainty,
or nonlinear effects beyond the local approximation.

#### 8.2 Prediction-Weighted Observability

Not every poorly observed state direction is
important to an engineering prediction.

Let $\mathbf G$ be the Jacobian of a chosen
future prediction with respect to the
full-space optimization coordinates.

For a feasible perturbation,

```math
\delta\mathbf g
\approx
\mathbf G\mathbf Z\mathbf v.
```

The measurements respond according to

```math
\delta\mathbf r
\approx
\mathbf B\mathbf Z\mathbf v.
```

A prediction-weighted sensitivity ratio is

```math
\boxed{
\eta(\mathbf v)
=
\frac{
\|\mathbf G\mathbf Z\mathbf v\|_2
}{
\|\mathbf B\mathbf Z\mathbf v\|_2
}.
}
```

This ratio is defined only when its denominator
is nonzero.

Directions invisible to the measurements require
separate treatment.

In particular, if

```math
\mathbf B\mathbf Z\mathbf v=\mathbf0
```

but

```math
\mathbf G\mathbf Z\mathbf v\ne\mathbf0,
```

then the available measurements cannot distinguish
a perturbation that affects the future prediction.

Such a direction is especially important for
state qualification.

Conversely, a poorly observed direction that
has negligible influence on relevant future
predictions may be tolerable.

This distinction motivates prediction-oriented
rather than exclusively state-oriented
measurement sufficiency.

#### 8.3 A Local Prediction-Error Bound

Suppose the feasible perturbation coordinates
satisfy an ellipsoidal uncertainty constraint

```math
\mathbf v^T\mathbf W\mathbf v\le1,
```

where $\mathbf W$ is positive definite.

For a scalar linearized prediction with
sensitivity row $\mathbf g^T$, the worst-case
absolute prediction error over this ellipsoid is

```math
\boxed{
B_g
=
\sqrt{
\mathbf g^T
\mathbf W^{-1}
\mathbf g
}.
}
```

This is exact for the stated linearized prediction
and ellipsoidal uncertainty set.

It is not automatically a rigorous bound for
the nonlinear forward model.

Nonlinear remainder estimates or direct nonlinear
validation are needed before using it as a
guaranteed engineering certificate.

---

### 9. Robust State Qualification

#### 9.1 Dependence on Future Controls

The current synthetic qualification metric
depends on the specified future boundary conditions.

More explicitly, write

```math
Q(t_k,H;\mathbf q^+),
```

where $\mathbf q^+$ is the future
boundary-control history.

A state may qualify under one operating
scenario but fail under another.

#### 9.2 Scenario-Robust Qualification

Let $\mathcal U$ denote an admissible
set of future boundary-control scenarios.

A proposed robust qualification metric is

```math
\boxed{
Q_{\mathrm{robust}}(t_k,H)
=
\sup_{\mathbf q^+\in\mathcal U}
Q(t_k,H;\mathbf q^+).
}
```

This extends the current single-scenario
synthetic experiment.

The supremum must not be described as evaluated
unless the scenario set has actually been
searched or rigorously bounded.

#### 9.3 Qualification Without Known Truth

Operational qualification generally cannot
compare the reconstructed state with an
unknown true state.

Instead, suppose a defensible uncertainty
set $\mathcal X_k$ contains all initial
states considered possible given the
measurements and modeling assumptions.

Let

```math
\hat{\mathbf x}_k
```

be the reconstructed state.

For each candidate state
$\mathbf x\in\mathcal X_k$ and future
control history $\mathbf q^+\in\mathcal U$,
compare the nonlinear forecasts initialized
from $\hat{\mathbf x}_k$ and $\mathbf x$.

A prospective robust criterion is

```math
\boxed{
\sup_{\mathbf x\in\mathcal X_k}
\sup_{\mathbf q^+\in\mathcal U}
\max\left(
\frac{E_p(\hat{\mathbf x}_k,\mathbf x;\mathbf q^+,H)}
{\varepsilon_p},
\frac{E_{\dot m}(\hat{\mathbf x}_k,\mathbf x;\mathbf q^+,H)}
{\varepsilon_{\dot m}}
\right)
\le1.
}
```

Here the errors compare the two nonlinear
forecast trajectories over the complete horizon.

This criterion would provide a deterministic
conditional guarantee only if the true state
belongs to $\mathcal X_k$ and the forward-model
assumptions hold.

Constructing credible uncertainty sets and
evaluating or bounding this supremum are
future research tasks.

The current R5.1 implementation does not
provide this operational certificate.

#### 9.4 Statistical Qualification

A different research direction is to define
a probability-based criterion, such as

```math
\Pr\left(
Q(t_k,H)\le1
\mid\mathcal Y_k
\right)
\ge1-\alpha.
```

Such a statement requires a specified
probabilistic model for the unknown state,
measurement errors, and any relevant
model or control uncertainty.

It is not interchangeable with a
deterministic worst-case guarantee.

No calibrated statistical qualification
claim is established by the current
synthetic experiments.

---

### 10. Verification and Engineering Acceptance

#### 10.1 Separate Numerical and Engineering Questions

The implementation should distinguish:

**Numerical correctness**

- Discrete residual consistency.
- Tangent-linear derivative verification.
- Adjoint transpose consistency.
- Constraint feasibility.
- SQP stationarity.
- Global mass conservation.
- Nonlinear restart reproducibility.

**Engineering reliability**

- Forecast pressure error.
- Forecast mass-flow error.
- Sensitivity to initial-state deficiencies.
- Sensitivity to measurement noise.
- Sensitivity to future operating scenarios.
- Qualification time and persistence.
- Uncertainty coverage when truth is unavailable.

Passing numerical tests is necessary for
trustworthy engineering results, but does
not by itself establish engineering qualification.

#### 10.2 State Qualification Time

For a fixed test configuration, the sampled
qualification time remains

```math
\boxed{
t_Q
=
\min\left\{
t_k:
Q(t_j,H)\le1
\text{ for every subsequently tested }t_j,
\text{ with successful reconstruction and restart}
\right\}.
}
```

If no such candidate exists, the
result is NOT QUALIFIED.

The result must be reported together with:

- The analysis-time sampling grid.
- The restart horizon.
- Pressure and mass-flow tolerances.
- Measurement availability and noise assumptions.
- Initial-state deficiency.
- Future boundary-control scenario.
- Forward-model and discretization provenance.
- Numerical convergence and conservation diagnostics.

#### 10.3 Central Research Hypothesis

The central research hypothesis is:

> A pipeline state need not be uniquely determined
> in every coordinate to support reliable engineering
> predictions. Qualification depends on whether the
> uncertainty remaining after reconstruction can
> materially affect the subsequent simulated behavior.

The research program tests this hypothesis by
combining measurement-information analysis,
discrete sensitivity propagation, and independent
nonlinear forecast restarts.

The current R5.1 experiment validates the
qualification workflow but does not yet establish
a nontrivial physical qualification time.

R5.2 will investigate initially deficient states,
expanded initial-state unknowns, controlled
measurement noise, and persistent qualification
under nonlinear forecast validation.

