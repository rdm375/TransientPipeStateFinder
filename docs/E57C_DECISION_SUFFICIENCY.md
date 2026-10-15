# E5.7c — Decision sufficiency without complete identifiability

## Research question

Can an engineering decision be reliable when the full
pipeline state is not locally identifiable?

This experiment tests a local mathematical proposition,
not a universal operational guarantee.

## Mathematical framework

Let A denote the constraint Jacobian and B the standardized
measurement Jacobian at a feasible nominal reconstruction.

Let Z be an orthonormal basis for ker(A).

The feasible measurement operator is

    S = B Z.

Its right singular vectors divide the feasible tangent
space into resolved and unresolved directions.

Let N span the unresolved feasible directions.

For an engineering forecast q(t), the unresolved sensitivity is

    g(t)^T N,

where g(t) is the gradient of q with respect to the
statefinder coordinates.

Incomplete identifiability means

    rank(B Z) < dim(ker(A)).

Decision-insensitive unresolved directions satisfy

    g(t)^T N = 0

at every forecast sample relevant to the decision.

Under unbounded unresolved uncertainty, nonzero sensitivity
at any sample prevents a satisfactory certificate at that
sample in the current local linearized model.

## Experimental design

### Case A — Two unknown initial pressures

Begin with the existing two-cell independent-initial-pressure
fixture.

Allow two independent initial pressure coordinates to vary
while keeping prescribed boundary controls fixed.

Generate both the truth and the reference trajectories
through independent calls to the original nonlinear
scheduled forward solver.

Verify the actual constraint nullity numerically.

Do not assume that two released coordinates automatically
produce a two-dimensional feasible tangent space.

### Case B — One measurement

Provide one standardized outlet-pressure observation.

Compute:
- Constraint rank and feasible nullity.
- Rank of BZ.
- Resolved and unresolved bases.
- Reconstructed initial-state errors.
- Forecast sensitivity along every unresolved direction.
- Decision classification with unbounded uncertainty.

Require the observed rank to be less than the feasible
dimension before labeling the case partially identifiable.

### Case C — Decision sufficiency

Evaluate a specified minimum-outlet-pressure decision.

First use the genuine nonlinear pipe forecast and its
exact tangent-linear sensitivities.

If the unresolved direction has nonzero sensitivity,
report that finding. Do not impose a numerical cutoff
to manufacture exact decision insensitivity.

Evaluate a separately labeled bounded-unresolved case,
with an explicitly justified coordinate-space radius.

A satisfactory decision with bounded unresolved uncertainty
is a different result from a satisfactory decision with
unbounded but exactly decision-insensitive directions.

### Case D — Controls

Include:
- No measurements.
- One measurement.
- Multiple independent measurements where available.
- At least one deliberately inadequate pressure limit.
- Conservation and tangent-residual checks.
- Independent nonlinear finite-difference checks.

## Acceptance criteria

1. Truth and reference are independently forward-propagated.
2. At least two feasible tangent directions are demonstrated.
3. At least one case has incomplete measurement rank.
4. Unresolved forecast sensitivities are explicitly reported.
5. Unbounded uncertainty is never silently replaced by
   a finite radius or numerical sensitivity cutoff.
6. Decision classifications are checked against their
   pressure-uncertainty envelopes.
7. All existing regression tests continue to pass.
8. Results clearly distinguish mathematical examples from
   claims about the physical pipeline model.

## Interpretation

A fully identified state is sufficient to remove local
measurement-null directions, but it may not be necessary
for a particular engineering decision.

E5.7c will determine whether this distinction can be
demonstrated with the current nonlinear pipeline fixture.

A failure to find an exactly decision-insensitive unresolved
direction is a legitimate scientific result and must be
reported without weakening the definition of sensitivity.
