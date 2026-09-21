# R7O3 mathematical and implementation audit

Date: 2026-09-21. Source package: `FIN_R7O3_TARGETP_HANDOFF_20260920`.
Global acceptance requires the completed gates in `verification.json`; an
unfinished replay or this analytic note alone does not establish coverage.

## 1. Exact model and eigenvalue convention

The supplied model has twelve labels j=0,...,11 and four normalized features

```text
F_j = (sqrt(lambda3/6) cos(pi j/2),
       sqrt(lambda4/6) cos(2 pi j/3),
       sqrt(lambda5/6) cos(5 pi j/6),
       sqrt(lambda6/12) (-1)^j).
p_j = exp(J3 c3_j+J4 c4_j+J5 c5_j+J6 (-1)^j) / normalization,
J3,J4,J5,J6 >= 0,
M4 = Cov_p(F),                  tau = 67/250.
```

The target means **at most one eigenvalue of M4 is strictly larger than tau**.
Equivalently, `lambda2(M4)<=tau` when eigenvalues are ordered decreasingly.
The source theorem's phrase “ordered increasingly” is inconsistent with its
order-free statement and proof. The accepted convention corrects this wording;
the original archived file is preserved.

Spectral intervals are those of the accepted finite operator, not arbitrary
free spectral parameters. The source's rounded feature intervals are checked
to contain a reconstruction from the accepted spectral inputs and exact
trigonometric values. No kernel/gain provenance is supplied by this audit.
The current strict-spectrum provider is also rerun and its exact Laplacian
intervals are compared with those frozen inputs.

## 2. Centered moment and min–max argument

For a fixed rational rank-three matrix B of shape 4 by 3, and a fixed c in Q^3,
set Z=B^T F. The identity

`E[(Z-c)(Z-c)^T] - Cov(Z) = (E[Z]-c)(E[Z]-c)^T >= 0`

holds with no optimality requirement on c. If throughout a parameter cell

`K = tau B^T B - E[(Z-c)(Z-c)^T] > 0`,

then `v^T(tau I4-M4)v>0` on the three-dimensional subspace range(B).
Courant–Fischer implies that the second largest eigenvalue is strictly below
tau on that cell. The accepted global statement is the conservative non-strict
ceiling after joining all inherited proof types and tails.

The basis need not be orthonormal. The Gram matrix must remain B^T B, not I3.
Exact minors prove rank, and the source center is retained unchanged. Floating
eigenvectors or optimizers are not called during this intake's formula replay.

## 3. Shared fields and exact aggregation

Use `r=exp(-2J3)`, `s=exp(-3J4/2)`, `t=exp(-J5/2)`, `y=exp(-2J6)` and A=sqrt(r).
Relative to the anchor label j=0, exact aggregation in label order
`[0,4,6,2,3,5,1]` gives multiplicities `[1,2,1,2,2,2,2]` and weights

```text
1, 2 s t^3, A^2 t^4, 2 A^2 s t,
2 A t^2 y, 2 A s t^(2-sqrt(3)) y, 2 A s t^(2+sqrt(3)) y.
```

`analytic_checks.py` derives these powers in Q(sqrt(3)) from the twelve label
cosines. They share the same four variables and one normalization denominator.
The denominator is at least one because of the anchor weight. Treating these
seven weights as unrelated free variables would be a different relaxation.

## 4. Validated interval jets and Taylor remainder

The implementation differentiates in `(A,u,v,y)=(sqrt(r),1-s,1-t,y)`.
The transformed interval contains the image of the original cell; its exact
rational midpoint and radii are used. Feature parameters, B and c are constant
with respect to these four derivatives.

The implemented product and inverse rules are the usual second derivative
identities, including both mixed gradient terms. The only noninteger power
is `z=(1-v)^sqrt(3)`; the other weights use z or 1/z. Its value is enclosed
using `19/11 < sqrt(3) < 26/15`, checked by exact squared rational comparisons.
Rational-power endpoints are verified by integer powers: a floating initial
guess is only a locator and cannot certify an endpoint.

For `0<t<=1`, the exact derivatives satisfy

```text
dz/dv = -sqrt(3) t^(sqrt(3)-1),
d²z/dv² = sqrt(3)(sqrt(3)-1) t^(sqrt(3)-2).
```

Since `0<sqrt(3)-1<1` and `-1<sqrt(3)-2<0`, valid conservative enclosures are
`sqrt(3)*t <= |dz/dv| <= sqrt(3)` and
`1 <= t^(sqrt(3)-2) <= 1/t_lo`. These bounds are for derivatives of the same
true irrational-power function; they are not derivatives of an approximating
rational exponent. The sign of the first derivative is retained.

Every centered-moment entry f has the enclosure

```text
f(m) + grad(f)(m) dot (X-m)
  + [-R,R],
R = (1/2) sum_ij sup_X |partial_ij f| radius_i radius_j.
```

The rectangle is convex, so the integral/mean-value Taylor remainder is valid.
Outward enclosures of midpoint values and gradients only enlarge this bound.
Positive definiteness is established by strict Sylvester minors or Gershgorin
lower bounds. The finalizer also reconstructs K from the newly saved moments
and checks positivity using exact rational interval determinant arithmetic;
it does not accept an `ok=true` flag or saved PD diagnostics alone.

## 5. Arithmetic and independence of verification

The full fresh replay uses the inspected R7O3 jet formulas with the previously
audited `fin_r7o2_review/intervals_fast.py` backend. Each elementary binary64
operation is enlarged with nextafter in the appropriate directions. Public
endpoints are exact Fractions, avoiding loss of direction when computing
midpoints and Taylor radii. Nonfinite values and zero-containing denominators
are rejected. The original rational-grid rounding initializer is disabled
only in the in-memory fast adapter, not edited in the archived source.

This is an independent arithmetic implementation, not a second independently
written implementation of every jet formula and not proof-assistant
formalization. The analytic derivative rules are inspected separately.

A separate process replays 21 selected source witnesses on the original
outward rational 10^-9 backend. The sample includes every one of the 14
source Gershgorin certificates and the reported weakest-margin leaf 10878.
Both the inequalities and exact stored rational bounds match in this sample.
Do not describe this as a rational replay of all 12,425 leaves.

An independent direct twelve-label exponential calculation also checks all
16 corners and the midpoint of each sampled cell: 357 points and 2,142
moment entries. These are cross-formulation controls, not domain coverage.

## 6. Coverage, inherited proofs and whole-domain conclusion

Each of the 5,432 original residual parents is matched by exact cell and path
to the accepted R7N input. An independent recursive traversal verifies each
split axis, split location, both children and every terminal. Duplicate IDs,
missing leaves, residual nodes and positive-volume overlap/gap are rejected.
Every active leaf is used exactly once, and there are 12,425 such leaves.

The imported compact partition is compared cell-for-cell with a reconstruction
from the original accepted R7N checkpoints. It has 18,663 original cells:
13,231 previously accepted SAFE cells and 5,432 former residual parents.
Expanding the latter yields 25,656 terminal cells. The exact tree reconstructs

`[1/900,1] x [1/128,1] x [1/9,1] x [1/1000000,1]`.

Equal volume alone is not used as the partition proof. Adjacent cells may
share faces; they do not overlap in their interiors.

Outside this hull at least one accepted tail applies: FR1 in r, s or t,
or FR42 in y. Exact lower faces remain in the compact hull, so no interface
is dropped. The tail premise hashes and accepted independent audit artifacts
are checked. Their sharper threshold

`sigma=(2 lambda3 (lambda4+lambda5)-lambda4 lambda5)/(24 lambda3)`

is independently enclosed below 67/250. The prior SAFE and tail results are
reused as accepted premises, not advertised as newly re-proved here.

## 7. Consequence and explicit exclusions

For the Cartesian C4 Hessian `H4=I4/g-M4`, the global ceiling implies
`index_negative(H4)<=1` for supplied `0<g<=250/67`. At the endpoint the
non-strict bound does not exclude zero modes. Below the endpoint it supplies
three positive transverse Hessian directions, not positivity of all four.

The source's smallest PD-test margin is not a normalized spectral gap; no
improved global constant is extracted from it. Target S is not established.
Neither an everywhere nor a stationary-only full-X7 theorem follows, and
existing full-X7 counterexamples remain valid. Global energy ordering,
minimizing-orbit uniqueness, gain/clock origin, selector/QW-2191, laboratory
evidence, legacy completion/role transfer, SM/GR, L_total and ToE closure are
not consequences of this intake.

## 8. Reproducibility corrections

The source `verify.py` checks stored formula-completion and mutation ledgers;
it does not recompute all moment formulas. Its success alone is insufficient.
The source replay README uses named command-line options, but the actual
`replay_clean_math_multi.py` reads positional arguments. It also recommends
removing a checkpoint directory. The intake instead exports working commands
that write only into `fin_r7o3_review`, leaving the source and its checkpoints
untouched. These are reproducibility defects, not mathematical refutations.

The 12 producer mutation results are retained as historical evidence. The
intake runs 15 new controls plus the six inherited arithmetic/geometry tests;
it does not relabel all producer tests as independently repeated end-to-end
tests. In particular, fabricated global completion cannot bypass the required
full leaf-count, identity and inequality gates.

All 12,425 certificates declare historical checker hash
`1f13e7d768a7e9a4f3888d19749118096f48b4413100ffd322b9db95fea92da0`,
which is not the shipped checker's hash
`88bbd4ea7a160ba3e9986ee0d5577c61e5f174be8032993a698f5bac33f66728`.
The source does not establish that these are byte-identical implementations.
That metadata link is not accepted as current-code provenance. Instead the
intake replays the inspected shipped checker, records its actual hash, and
exports new bound records tied to the unchanged source certificates. The
fresh proof does not depend on reconstructing the historical checker.
