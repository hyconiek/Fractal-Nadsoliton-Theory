# Consolidated accepted rank-seven results

Repository integration audit: 2026-09-19. This document supersedes unsupported
promotion language in the imported handoffs, not their raw historical content.

## A. Conditional mathematical model

All new classical results concern the supplied strict finite operator,
`A7=X7 X7^T`, and the declared gain-dependent variational model. C4 consists
only of the cosine/alternating columns. The full probability distribution is
12-dimensional with an 11-dimensional simplex tangent. C4, X7, and a fixed
phase torus are distinct domains. The term “physical domain” in the imported
code means the specified nonnegative shared-field family, not demonstrated
laboratory physics.

The accepted upstream spectral enclosures are checked against the current
repository provider. All archives remain available with their original hashes.

## B. Baseline mathematical progress accepted in scope

1. R7P-019/020: the corrected one-dimensional resolvent has a unique minimum
   in the supplied rational r box near 0.6608386. The derivative polynomial
   has opposite strict signs on the two complementary intervals and positive
   derivative throughout the root box. Its parameter is `sech`, not `exp`.
2. R7P-026–031: local equal-energy, stability, saddle and simple-fold results
   are accepted in their declared boxes/gains. The equal-energy event is near
   3.71834489812; the fold is near 3.51564471684. These are not first-global-
   transition or connected-continuation theorems.
3. R7P-036: a genuine full-X7 stationary point at exact g=5 has two negative
   Hessian directions and no zero directions. The omitted stationary residuals
   vanish by the period-four reflection symmetry. The two 45 blocks each have
   negative determinant, while the other blocks are positive. This refutes
   the unrestricted stationary-only index-one conjecture, not a conjecture
   restricted to gains near the local equal-energy event.
4. R7P-055: the global **boundary-Ising** covariance ceiling is accepted. The
   physical closure is represented by weights `(1,2st^3,rt^4,2rst)` on the
   unit cube; the 436-leaf checked cover includes a separately paid local
   double-root neighborhood. It is not the full finite-J6 four-dimensional
   theorem. The characteristic criterion requires the certified positive
   shifted trace coefficient.
5. R7P-063: the shared-field intraparity covariance W_par obeys the ceiling.
   Both Weyl cases, the parity bound, dangerous odd-sector reduction, and
   conservative dominant-mass certificate are retained. **This is not a
   safe-region theorem for M4=W_par+bb^T by itself.**
6. R7P-068 and the later centered boxes: the local characteristic-cone
   argument is accepted where its interval signs and both Schur forms pass.
   Integral Taylor uses exact anchor identities and enclosed derivatives;
   positivity in the nonnegative cone is not claimed in all ambient directions.

The baseline tests were re-executed, not merely copied from the Python 3.13
logs. Some assertions read saved artifacts; such assertions are regression
checks only. Mathematical acceptance also uses the inspected analytic arguments
and the additional independent computations described below.

## C. Accepted FR1 projected tails and parity inequality

For nonnegative shared J3,J4,J5,J6, the following tails satisfy
`lambda2(M4)<=sigma`, with strict margin in the projected bounds:

- `exp(-J3)<=1/30`;
- `exp(-3J4/2)<=1/128`;
- `exp(-J5/2)<=1/9`.

The audit rederived the bounds by enumerating normalized C4 feature distances
after orthogonal projection, without using the imported closed-form coefficients.
For the J3 tail remove the line through F0 and F4=F8. For J4 remove the
3/5 line, bound the remaining alternating-coordinate variance by lambda6/12,
and bound the other compressed variances by anchored second moments. For J5
remove the line through F0 and F5=F7; bound the remaining label groups by
t, t², t³, t⁴. The slow t^(2-sqrt(3)) states then contribute exactly zero
to this compression. Courant--Fischer and positivity justify the resulting
second-eigenvalue bound. See `FR1_independent_geometry.json`.

The stronger parity inequality is also accepted:
`q/(1-q)>=cosh(J3)`. In the parity partition difference, the relevant even
Taylor coefficients are proportional to `4^n+2-2*3^n`, which vanish at
n=1,2 and are positive for n>=3. The odd sinh difference is nonnegative.
Nonnegative J4 and J6 can only strengthen this comparison. With
a=exp(-J3), the equivalent bound is `q(1+a)^2>=1+a^2`.

## D. Accepted FR42 large-J6 tail

The complete supplied 637-leaf boundary cover at `sigma-1/200000` was
recomputed and its binary partition, exact leaf bounds, and local quarantine
inclusions checked. Its 47 local leaves are covered by independently replayed
FR9/FR16 boxes. Their odd-mass radius exceeds `1/1000001`.

The code's assertion that C4 pair distances are cyclic was not used. All 66
feature pairs were explicitly enclosed, proving the sufficient diameter bound
D²<5. Therefore the one-sided covariance perturbation is at most 5e I.
For y=exp(-2J6)<=1/1000000, the parity bound gives

`e<=1/1000001`, and `5/1000001<1/200000`.

Thus the declared entire large-J6 tail satisfies the ceiling. This proof
does not require the unbundled FR2/FR3 proof objects, provisional FR20/FR21
reserve, or a numerical search result.

## E. Accepted FR223 union, with new subdivision certificates

The exact domain here uses `Fraction(str(value))` for the JSON-serialized
mask endpoints, without the 1.02 navigation buffer. A raw one-box replay
accepts 10 centered and 89 shifted masks. It rejects the seven whole-box
attempts named FR32, FR48, FR52, FR54, FR56, FR58 and FR60.

The audit subsequently certifies those same seven domains by subdivision:

| Original mask | Certified subrectangles |
|---|---:|
| FR32 | 2 |
| FR48 | 2 |
| FR52 | 2 |
| FR54 | 5 |
| FR56 | 2 |
| FR58 | 2 |
| FR60 | 3 |

The recorded split trees cover each original rectangle with no unresolved
leaves. This supplies a new proof object missing from the original one-box
claims. The accepted geometric union is therefore the original 106-mask union
represented by 99 direct certificates and 18 subrectangle certificates.
It is **not** its convex hull and **not** its buffered enlargement.

Each shifted leaf pays `c2>0` and `P<=0 or P1>=0`. P/P1 refer to the fixed
threshold sigma of the 3-by-3 Schur matrix. Only the threshold/inertia decision
is transferred to M4. The withdrawn `2e/25` physical-gap claim is not restored.
Outward rational rounding to 10^-60 after arithmetic operations only enlarges
enclosures and cannot create a false positive sign certificate.

## F. Corrected local phase certificates

The imported local phase programs froze some already-rounded trigonometric
constants and used a floating eigenvalue perturbation calculation. Their
reported flags alone are insufficient for the claimed exact fixture. The
replacement audit checker instead uses:

- Exact decimal amplitudes `(0.1131879146,0.1698528641,0.2269339093,-0.3380663037)`;
- interval pi and square roots in the actual twelve-label field;
- full log-mgf derivatives, and quartic derivatives directly from
  `E[h^3]/6+E[h^4]/24` (the second moment and its square are phase-independent
  by exact Fourier orthogonality);
- a rational proposed preconditioner with verified nonsingularity, strict
  Krawczyk inclusion and contraction on each rational radius-10^-7 box;
- an invertible rational change of basis followed by interval LDL inertia;
- pairwise separation of the boxes on the torus, including seams.

For each model there are at least 60 distinct certified local roots, with
negative-index counts `(12,24,18,6)` for indices `(0,1,2,3)`. The source
pairing has matching indices. No homotopy without additional roots, global
phase census, uniform amplitude robustness, or global topology equivalence
is inferred. The 1272 unresolved historical complement leaves remain
unresolved; their existence is not erased by local root recertification.

## G. Energy, angular instability, and the quantum comparison

The rational probability witness has exactly unit sum and positive entries.
Recomputation with the original rational spectral endpoints gives strictly
negative rank-seven energy at exact g=3.71835 (about -2.30816e-6). Combining
this with the **existing** ST448 lower result through A7<=A_full yields the
conditional-on-that-upstream-theorem bracket `[2.8934,3.71835]`. This audit
does not claim a new replay of the complete ST448 cover or identify the
first attaining orbit.

The pure-alternating angular instability is separately rechecked using exact
endpoint inputs on `[0.41421132290,0.41421132293]`. The defining function has
opposite signs at the endpoints and is strictly increasing for x>0; the 4/5
competitors remain negative. The 3-sine direction is below the 3-cosine
direction before its crossing. This is a constrained-sphere result, not the
radial transition.

The historical discord gap `27234855667/12500000000000` is exactly reproduced
from the existing W-spectrum provider. The apparent 5e-16 discrepancy arose
after independent rounding and subtraction of Laplacian intervals, which loses
their common row-sum dependence. No correction to the historical discord
theorem is necessary. The positive-loading interval [0.049,0.05] scaling
argument is accepted only under its existing canonical/mixed-family assumptions;
it is not a causal map from discord to localization.

## H. Cooperativity and passive results

The external ingredient was checked in the primary source: Ginibre,
*General formulation of Griffiths' inequalities*, CMP 16 (1970), Proposition 3
and Example 4. They cover real positive-definite functions on the relevant
product of finite cyclic groups with Haar measure. The four real characters
and their nonnegative supplied couplings satisfy these hypotheses. Therefore
nonnegative means/covariances and the specified isotone fixed-point iteration
are accepted in this family. This is application of an existing theorem, not
new general physics. [Primary paper](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-16/issue-4/General-formulation-of-Griffiths-inequalities/cmp/1103842172.pdf).

The already-established Hodge 11+55 identity and independent-walker generator
remain valid. The equilibrium alternating mode's excess kurtosis is -2/N,
not Gaussian at finite N. Passive Schur/Stieltjes statements retain their
symmetric passive assumptions. The exact number 28 of distinct positive
cycle eigenvalues and strict rank certification of every visible memory
residue are not newly established here from stored numerical data.

## I. Not promoted or still open

- Global positive-orthant 4D ceiling outside the accepted union and tails.
- Full-7D stationary exhaustion, gain-restricted stationary index conjectures,
  global minimizing-orbit uniqueness, and physical saddle trajectories.
- FR2/FR3, unbundled intermediate reserve/tube claims and provisional FR20/21
  as standalone replayable theorems. No imported sketch substitutes for its
  missing finite cover.
- Exact R7P-079 alignment digits solely from their stored summary; the
  referenced R7P-077 file is absent and the specialized error-bound generator
  is not supplied as a replayable certificate.
- Global phase exhaustion and amplitude-uniform phase topology.
- Active-gain provenance, QW-2191/selector discharge, dimensional units,
  laboratory apparatus/evidence, legacy completion/role transfer, SM/GR,
  L_total, or ToE closure.

The next research step must start from this corrected frontier, not from the
unqualified raw-mask list or an unsupported global summary.
