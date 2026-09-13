# FIN rank-7 — Current Research Handoff
## Post-audit continuation through certified coexistence, simple fold, stationary index-2 counterexamples, and boundary-Ising reconstruction

**Handoff date:** 2026-09-14  
**Project:** Fractal Information Theory / FIN rank-7 active-gain landscape  
**Purpose:** give a receiving AI agent a self-contained, epistemically typed continuation state after the earlier full-chat handoff and after restoration of the missing `fin_handoff_audit` package.

This file is **not** a new official FIN release and does **not** promote the conditional active-gain model to physical FIN. It is a research handoff with explicit proof status, failed conjectures, replay state, certificate paths, and recommended next work.

---

# 0. Non-negotiable guardrails

The following distinctions must remain explicit.

1. The conditional active-gain functional is
   \[
   V_g(p)=D(p\|u)-\frac g2(p-u)^T A_7(p-u)
   \]
   or equivalently its dual formulation. The source, sign and physical magnitude of `g` are **not strict-derived FIN outputs**.
2. No result below derives physical units, a physical clock, apparatus, selector, legacy role transfer, Standard Model, GR, `L_total`, or a ToE.
3. Full strict-operator results and rank-7 signed-mediator results must not be conflated.
4. Reflection-fixed / four-amplitude results do not automatically control the full seven-coordinate dual Hessian.
5. The universal full-seven-coordinate statements
   - “`index(H7)<=1 everywhere`” and
   - “`index(H7)<=1 at every stationary point for all gains`”
   are both **false**.
6. A narrower gain-restricted stationary classification near the coexistence event remains meaningful and open.
7. The correct one-dimensional resolvent coordinate is based on `sech`, not the earlier erroneous exponential substitution.

---

# 1. Starting point and replay status

The campaign started from:
- the previous full handoff `FIN_full_chat_research_handoff_pre_and_post_Discord.md`,
- the master plan `FIN_Post_Handoff_Research_Master_Plan_EN.md`,
- the restored `fin_handoff_audit/` package,
- the strict spectral input already accepted by the audit.

The restored audit files matched their historical SHA-256 records. A fresh replay of the **rank-7 intake core** was then performed:

- fresh rank-7 audit tests: **19/19 PASS**;
- regenerated `exact` layer: exactly equal to supplied `results.json`;
- regenerated manifest: exactly equal to supplied manifest.

The historical **56 inherited regressions** from other quantum/source packages were not all freshly replayed in the same local runtime because their complete local source trees were not part of the upload. They remain historical PASS, not freshly re-certified here.

The current follow-up regression suite itself now contains **38 tests**, all PASS, and the read-only verifier reports `VERIFY_PASS`.

Current verifier summary:
- tasks tracked: `128`,
- claims tracked: `31`,
- current local regression tests: `38/38 PASS`.

---

# 2. Coordinate systems frozen by the campaign

Strict positive Fourier/Laplacian values used by the rank-7 model:

\[
\lambda_3\approx1.9614068619764449,
\quad
\lambda_4\approx2.199568849333209,
\]
\[
\lambda_5\approx2.2986062720790956,
\quad
\lambda_6\approx2.3421820411463004.
\]

The four-amplitude reflection-even chart uses natural coordinates

\[
s=(s_3,s_4,s_5,s_6)
\]

with quadratic dual term

\[
\frac{\|s\|^2}{2g}.
\]

The full rank-7 dual has seven real coordinates:
- cosine/sine pair for `k=3`,
- cosine/sine pair for `k=4`,
- cosine/sine pair for `k=5`,
- one real alternating `k=6` coordinate.

The reflection-even four-amplitude chart is therefore a **proper invariant subspace** of the full seven-coordinate dual.

---

# 3. R7P-014 — exact stationary primal/dual inertia transfer

A key theorem established in the campaign:

> At a common **interior stationary point**, the full seven-coordinate dual Hessian `H7` and the 11-dimensional primal tangent Hessian have exactly the same numbers of negative and zero directions. The primal tangent Hessian has four additional positive directions.

This comes from two Schur complements of one common joint Hessian, not from sample-wise eigenvalue comparison.

Consequences:

- a certified stationary full-`H7` index-2 point is also genuinely index 2 on the primal tangent space;
- a certified full-`H7` local minimum transfers to a primal local minimum;
- this theorem applies **at stationary points**, not arbitrary fields.

Proof artifact:
`proofs/R7P-014_stationary_inertia_transfer.md`.

---

# 4. Corrected face/resolvent programme

The earlier handoff used an incorrect exponential coordinate in one place. The audit forced the corrected parameterization.

On the relevant one-dimensional face, the correct reduced coordinate is related by a `sech` law. The campaign independently reconstructed the accepted face certificates and then certified the unique one-dimensional resolvent minimum.

## R7P-019 — unique corrected 1D minimum

Certified root interval:

\[
r\in[0.6608385,\,0.6608387].
\]

Mapped natural amplitude:

\[
s_3\in[1.703818686,\,1.703819391].
\]

Conservative resolvent interval:

\[
S_{\min}\in[0.05754860,\,0.05755032].
\]

The reduced derivative numerator has degree five. Bernstein sign arguments certify one stationary point in the admissible interval, and local derivative sign certifies uniqueness.

## R7P-023 — parity first-order splitting

The two fixed-`t` first-order parity-mixing eigenvalue shifts at the dangerous double root are both interval-certified **strictly negative**.

The earlier reported reoptimized-envelope coefficient

\[
0.1312828584
\]

remains a **numerical seed**, not an upgraded theorem.

---

# 5. Certified local coexistence event — R7P-025 to R7P-029

This is one of the main new results after the restored audit.

## 5.1 R7P-026 — unique local equal-energy root

The five-equation system

\[
s/g-\mathbb E[C_4]=0,
\qquad
\Phi_g(s)=0
\]

was treated in five unknowns `(s3,s4,s5,s6,g)` with validated outward interval arithmetic and a 5D Krawczyk inclusion.

Certified center:

\[
s_3 = 1.81990358128008268\ldots
\]
\[
s_4 = 1.91398955466872407\ldots
\]
\[
s_5 = 1.91456913254684809\ldots
\]
\[
s_6 = 1.36720328019550390\ldots
\]
\[
\boxed{
g_{\rm eq}
=
3.71834489812038751367\ldots
}
\]

The certified root box uses radius approximately \(10^{-9}\) in each coordinate, and the Krawczyk image lies strictly inside that box with large inclusion margin.

Scope:

\[
\boxed{
\text{unique local reflection-even equal-energy stationary event}
}
\]

within the certified five-dimensional box.

It is **not** a theorem that this is the first global transition over the entire simplex or the only disconnected equal-energy branch.

Certificate:
`certificates/R7P-026_equal_energy_event.json`.

---

# 6. R7P-027 — full seven-coordinate stability of the localized coexistence root

The localized equal-energy root is not merely stable in the four-amplitude chart.

Exact reflection block decomposition:

\[
H_7
=
H_4
\oplus
H_{\sin},
\]

where:
- `H4` acts on `(cos3,cos4,cos5,k6)`,
- `Hsin` acts on `(sin3,sin4,sin5)`.

Interval LDL pivots are strictly positive throughout the entire R7P-026 root box for both blocks.

Hence:

\[
\boxed{
\operatorname{inertia}(H_7)=(0,0,7)
}
\]

and by R7P-014:

\[
\boxed{
\operatorname{inertia}(H_{\rm primal,tan})=(0,0,11).
}
\]

Therefore this is a genuine **local minimum** in the full rank-7 dual and in the primal tangent problem.

Certificate:
`certificates/R7P-027_localized_full_stability.json`.

---

# 7. R7P-028 — local crossing transversality

Along the stationary localized branch, the equal-energy crossing is transverse.

Certified branch-energy derivative interval:

\[
\boxed{
\frac{d\Phi_{\rm loc}}{dg}
\in
[-0.452413733,\,-0.452413730].
}
\]

Thus the local equal-energy event is not a tangent touch.

This proves a **locally unique crossing direction** in the certified neighborhood.

It does not alone give a globally connected branch all the way to the fold.

Certificate:
`certificates/R7P-028_crossing_transversality.json`.

---

# 8. R7P-029 — certified barrier saddle near coexistence

At the declared gain close to coexistence, the barrier state is interval isolated and its full seven-coordinate Hessian is certified.

Approximate saddle coordinates:

\[
s_{\rm sad}\approx
(0.94095706536,\,
1.00143940043,\,
0.96210885166,\,
0.68641498258).
\]

Full Hessian:

\[
\boxed{
\operatorname{inertia}(H_7)=(1,0,6).
}
\]

So the barrier is a genuine full-seven-coordinate index-1 saddle.

Certified dual barrier:

\[
\Phi_{\rm sad}
\in
[0.0465554318,\,
 0.0465554868].
\]

Barrier over uniform is therefore \(>0.0465554\).

A localized root at the same declared gain is also isolated; subtracting the localized objective gives a positive barrier interval around the same value.

Important scope:

- this is a local barrier saddle;
- it is **not** yet proven to be the globally lowest mountain pass;
- it is not a dynamical transition rate.

Certificate:
`certificates/R7P-029_barrier_saddle.json`.

---

# 9. R7P-030 / R7P-031 — certified simple fold

The reported numerical fold near \(g\approx3.5156447\) was rebuilt using the correct augmented 9-variable system:

\[
F(s,g)=0,
\]
\[
D_sF(s,g)v=0,
\]
\[
\|v\|^2=1.
\]

The validated 9D Krawczyk inclusion isolates one augmented root.

Certified gain interval:

\[
\boxed{
g_{\rm fold}
\in
[3.5156447068395917,\,
 3.5156447268395917].
}
\]

Approximate fold coordinates:

\[
s\approx
(1.36431142818,\,
1.43310270143,\,
1.40804570126,\,
1.00568207732).
\]

Fold nullvector approximately:

\[
v\approx
(0.50736753968,\,
0.52868565531,\,
0.55376069407,\,
0.39549810524).
\]

Certified simple-fold coefficients:

\[
\boxed{
v^T F_g
\in
[-0.2125716401,\,-0.2125716260]
<0
}
\]

and

\[
\boxed{
D^3\Phi[v,v,v]
\in
[0.1189825275,\,0.1189826671]
>0.
}
\]

The reflection-even Hessian has one null direction and three positive transverse directions.

The sine/reflection-odd block is strictly positive.

Therefore in the **full seven-coordinate model**:

\[
\boxed{
\operatorname{inertia}(H_7)=(0,1,6)
}
\]

at the fold.

This is a certified **simple stationary saddle-node**, local in the conditional rank-7 model.

It is not yet a theorem that it is the first fold globally over all stationary orbits.

Certificate:
`certificates/R7P-031_simple_fold.json`.

---

# 10. The universal stationary index-one conjecture is false

This is the central falsification result of the current campaign.

The original audit had already given a nonstationary full-seven-coordinate field with at least two negative Hessian directions. That refuted the **everywhere** index-one bound but left a stationary-only theorem logically open.

The follow-up campaign closed that loophole.

## 10.1 R7P-036 — stationary index-2 counterexample at exact g=5

In the exact invariant two-harmonic family

\[
h_j=J\cos(\pi j/2)+K(-1)^j,
\]

at exact

\[
g=5,
\]

interval-Krawczyk isolates a unique full-`X7` stationary root near

\[
J\approx1.24494819516,
\qquad
K\approx0.77955562099.
\]

On the certified box, the full Hessian decomposes into:
- a positive `k3_s` singleton,
- a positive `(k3_c,k6)` block,
- two disjoint `(k4,k5)` 2x2 blocks.

Each `(k4,k5)` determinant is strictly negative throughout the root box, below approximately

\[
-0.0202278506.
\]

Hence each contributes exactly one negative direction.

Therefore:

\[
\boxed{
\operatorname{inertia}(H_7)=(2,0,5).
}
\]

By R7P-014, the primal tangent stationary point also has Morse index 2.

At the same stationary root, the four-amplitude restricted Hessian has only index 1. The second negative direction is the sine `(k4,k5)` partner.

This proves:

\[
\boxed{
\text{“Every stationary full-H7 point has index}\le1\text{ for all gains” is false.}
\]

Certificate:
`certificates/R7P-036_stationary_index2_witness.json`.

Proof:
`proofs/R7P-036_stationary_index2_counterexample.md`.

---

# 11. R7P-038 — edge-reflection bifurcation strengthens the mechanism

A second, structurally distinct reflection-fixed family was analyzed.

A simple fold occurs near

\[
\boxed{
g\approx4.35221467865556.
}
\]

At the fold:

\[
\boxed{
\operatorname{inertia}(H_7)=(1,1,5).
}
\]

At exact gain

\[
g=4.36,
\]

two distinct stationary reflection-fixed representatives are interval isolated:
- one full-`H7` index-1 branch,
- one full-`H7` index-2 branch.

Thus the extra negative direction is not an isolated accident of the two-harmonic \(g=5\) example. A full-seven-coordinate index-2 branch is created by a certified reflection-fixed bifurcation.

Artifacts:
- `certificates/R7P-038_edge_reflection_fold.json`,
- `certificates/R7P-038_edge_reflection_branches_g4p36.json`.

---

# 12. R7P-040 — revised stationary theorem target

The campaign explicitly retires the unrestricted theorem

\[
\operatorname{index}(H_7)\le1
\quad
\text{for all stationary points/all gains}.
\]

Do not spend further proof effort trying to restore it.

The meaningful surviving question is gain/domain restricted:

> Is every relevant stationary point near the certified coexistence regime, on a declared gain interval around \(g_{\rm eq}\), of index 0 or 1?

Current evidence:
- certified localized coexistence root: index 0;
- certified coexistence barrier saddle: index 1;
- current numerical reconnaissance at gains roughly 3.72, 3.8, 4.0, 4.2, 4.3 found only indices 0 and 1;
- first certified higher-index reflection bifurcation found only around \(g\approx4.3522\).

This suggests a possible **gain-restricted topology change** rather than a universal Morse-index theorem.

Proof/status note:
`proofs/R7P-040_stationary_theorem_decision.md`.

---

# 13. Boundary-Ising reconstruction — R7P-041 to R7P-045

The compactified even-parity boundary model was rebuilt from first principles.

## 13.1 Exact four-state model and degeneracies

On even labels \(j=0,2,4,6,8,10\), define binary variables \(A,Y\in\{\pm1\}\).

The four states `(++,+-,-+,--)` correspond to label sets

\[
\{0\},\quad
\{4,8\},\quad
\{6\},\quad
\{2,10\}
\]

with multiplicities

\[
(1,2,1,2).
\]

The exact two-spin field is

\[
H_A A+H_Y Y+KAY,
\]

with

\[
H_A=J_3+J_5/4,
\]

\[
H_Y=3J_4/4-(\log2)/2,
\]

\[
K=3J_5/4.
\]

The degeneracy shift

\[
-(\log2)/2
\]

is mandatory.

Using

\[
X=e^{2J_3},\qquad
Y=e^{3J_4/2},\qquad
Z=e^{J_5/2},
\]

unnormalized probabilities are exactly

\[
\boxed{
(XYZ^4,\ 2XZ,\ Y,\ 2Z^3).
}
\]

A uniform-degeneracy two-spin model is therefore a different model.

---

# 14. Exact positive interior and corrected probability closure

For positive probabilities \(p_1,\ldots,p_4\), the physical domain \(X,Y,Z\ge1\) is exactly characterized by

\[
g_1=p_1p_4-p_2p_3\ge0,
\]

\[
g_2=4p_1p_3-p_2p_4\ge0,
\]

\[
g_3=p_1p_2^2-p_3p_4^2\ge0.
\]

Exact factorizations:

\[
g_1 S^2=2XYZ(Z^6-1),
\]

\[
g_2 S^2=4XZ^4(Y^2-1),
\]

\[
g_3 S^3=4YZ^6(X^3-1).
\]

For positive probabilities, these inequalities are also sufficient.

However:

> The weak semialgebraic inequalities are only an **outer relaxation on the zero-probability boundary**.

The true exponential-family closure has only nontrivial infinite-parameter supports:
- subsets of \(\{1,2\}\),
- subsets of \(\{1,3\}\),
- vertex \(\{1\}\),

with explicit ratio restrictions.

Therefore a future Bernstein cover over merely \(g_1,g_2,g_3\ge0\) must separately discharge spurious zero-probability boundary points.

This is a major correction to the older handoff proof architecture.

---

# 15. Boundary covariance invariants

For the exact four feature points, define

\[
P(t)=\det(tI-M)
=
t^3-e_1t^2+e_2t-e_3.
\]

The determinant is exactly

\[
\boxed{
e_3=
\frac38\lambda_3\lambda_4\lambda_5
p_1p_2p_3p_4.
}
\]

The trace coefficient is the pair-distance formula

\[
e_1=\sum_{i<j}p_ip_jd_{ij}
\]

with

\[
d_{12}=d_{34}=\frac{3(\lambda_4+\lambda_5)}8,
\]

\[
d_{13}=\frac{2(\lambda_3+\lambda_5)}3,
\]

\[
d_{14}=d_{23}
=
\frac{16\lambda_3+9\lambda_4+\lambda_5}{24},
\]

\[
d_{24}=\frac{4\lambda_3+\lambda_5}{6}.
\]

The exact formula for `e2` is stored in the proof artifact.

---

# 16. Exact eigenvalue-count criterion at sigma*

Let

\[
\sigma_*
=
\frac{
2\lambda_3(\lambda_4+\lambda_5)
-
\lambda_4\lambda_5
}{
24\lambda_3
}.
\]

Shift the characteristic polynomial:

\[
Q(z)=P(\sigma_*+z)
=
z^3+c_2z^2+c_1z+c_0.
\]

If

\[
c_2=P''(\sigma_*)/2>0,
\]

then exactly two covariance eigenvalues lie above \(\sigma_*\) iff

\[
P(\sigma_*)>0,
\qquad
P'(\sigma_*)<0.
\]

Equivalently:

\[
\boxed{
\lambda_2\le\sigma_*
\iff
P(\sigma_*)\le0
\quad\text{or}\quad
P'(\sigma_*)\ge0.
}
\]

For the strict four-state feature geometry, \(c_2>0\) globally.

This is certified by a minimum-enclosing-ball trace bound.

Correct exact radius:

\[
\boxed{
R^2
=
\frac{
16\lambda_3\lambda_4
+9\lambda_4^2
+10\lambda_4\lambda_5
+\lambda_5^2
}{
96\lambda_4
}.
}
\]

The earlier numerical value was correct; an earlier symbolic formula was not.

Strict intervals give

\[
3\sigma_*-R^2>0
\]

with margin about

\[
0.004758950459.
\]

---

# 17. Exact double-root identities — R7P-045

Define

\[
t_*^2
=
\frac{
(2\lambda_3-\lambda_4)
(2\lambda_3-\lambda_5)
}{
4\lambda_3^2
}.
\]

At

\[
p=
\left(
\frac{1+t_*}{6},
\frac{1+t_*}{3},
\frac{1-t_*}{6},
\frac{1-t_*}{3}
\right),
\]

one has exactly

\[
P(\sigma_*)=0,
\qquad
P'(\sigma_*)=0.
\]

Hence two covariance eigenvalues equal \(\sigma_*\).

The third eigenvalue is exactly

\[
\boxed{
\rho=
\frac{\lambda_4\lambda_5}{24\lambda_3}
}
\]

and

\[
\sigma_*-\rho
=
\frac{
\lambda_3\lambda_4+
\lambda_3\lambda_5-
\lambda_4\lambda_5
}{
12\lambda_3
}
>0.
\]

This point is a genuine positive-probability physical boundary point with \(Y=Z=1\) and \(X>1\), not a spurious simplex-boundary point.

Proof:
`proofs/R7P-041_045_boundary_ising.md`.

---

# 18. What is now certified versus what remains numerical

## Certified / exact / validated

- rank-7 audit intake core replay;
- coordinate/factorization consistency;
- D12 action checks;
- primal/dual stationary inertia transfer;
- corrected 1D face/resolvent minimum uniqueness;
- fixed-`t` parity first-order shifts are negative;
- local equal-energy event R7P-026;
- localized full-H7 minimum R7P-027;
- local crossing transversality R7P-028;
- barrier full-H7 index-1 saddle R7P-029;
- simple full-H7 fold R7P-031;
- stationary full-H7 index-2 counterexample at exact g=5 R7P-036;
- edge-reflection higher-index fold/branches R7P-038;
- exact boundary-Ising four-state model, positive interior domain, characteristic invariants, eigenvalue-count criterion, exact double root R7P-041..045.

## Still numerical or incomplete

- globally connected certified continuation from fold to coexistence;
- exhaustive stationary atlas on a full gain interval;
- exact gain at which full-H7 index-2 stationary branches first appear;
- global boundary-Ising theorem \(\lambda_2\le\sigma_*\) over the complete admissible boundary domain;
- intraparity `W` closure depending on that theorem;
- off-face four-amplitude curvature ceiling;
- complete 60-root phase theorem;
- global rank-7 minimizer orbit / global transition.

---

# 19. Most important failed conjectures / do-not-repeat list

Do **not** reinstate any of the following without genuinely new hypotheses.

### Refuted
1. Full seven-coordinate Hessian has index <=1 everywhere.
2. Full seven-coordinate Hessian has index <=1 at every stationary point for all gains.
3. Reflection-even H4 Morse index determines full H7 Morse index.
4. A semialgebraic closure using only \(g_1,g_2,g_3\ge0\) exactly describes the zero-probability boundary of the boundary-Ising exponential family.

### Previously refuted and still refuted
5. Symmetry alone exhausts global rank-7 minimizers.
6. One-step reflection averaging proves rank-7 globality.
7. Turning on \(s_4,s_5\) always lowers the relevant second curvature.
8. Simple global monotonicity shortcuts for the resolvent / parity mixture.
9. Passive strict memory generates active gain.

---

# 20. Updated local branch picture

The best currently justified local branch narrative is:

1. A reflection-even localized stationary pair is born in a certified simple saddle-node at
   \[
   g_{\rm fold}
   \in
   [3.51564470684,\,
    3.51564472684].
   \]
2. At the fold the full H7 has one zero and six positive directions.
3. Along the upper/localized branch there is a certified local equal-energy event at
   \[
   g_{\rm eq}
   =
   3.7183448981203875\ldots
   \]
   with a full-H7 local minimum.
4. At essentially the same gain an interval-isolated barrier saddle has full H7 index 1.
5. The crossing is transverse.
6. At larger gain, additional stationary topology appears:
   - an edge-reflection full-H7 fold near \(g\approx4.35221468\),
   - index-2 branches already certified by \(g=4.36\),
   - a stationary index-2 two-harmonic witness at exact \(g=5\).

This strongly suggests a sequence of **gain-dependent Morse topology changes**, not one universal index-one stationary theorem.

What is **not yet certified** is a continuous, gap-free branch tube joining every point from the first fold to coexistence and then onward.

---

# 21. Immediate recommended next work

The receiving agent should prioritize the following.

## Priority A — R7P-032: certified local branch diagram
Use validated continuation with overlapping tubes from the simple fold toward the equal-energy event.

Goal:
- certify connected branch segments;
- explicitly leave uncovered intervals if overlap fails;
- label fold, equal-energy event, uniform spinodal separately.

Do not draw a continuous certified hysteresis loop across missing joins.

## Priority B — refine the gain-restricted stationary theorem
Now that universal index-one is false, define a finite interval, for example below the certified edge-reflection index-2 bifurcation, and ask:

\[
\text{Are all stationary points in a declared domain/indexed atlas of index 0 or 1?}
\]

This requires bounded stationary covering, not more multistart anecdotes.

## Priority C — continue boundary-Ising R7P-046..055
The exact algebra is now solid enough to build the real proof specification.

Next:
- enumerate true zero-probability strata using the exponential-family closure;
- produce relaxed-domain negative controls;
- define exact Bernstein charts;
- cover only the physically admissible domain;
- isolate the double-root equality neighborhood separately.

## Priority D — only then propagate to intraparity W / off-face curvature
Do not use the older “boundary theorem essentially done” claim. R7P-041..045 reconstructed the algebra but did **not** yet prove the global admissible-domain inequality.

## Priority E — phase census later
The 60-root quartic/full phase correspondence remains promising but is lower priority than the now sharply defined stationary/boundary questions.

---

# 22. Key machine-readable artifacts

Current campaign root:

`/mnt/data/fin_rank7_followup/`

Core state:
- `STATE_MAP.md`
- `WORKLOG.md`
- `CLAIMS.json`
- `TASKS.json`
- `environment.json`
- `verify.py`

Key proof notes:
- `proofs/R7P-014_stationary_inertia_transfer.md`
- `proofs/R7P-023_parity_mixing_asymptotics.md`
- `proofs/R7P-033_full7_everywhere_counterexample.md`
- `proofs/R7P-034_two_harmonic_reduction.md`
- `proofs/R7P-036_stationary_index2_counterexample.md`
- `proofs/R7P-040_stationary_theorem_decision.md`
- `proofs/R7P-041_045_boundary_ising.md`

Key certificates:
- `certificates/R7P-017_021_face_certificate.json`
- `certificates/R7P-023_parity_first_order.json`
- `certificates/R7P-026_equal_energy_event.json`
- `certificates/R7P-027_localized_full_stability.json`
- `certificates/R7P-028_crossing_transversality.json`
- `certificates/R7P-029_barrier_saddle.json`
- `certificates/R7P-031_simple_fold.json`
- `certificates/R7P-036_stationary_index2_witness.json`
- `certificates/R7P-038_edge_reflection_fold.json`
- `certificates/R7P-038_edge_reflection_branches_g4p36.json`
- `certificates/R7P-041_boundary_ising.json`
- `certificates/R7P-042_boundary_ising.json`
- `certificates/R7P-043_boundary_ising.json`
- `certificates/R7P-044_boundary_ising.json`
- `certificates/R7P-045_boundary_ising.json`

Key numerical/result records:
- `results/R7P-025_candidate_records.json`
- `results/R7P-030_fold_candidate.json`
- `results/R7P-034_035_two_harmonic_stationary_atlas.json`
- `results/R7P-039_restricted_full_mismatch_table.json`
- `results/R7P-041_045_boundary_ising.json`

Current source modules:
- `src/model.py`
- `src/derivatives.py`
- `src/intervals.py`
- `src/coexistence_certificate.py`
- `src/fold_certificate.py`
- `src/two_harmonic.py`
- `src/two_harmonic_certificate.py`
- `src/edge_reflection_certificate.py`
- `src/boundary_ising.py`
- `src/face_certificate.py`
- `src/face_api.py`
- `src/parity_asymptotics.py`

Tests:
- `tests/test_B.py`
- `tests/test_C_face.py`
- `tests/test_D_coexistence.py`
- `tests/test_D_fold.py`
- `tests/test_E_edge_reflection.py`
- `tests/test_E_stationary_counterexample.py`
- `tests/test_F_boundary_ising.py`
- `tests/test_schema_and_resume.py`

Current regression status:
\[
\boxed{38/38\ \text{PASS}}
\]
with `VERIFY_PASS`.

---

# 23. Minimal receiving-agent brief

If only one paragraph is retained, retain this one:

> The restored rank-7 audit was freshly replayed at the 19-test core level and the exact layer reproduced exactly. The campaign then certified a local reflection-even equal-energy event at \(g=3.7183448981203875\ldots\), proved the localized state is a full-seven-coordinate minimum, certified an index-1 barrier saddle and transverse energy crossing, and interval-certified a simple full-H7 saddle-node at \(g\approx3.51564471684\). The major falsification is that stationary full-H7 Morse index is **not** universally <=1: an exact-gain \(g=5\) two-harmonic stationary root has certified full-H7 index 2, and a separate edge-reflection fold near \(g\approx4.35221468\) creates index-2 branches by \(g=4.36\). Therefore the next stationary theorem must be gain/domain restricted. Independently, the compactified boundary-Ising algebra was reconstructed exactly, including the mandatory degeneracy shift, exact positive-interior semialgebraic domain, the fact that its zero-probability closure is strictly smaller than the naive weak relaxation, exact covariance invariants, a sufficient eigenvalue-count criterion, and exact double-root identities at \(\sigma_*\). The next high-value work is (i) certified continuation from fold to coexistence, and (ii) a physically admissible boundary-Ising cover R7P-046..055. None of this sources active gain or closes FIN physically.

---

# 24. Restart instructions

1. Run the current regression suite before new work.
2. Run `verify.py` and require `VERIFY_PASS`.
3. Read `STATE_MAP.md`, `WORKLOG.md`, `proofs/R7P-040_stationary_theorem_decision.md`, and `proofs/R7P-041_045_boundary_ising.md`.
4. Do **not** reopen the universal stationary index-one conjecture.
5. For branch work, start at R7P-032 with the certified fold and coexistence boxes as endpoints.
6. For boundary work, start at R7P-046; do not launch a bulk Bernstein cover until the true boundary strata and physical closure are encoded.
7. Preserve strict/conditional/physical guardrails in every exported claim.
