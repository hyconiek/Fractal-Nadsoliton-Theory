# R7N accepted results after independent repository intake

Date: 2026-09-20. Primary source: `FIN_R7N_HANDOFF_20260920`.
This register supersedes unqualified promotion language, not the immutable
source artifacts. The proofs are mathematical/interval computer-assisted
arguments, not proof-assistant formalizations or physical experiments.

## 1. Exact fixed-fixture phase census

Let the four amplitude literals be exact decimal rationals:

`r3=0.1131879146`, `r4=0.1698528641`, `r5=0.2269339093`,
`z6=-0.3380663037`.

For j=0,...,11 define

`h_j(phi)=sum_{k=3,4,5} r_k cos(2*pi*k*j/12+phi_k)/sqrt(3)
           +z6*(-1)^j/sqrt(12)`.

Write `K(phi)=log(mean_j exp(h_j(phi)))`. The uniform mean of h is zero,
and the quartic cumulant truncation is

`K4(phi)=M2/2+M3/6+(M4-3 M2^2)/24`,
where `Mn=mean_j h_j^n`.

**Accepted theorem.** Each of K4 and K has exactly 60 distinct critical points
on `(R/(2*pi Z))^3`. All are nondegenerate. Their phase-Hessian negative-index
histograms are, separately for each function,

| Negative index | 0 | 1 | 2 | 3 |
|---|---:|---:|---:|---:|
| Number of critical points | 12 | 24 | 18 | 6 |

Thus each has 12 local minima, 42 saddles and 6 local maxima as a phase
function. Do not silently identify these signs with the Hessian of a different
energy function, a constrained sphere problem, or the full Cartesian X7 dual.

### Proof chain checked in this intake

1. Recompute the local root certificates using the corrected exact-fixture
   field, true interval pi/square roots, strict Krawczyk inclusion and
   interval inertia checks.
2. At each supplied rational preconditioner verify nonsingularity and
   `sup ||I-A Hessian||_infinity<1` on the larger convex lifted phase box.
   The already-isolated small root lies inside. Uniform injectivity therefore
   proves there is no second root in that collar.
3. Verify pairwise collar separation on the torus. Quartic radii are 0.05;
   full radii range from 0.0003 to 0.0015. Recomputed largest contraction
   bounds are approximately 0.919261 and 0.997971, respectively, strictly
   below one. These decimals are explanatory, not substituted for exact bounds.
4. Recompute every gradient exclusion on the exact normalized cube `[0,1]^3`.
   For K4 there are 27,272 such cells and 640 collar cells. For K there are
   79,633 such cells and 864 collar cells.
5. Independently verify that these rectangles form complete partitions of
   the normalized torus, and that every root-cell interval is contained in
   the appropriate certified collar after a valid periodic translation.

Existence of 60 distinct local roots gives the lower count. Exhaustion of
the complement and collar uniqueness give the matching upper count.
This closes the **fixed-fixture** phase-exhaustion atom that was open in
the 2026-09-19 audit.

## 2. Validated K16/K20 surrogate route

The full-gradient exclusion does not equate a truncated polynomial to the full
log-mgf. It uses a uniform error bound for real phases.

The audit reconstructed the Fourier coefficients independently via

`E(alpha)=mean exp(alpha h)=1+sum E_n alpha^n`,
`log E(alpha)=sum L_n alpha^n`,
`L_n=E_n-(1/n)sum_{k=1}^{n-1} k L_k E_{n-k}`.

The E_n are built by exact finite-character selection on Z12, with interval
enclosures of the exact amplitude normalizations. This differs from the
source's unscaled moment/cumulant recursion. The supplied retained coefficient
intervals enclose the independently reconstructed values.

For the analytic remainder let

`H=(r3+r4+r5)/sqrt(3)+abs(z6)/sqrt(12)`.

On `|alpha|<=R`, with R>1 and RH<pi/2, the real part of E(alpha) is at least
`exp(-RH) cos(RH)>0`. Hence the phase derivative of log E is analytic there,
with magnitude bounded by

`M=R*(r5/sqrt(3))*exp(2RH)/cos(RH)`.

Cauchy's coefficient estimate gives the gradient tail after degree N bounded
by `M R^(-N-1)/(1-1/R)`. The largest amplitude is r5 in this fixed fixture.
Add the componentwise absolute Fourier sum of every omitted resonance.

Recomputed total gradient-error bounds are approximately:

- K16, 25 retained resonance pairs: `2.383381984359391e-8`;
- K20, 45 retained resonance pairs: `1.9896999978556874e-10`.

Every accepted full-gradient cell has one surrogate-gradient interval strictly
outside plus/minus the corresponding **recomputed** bound. All 54,341 K16
and 25,292 K20 exclusions passed. K20's direct and adaptive layers replace
the saved K16 residual without a gap. The proof uses derivatives with respect
to phi evaluated at phi=2*pi*z; zero exclusion is unaffected by that fixed
positive coordinate scaling.

No global alpha-homotopy theorem or absence of bifurcations along a quartic-
to-full interpolation is inferred merely from the equal endpoint counts.

## 3. Fixed-sign symmetry classification

The real-space action `j -> eps*j+a` induces
`phi_k -> eps*(phi_k+2*pi*k*a/12)` and sends z6 to `(-1)^a z6`.
Therefore the fixed negative-z6 fixture permits the 12-element subgroup
with even a and eps=+/-1, not all 24 D12 transformations inside one sign.

Numerical nearest-neighbor proposals were used only to locate a target collar.
The image of the entire small root box was then interval-certified inside
that collar. Function invariance and collar uniqueness make the resulting
permutation exact. The permutation action closes as a group.

The quartic roots form nine such orbits: eight of size six and one of size
twelve. This is a symmetry classification of the supplied fixture, not a
physical selector or a theorem for variable amplitudes.

## 4. Target P: accepted partial result only

The exact threshold is tau0=67/250. The compact hull is

`[1/900,1] x [1/128,1] x [1/9,1] x [1/1000000,1]`

in the nonnegative shared-field coordinates `(r,s,t,y)`. Accepted upstream
tails handle the complementary unbounded field regions. The current result
does not cover the entire hull.

The independent audit reconstructs the final partition from the source chain
and checks every split from exact rectangle geometry. It contains 13,231
accepted leaves and 5,432 unresolved leaves. On accepted leaves the replay
recomputes, as applicable:

- the covariance trace bound from exact outer weight boxes;
- the second elementary covariance invariant using triangle Gram areas and
  independently maximized triple-probability bounds;
- positive definiteness of a rank-three compression, using the stored rational
  basis only as a proposal and checking its rank independently;
- a centered-second-moment upper bound with a supplied rational center,
  which dominates the projected covariance for any such center.

The actual irrational powers in the weights are newly enclosed with interval
sqrt(3), rather than inferred from saved midpoint probabilities. Rational
weight and geometry bounds are rounded outward, and matrix signs are checked
with interval arithmetic. A numerical eigenvalue or saved Boolean is not the
proof decision.

The unresolved compact-volume fraction is exactly the recorded rational

`23991168798917221399653867639993952899253335329138479599153946323 /
 66072271890625000000000000000000000000000000000000000000000000000`,

approximately 0.3631049472406795. The complement therefore has relative
volume approximately 0.6368950527593205. This is coordinate volume, not
probability of truth, physical likelihood or a global error bar.

Target P remains globally open. The sharper sigma ceiling remains open.
The implication from a future global Target-P proof to a C4 index bound for
`0<g<=250/67` is valid but its global antecedent has not been established.
No partial tau0 leaf is silently promoted to sigma.

## 5. Scientific nonconclusions retained

- No amplitude-uniform phase census or identification of the exact fixture
  with an exact angular/radial coexistence state.
- No globally certified full-seven-coordinate stationary atlas, unique
  minimizing orbit, or improved global energetic bracket from this campaign.
- The existing g=5 stationary index-two counterexample is not superseded.
- No causal discord/localization map, intrinsic active gain, dimensional clock,
  selector/QW-2191 discharge, laboratory implementation, legacy completion or
  role transfer, SM/GR, L_total, or ToE closure.

## 6. Verification and provenance boundaries

All 740 original manifest entries match. The three absolute historical input
paths have matching local replacements verified by SHA-256. The supplied
portable verifier is a smoke/geometry/stored-record checker, not a full fresh
formula replay. Its fixed timeout and its cache-file manifest dependency are
portability issues, not mathematical refutations. The local audit uses new
relative-path checkers and preserves all originals.

The historical 119+7 tests are retained as baseline evidence and were not
counted again as newly executed tests in this intake. The new evidence is
the independently recomputed coefficient, root, collar, partition, full leaf,
symmetry, and partial-domain proof chain described above.
