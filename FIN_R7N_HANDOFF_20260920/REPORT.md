# FIN rank-seven R7N campaign report

Date: 2026-09-20

## 1. Executive scientific change

The campaign produced one primary exact fixed-fixture result and one substantial but incomplete global-curvature result.

The major new theorem is a complete phase census. For the exact declared fixed amplitudes and fixed negative alternating component, both the quartic phase function and the full log-mgf phase function have **exactly 60 critical points on the full three-torus**. In both cases the negative-index histogram is `(12,24,18,6)` for indices `0,1,2,3`. This upgrades the previous state from sixty locally certified roots plus an unresolved complement to genuine global exhaustion.

The four-amplitude practical ceiling at `tau0=67/250` remains open globally. A bounded proof campaign certified a large subset of the compact hull but retained a rigorously serialized residual occupying about `36.31049472406795%` of the compact-hull volume. No admissible counterexample was found. The sharper `sigma` target was not re-entered after this practical-target stop.

No new full-X7 stationary/global-minimizer or energetic theorem is claimed.

## 2. Audited baseline and execution state

The September 19 consolidated audit was imported read-only and its evidence hashes were rechecked. All 23 files listed in the consolidated verification retain their expected SHA-256 values.

A fresh isolated replay executed all 29 supplied baseline test files: **119/119 tests passed**. Seven integration controls also pass in the current NumPy 2.x environment after a compatibility-only normalization of `repr(np.float64(x))` before `Fraction` parsing. This normalization does not change mathematical inputs, thresholds, or proof inequalities.

The canonical safe union contains 99 direct sigma-certified domains plus 18 repair leaves covering the seven historically failed whole masks. The navigation factor 1.02 is not used as proof geometry.

## 3. Target P — rational practical ceiling

Target P is

`lambda2(M4(J)) <= 67/250`

for all nonnegative shared fields `J3,J4,J5,J6`.

The exact threshold implication `sigma < 67/250` was checked, as was the rational gain endpoint `g<=250/67`. This implication remains conditional on global Target-P closure.

The compact-hull campaign used a hierarchy of sufficient tests: imported sigma-safe mask containment, PSD trace bounds, an `e2`/Cauchy--Binet test, cell-adaptive rank-three compression, centered-moment compression, and directed refinement in the `t` coordinate. The two directed `t` passes reduced the certified residual from about 54.46% to 41.81%, then to **36.31%** of compact-hull volume.

The final residual contains 5,432 cells. At their centers, no numerical `lambda2` exceeded `0.268`; the largest observed center value was approximately `0.25392497687445353`. These center values are navigation evidence only and are not used as proof exclusions.

The second directed pass had materially smaller relative gain than the first, and subsequent pairwise compression probes had marginal hit rate. The anti-loop policy therefore stops this architecture. The next useful mathematical atom is a physical-coupled matrix enclosure retaining the common `(r,s,t,y)` dependencies, especially the `t`, `t^3`, `t^4`, and `t^(2+/-sqrt(3))` structure. A third global `t` split is not justified.

**Target P remains globally unresolved.**

## 4. Target S — sharp sigma ceiling

The sharp target remains open. Tau0-only leaves were not promoted to sigma. After the practical target retained a large residual, the campaign prioritized the independent phase lane instead of spending the bounded resource budget on a still harder sharp cover.

No new equality classification or sharp global gap is claimed.

## 5. Quartic phase theorem

The phase chart uses normalized coordinates `z=phi/(2*pi)` on `[0,1]^3` with periodic seams, while evaluations retain interval pi and the proper derivative scaling.

The old complement partition was audited: it contained 396 safe leaves, 1,272 unresolved leaves, and zero removed root-neighborhood leaves. Therefore the archived `ROOT_RADIUS=0.05` did not itself erase annuli in the stored run, but it could not be reused as a uniqueness theorem.

All sixty known quartic roots were locally recertified. Larger uniqueness collars of radius `0.05` rad were then proved by a uniform preconditioned-Jacobian contraction/injectivity argument. The old safe cells were revalidated with the corrected fixture.

A root-aware normalized-torus cover then exhausted the entire complement. The final cover has:

- 27,272 gradient-exclusion leaves,
- 640 leaves contained in certified root collars,
- zero unresolved leaves.

The fixed-sign symmetry audit permits the 12-element subgroup consisting of even translations and all corresponding reflections; odd translations flip the nonzero alternating component and are not symmetries of the same fixture. The sixty roots form nine orbits: eight of size six and one of size twelve.

Therefore the exact fixed quartic fixture has **exactly 60 critical points** with negative-index histogram `12,24,18,6`.

## 6. Full log-mgf phase theorem

All sixty previously known full-function roots were recertified locally. Larger full-function uniqueness collars were proved, with certified radii ranging from `0.0003` to `0.0015` rad. Distinct root centers have minimum torus `L_inf` separation about `0.7131676442373065` rad, and the smallest separation after subtracting both collar radii is still about `0.7123676442373066` rad. The collars are therefore pairwise disjoint by a very large margin.

Direct interval evaluation of the full gradient was too dependency-inflated near the thin annuli around roots. The successful route was a rigorously bounded cumulant/Fourier surrogate.

### 6.1 K16 layer

The order-16 surrogate retained 25 resonances. Analyticity of the log-mgf in the stated complex-alpha disk plus explicit omitted-resonance accounting gives a uniform full-gradient error bound of approximately

`2.383381984359391e-8`.

The K16 baseline cover certified 54,341 gradient-exclusion leaves and passed 5,382 residual cells to the next layer.

### 6.2 K20 layer

Increasing to order 20 and retaining 45 resonances reduced the rigorous uniform gradient error to approximately

`1.9896999978556874e-10`.

Without subdivision K20 certified 2,887 of the 5,382 K16 residual cells. An adaptive K20 cover then exhausted the remaining 2,495 parents with 22,405 further gradient-exclusion leaves and 864 root-collar leaves, leaving **zero unresolved cells**.

The full proof therefore contains 79,633 gradient-exclusion leaves across K16/K20 plus 864 root-collar leaves.

## 7. Independent replay of the full census

The complete full-function proof was checked independently of search order.

The geometry/mutation audit verifies the baseline partition, exact K16-to-K20 handoff, exact binary coverage of all 2,495 adaptive parents, every stored gradient sign relative to its declared epsilon, root-collar containment, all sixty root IDs, and the local index histogram. Deliberate mutations of a leaf, boundary, gradient sign, collar radius, source amplitude, and forbidden odd-translation symmetry are rejected.

A second layer recomputed every gradient-exclusion inequality from the surrogate formulas:

- K16: **54,341/54,341 passed**,
- K20: **25,292/25,292 passed**,
- total: **79,633/79,633 passed**, zero failures.

Thus the exact fixed full log-mgf fixture has **exactly 60 critical points** on the full phase three-torus, with the same `12,24,18,6` index histogram.

## 8. Quartic/full correspondence

The two separately exhaustive catalogs carry a one-to-one `quartic_id` label correspondence and have the same index histogram. The campaign does **not** claim a global homotopy theorem for

`K_alpha=(1-alpha)K4+alpha K_full`.

No branch tubes or complement control through every `alpha in [0,1]` were certified. The result is therefore two exact endpoint censuses plus labelled correspondence, not global topological equivalence.

## 9. Optional amplitude robustness

No amplitude-stability theorem was attempted after fixed-fixture exhaustion succeeded. Uniform margins over an amplitude box have not been paid. The exact decimal fixture is not promoted to an exact coexistence amplitude or to the different radial event.

## 10. Secondary full-X7 and energetic lane

The phase success does not provide a reduction theorem restricting all seven-coordinate stationary states or minimizers to the fixed amplitude/phase subspace. Likewise the partial four-amplitude curvature result does not transfer to X7.

For that reason no new expensive seven-dimensional stationary-complement grid was launched. The known `g=5` index-two stationary witness remains a regression guard against accidental universal index-one claims. The existing global-energy bracket is not sharpened in this campaign.

## 11. Verification summary

Final verification records:

- baseline: 119/119 tests in 29/29 files,
- integration controls: 7/7 pass,
- K16 formula replay: 54,341/54,341 pass,
- K20 formula replay: 25,292/25,292 pass,
- phase geometry/mutation audit: pass,
- cross-claim promotion/firewall audit: pass,
- consolidated intake hashes: 23/23 match.

Metadata/schema checks are not counted as mathematical proofs. The formula-level replay and coverage audits are recorded separately from provenance checks.

## 12. Final theorem status and nonconclusions

New exact fixed-fixture results:

1. Quartic phase fixture: exactly 60 critical points.
2. Full log-mgf phase fixture: exactly 60 critical points, independently replayed.

Open primary results:

1. Target P global four-amplitude ceiling at `67/250`.
2. Target S sharp global ceiling at `sigma`.

The campaign does not establish variable-amplitude phase exhaustion, a global quartic/full homotopy, full-X7 stationary exhaustion, global minimizer uniqueness, a physical gain/source law, laboratory calibration, Standard Model/gravity consequences, or Theory-of-Everything closure.

## 13. Ranked next atoms

1. **Physical-coupled Target-P enclosure:** preserve shared `(r,s,t,y)` dependencies, especially irrational `t` powers, to attack the remaining 36.31% compact residual.
2. **Optional phase amplitude stability:** only if a uniform amplitude box around the fixed fixture is scientifically useful; pay root, index, and complement margins jointly.
3. **Full-X7 exclusion primitive:** develop a genuinely seven-dimensional validated method before launching another stationary-complement campaign.

No automatic next campaign is started by this report.
