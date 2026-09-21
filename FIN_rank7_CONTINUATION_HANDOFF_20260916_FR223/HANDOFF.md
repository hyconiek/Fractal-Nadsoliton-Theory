# FIN rank-seven continuation handoff — FR223 checkpoint — 2026-09-16

## 0. Executive state

This package continues `FIN_rank7_CONTINUATION_HANDOFF_20260915_FR17.zip` and preserves the full FR17 tree underneath the current checkpoint.

**Main theorem status:** the physical positive-orthant four-amplitude ceiling

`lambda2(M4) <= sigma_*`

is still **OPEN globally in 4D**.

No certified physical counterexample has been found. The continuation after FR17 has instead enlarged the rigorous safe union by many anisotropic/off-centre rectangles and strengthened the global large-`J6` tail.

Current numerical navigation checkpoint: **FR223**.

- `sigma_* = 0.26744324422884014`
- FR223 best saved gap: `lambda2-sigma_* = -8.537388502427579e-06`
- no positive-gap point in the four saved FR223 runs
- FR223 is **NUMERICAL_SEARCH_ONLY**, not a theorem
- current compact search domain in `(r,s,t,y)`:
  - `r in [1/900,1]`
  - `s in [1/128,1]`
  - `t in [1/9,1]`
  - `y=exp(-2J6) in [1e-6,1]`

The old FR1--FR17 handoff is preserved verbatim as `HANDOFF_PRE_FR223.md`.

## 1. Baseline and historical status

The frozen historical R7P-001--R7P-128 ledger remains unchanged. Do not renumber or reopen it. The baseline verifier in the FR17 package reported:

- `VERIFY_PASS`
- 128 tasks
- 62 claims
- one explicitly declared historical `NONREPLAYED_INPUT`

The earlier scientific conclusions remain unchanged, including the fact that no C4/four-amplitude ceiling is transferred automatically to the full seven-coordinate Hessian.

For FR1--FR17 exact details, proof notes, tests and artifact map, read `HANDOFF_PRE_FR223.md` first.

## 2. Coordinates and target

Local equality-neighborhood coordinates used throughout the continuation are

- `x = r-r_*`, with `r=exp(-2J3)`;
- `u = 1-exp(-3J4/2)`;
- `v = 1-exp(-J5/2)`;
- `e = 1-q_even`.

The residual-search variables are `(r,s,t,y)`, with

- `s=exp(-3J4/2)=1-u`;
- `t=exp(-J5/2)=1-v`;
- `y=exp(-2J6)`.

The target remains `lambda2(M4)<=sigma_*` in the physical shared-field model.

## 3. Critical methodological correction after FR17

A post-FR17 exploratory direction attempted to infer a stronger linear reserve from the 3x3 Schur-reduced matrix. The proposed global inequality

`lambda2(M4) <= sigma_* - (2/25)e`

must **NOT** be used.

Reason: the Schur matrix was constructed specifically for the threshold `sigma_*`. Its sign at that threshold is equivalent to the original spectral decision, but its numerical eigenvalue distance from `sigma_*` is not the physical quantity `lambda2(M4)-sigma_*`. Direct 4x4 checks produced physical points violating the strengthened `2/25` reserve while remaining strictly below the true ceiling.

Therefore:

- the ordinary threshold test at `sigma_*` remains valid;
- the `2/25` reserve claim is invalid and retired;
- do not use Schur-eigenvalue distance as a physical spectral gap.

A weaker extreme-face reserve with coefficient `1/30` was developed in-session, but its standalone FR20/FR21 files are not present in the final persisted overlay. Replay it before relying on it as a portable theorem.

## 4. FR42 — strengthened global large-J6 tail

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `results/FR42_global_large_J6_tail_5x.json`
- `proofs/FR42_global_large_J6_tail_5x.md`

The global large-`J6` tail was strengthened from FR5

`y=exp(-2J6) <= 1/5,000,000`

to

`y <= 1/1,000,000`.

The shifted boundary cover at reserve `delta=1/200000=5e-6` closes with:

- 637 terminal leaves
- 0 unresolved leaves
- 440 `SAFE_A_NONPOS`
- 150 `SAFE_B_NONNEG`
- 37 handed to FR9
- 10 handed to FR16

The exact covariance perturbation bound is `<5e-6` on this tail.

**Consequence:** the currently searched residual is bounded by `y>=1e-6`. This is the active global tail used by FR223.

A later attempt to push the boundary reserve to roughly `delta=1e-5` timed out in the cover. Treat that as a computational timeout, **not** a mathematical failure.

## 5. Off-centre shifted-box machinery

The key post-FR17 technical change is `src/frontier_shifted_boxes.py`.

Unlike the equality-centred local checker, this module certifies arbitrary small rectangles in `(x,u,v,e)` using second-order interval AD and midpoint Taylor enclosure of the R7P-044 characteristic/inertia disjunction.

For a rectangle it proves safety when:

- `c2.lo > 0`, and
- either `P.hi <= 0` or `P1.lo >= 0`.

This removed the repeated problem where numerical maximizers merely stuck to artificial faces of equality-centred boxes.

Important: a checker rejection is only a **proof-method negative control**, never a physical counterexample.

## 6. Current local certified/search-mask union at FR223

The FR223 artifact preserves the complete active local mask geometry used by the residual search:

- **10 centered masks**
- **96 shifted masks**
- navigation buffer: **1.02**

Files:

- `results/FR223_ACTIVE_MASK_LEDGER.json`
- `results/FR223_ACTIVE_MASK_LEDGER.csv`
- source artifact: `results/FR223_post_FR222_residual_search.json`

The 2% navigation buffer is deliberately wider than the rigorous certificates and exists only to stop the optimizer from returning a floating-point point infinitesimally outside a certified wall. It must not be confused with a proof enlargement.

### Centered masks currently encoded in FR223

FR44, FR15, FR16, FR18, FR62, FR10, FR11, FR12, FR13, FR14.

### Shifted masks currently encoded in FR223

The ledger contains 96 rectangles, including the chain FR26, FR28, FR30, FR32, ... through FR222. Use the machine-readable ledger rather than reconstructing these from chat text.

Selected high-value persisted standalone result files include:

- FR26, FR28, FR30, FR31
- FR42, FR44, FR46, FR48
- FR118, FR122, FR126, FR128, FR130, FR132, FR133
- FR178, FR182, FR184, FR191
- FR200, FR206, FR210, FR216, FR220, FR222, FR223

Not every intermediate FR result/test file survived as a separate file in the working overlay. However, the active geometry through FR222 is preserved in the FR223 mask ledger. Before promoting a missing standalone FR into a publication-grade theorem, replay that rectangle with `src/frontier_shifted_boxes.py` (or the centered checker for centered masks) and freeze a dedicated result/test artifact.

A helper script `replay_FR223_active_union.py` is included. A full 106-mask aggregate replay exceeded the interactive execution limit during handoff construction; therefore no blanket aggregate PASS is asserted here. Use sharded replay if portability-grade verification of every stored mask is required.

## 7. Major post-FR42 local milestones

The following are especially useful landmarks in the current union. Exact current bounds for all masks are in `FR223_ACTIVE_MASK_LEDGER.json`.

### Strong parity chimneys / bridges

- FR118: mid-`u,v` parity chimney, persisted standalone.
- FR122: diagonal `u-v` bridge, persisted standalone.
- FR128: high-`u`, low-`v` parity chimney, persisted standalone.
- FR130: upper-`v` parity chimney, persisted standalone.
- FR182: triple-corner bridge, persisted standalone.
- FR200: triple-corner bridge with session range `e<=1/1800`.
- FR206: diagonal FR138--FR158 bridge with session range `e<=1/1500`.
- FR220: later triple bridge, persisted standalone.
- FR222: diagonal bridge between FR188 and FR140, persisted standalone.

### Important geometric extensions

- FR126: strong right extension of the FR32 corridor.
- FR178: large right-`x` extension of thin-`u` FR56.
- FR184: left-`x` bridge around FR100/102.
- FR210: high-`u` corridor extended all the way to `v=0` on its declared box.
- FR216: upper-`u` extension of FR214.

Again, treat the JSON ledger as the canonical machine-readable current geometry.

## 8. FR223 residual search — exact saved state

Artifact: `results/FR223_post_FR222_residual_search.json`.

Status: **NUMERICAL_SEARCH_ONLY**.

Solver:

- SciPy differential evolution
- `maxiter=350`
- `popsize=18`
- `tol=1e-10`
- polish enabled
- 4 fixed seeds saved

Three runs converged to the global-tail boundary `y≈1e-6` with gap about `-9.14605e-6`.

The best saved run was seed `223034`:

```text
x = 0.00014005601545497015
u = 0.00012451172122129872
v = 0.00016451613150714195
e = 0.000010200000366878254
y = 0.000011276252044436085
J3 = 0.4554094945064724
J4 = 0.00008301298229944107
J5 = 0.00032905933154065683
J6 = 5.696405816335486
lambda2 = 0.2674347068403377
gap = -8.537388502427579e-06
```

This best run hit the iteration limit, so it is navigation only. It is still comfortably below `sigma_*` and is not a counterexample.

## 9. Exact next atom — FR224

See `NEXT_ATOM_FR224.md` for the automated wall diagnostic.

The FR223 best point is almost exactly outside the **upper `v` buffered wall shared by FR100 / FR102 / FR184**. Relative violation of that buffered `v` face is only about `1.5e-8`.

It is also extremely close to:

- the old parity wall of FR28/FR40/FR110, and
- the lower-`x` wall of FR150/FR92,

but the strongest nearby geometry is FR184, whose parity radius is much larger than the old `e<=1e-5` boxes.

**Recommended FR224:** first try a direct upward-`v` extension of FR184 while holding its current `x,u,e` bounds. If a single wide rectangle fails from dependency overestimation, use one adjacent off-centre slab crossing the FR184 `v=1/6200` wall.

Then freeze FR224 and rerun the same fixed-seed residual search with the updated mask union.

If future residual searches stop hitting local walls and repeatedly return to `y=1e-6`, switch priority back to the FR42 global-tail upgrade using the enlarged local union as the quarantine set.

## 10. Scientific interpretation at the checkpoint

The post-FR17 campaign has repeatedly shown the same pattern:

1. adversarial search returns a negative gap;
2. the point lies essentially on an artificial face of the current certified union;
3. an anisotropic/off-centre rectangle certifies that face;
4. the next search moves to another face rather than revealing a stable positive-gap interior point.

This is evidence guiding the proof search, not a proof of global validity.

At FR223 there is still **no certified physical violation** and no stable numerical positive-gap candidate.

## 11. Guardrails / do not repeat

- Do not infer the global theorem from numerical absence of a violation.
- Do not interpret checker `FAILED` as a physical counterexample.
- Do not use the invalid post-FR17 `2/25` Schur-gap claim.
- Do not interpret eigenvalue distances of the threshold-specific Schur transform as physical spectral gaps.
- Do not convexify the union of local rectangles unless separately proved.
- Do not free the physical parity variables and interpret relaxed violations as physical.
- Do not assume global monotonicity in parity weight `q`.
- Do not assume transverse monotonicity in `J4,J5` from the extreme face.
- Do not transfer this C4/four-amplitude ceiling to the full seven-coordinate Hessian.
- Keep numerical search masks and rigorous certificate bounds conceptually separate.

## 12. Replay / portability notes

Start with:

1. `HANDOFF.md` — this file.
2. `CURRENT_STATE_FR223.json` — compact machine-readable current state.
3. `NEXT_ATOM_FR224.md` — exact next-step wall diagnostic.
4. `results/FR223_ACTIVE_MASK_LEDGER.json` — complete active local geometry.
5. `results/FR223_post_FR222_residual_search.json` — latest search.
6. `src/frontier_shifted_boxes.py` — current off-centre rigorous checker.
7. `proofs/FR42_global_large_J6_tail_5x.md` + corresponding JSON — current global `J6` tail.
8. `HANDOFF_PRE_FR223.md` — full FR1--FR17 historical continuation handoff.

For baseline verification, use the original replay instructions preserved from the FR17 package.

For missing standalone post-FR17 FR files, reconstruct the exact rectangle from `FR223_ACTIVE_MASK_LEDGER.json`, replay with the generic checker, then freeze a dedicated proof/result/test before citing it independently.

## 13. Package-integrity files

`MANIFEST_FR223.sha256` is generated at packaging time and covers the files in this FR223 handoff tree (excluding the manifest itself during hash generation).
