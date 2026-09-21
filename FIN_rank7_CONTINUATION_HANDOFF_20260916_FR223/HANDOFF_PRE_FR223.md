# FIN rank-seven continuation handoff — 2026-09-15

## 0. Identity and stopping state

This handoff continues the completed 128-task campaign from `HANDOFF_BASELINE_20260914.md`.
The historical R7P ledger is **frozen**: R7P-001 through R7P-128 remain terminal and are not reopened or renumbered here.

The continuation attacks the first ranked frontier from the baseline handoff:

> close or refute the residual compact core of the physical positive-orthant four-amplitude ceiling `lambda2(M4)<=sigma_*`.

Current stopping state: **the global 4D ceiling is still open**. No certified admissible counterexample has been produced. The continuation has, however, sharply reduced the residual domain and enlarged the rigorously certified equality-neighborhood geometry.

Baseline verifier replayed in this continuation tree:

- `VERIFY_PASS`
- tasks: 128
- claims: 62
- one historical source ZIP remains explicitly `NONREPLAYED_INPUT`, exactly as in the baseline handoff.

Current durable continuation test layer:

- **21 passed, 0 failed** across FR1 and FR8--FR14 sharded tests.
- See `results/FRONTIER_REPLAY_20260915.json` and `logs/frontier_replay_20260915.log`.

## 1. Baseline results that remain unchanged

The following conclusions from the 2026-09-14 handoff remain in force and are not upgraded by this continuation:

1. The exact boundary-Ising ceiling is globally certified on its declared boundary closure.
2. The shared-field intraparity ceiling is certified on its declared domain.
3. A stationary full-seven-coordinate index-at-least-two witness at `g=5` refutes the stationary-only universal index-one conjecture.
4. The full-seven-coordinate stationary atlas is discovery/saturation only, not exhaustion.
5. The rigorous global energetic bracket remains `g_global in [2.8934,3.71835]`; the attaining orbit is not identified.
6. The 60 quartic and 60 nearby full phase roots are locally isolated, but global phase exhaustion remains open.
7. No gain source, clock, selector, apparatus, or causal quantum-to-localization bridge has been derived.
8. No result in this continuation transfers the C4/four-amplitude ceiling to the full seven-coordinate Hessian.

## 2. Coordinates used by the continuation

The local equality-neighborhood coordinates are

- `x = r-r_*`, with `r=exp(-2 J3)` and signed `x`;
- `u = 1-exp(-3 J4/2)`;
- `v = 1-exp(-J5/2)`;
- `e = 1-q_even`.

The target remains the physical shared-field four-amplitude covariance ceiling

`lambda2(M4) <= sigma_*`.

The core local checker is the same second-order rational interval-AD plus R7P-044 characteristic/inertia Schur-cone criterion used by R7P-068. A checker rejection is retained as a **method negative control**, never as a physical counterexample unless an independent admissible eigenvalue violation is actually certified.

## 3. New durable results

### FR1 — projected global tails and stronger parity constraint

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR1_residual_tail_upgrade.md`
- `results/FR1_residual_tail_upgrade.json`
- `src/frontier_residual_tails.py`
- `tests/test_FR1_residual_tails.py`

New globally safe tails:

- `exp(-J3) <= 1/30`;
- `exp(-3 J4/2) <= 1/128`;
- `exp(-J5/2) <= 1/9`.

The `J5` result is a major strengthening of the old `2^-11` tail. It uses Courant--Fischer projection to remove the slow `j=5,7` states because they have the same C4 feature vector.

FR1 also proves the stronger exact physical parity constraint

`q/(1-q) >= cosh(J3)`,

or with `a=exp(-J3)`,

`q(1+a)^2 >= 1+a^2`.

Therefore every still-unresolved point must in particular satisfy

- `exp(-J3) > 1/30`,
- `exp(-3 J4/2) > 1/128`,
- `exp(-J5/2) > 1/9`,

in addition to lying outside all earlier certified boundary/face/local regions.

### FR2 — strengthened determinant reserve on the full extreme face

Status: **PROOF NOTE PRESERVED; standalone replay script missing in this continuation tree**.

Artifact: `proofs/FR2_extreme_face_reserve.md`.

On `J4=J5=0`, in the dangerous branch,

`R = det(sigma I-Mtilde)/det(sigma I-Wpar)`

is proved monotone nondecreasing in the physical parity weight `q`, so the minimum occurs at `J6=0`. The face reserve is strengthened from `9/50` to

`R >= 19/100`.

This is a face theorem only. It does not establish transverse monotonicity in positive `J4` or `J5`.

### FR3 — certified tube around the entire extreme face

Status: **PROOF NOTE PRESERVED; standalone replay script missing in this continuation tree**.

Artifact: `proofs/FR3_global_extreme_face_tube.md`.

A uniform face gap plus global covariance Lipschitz control yields a certified off-face tube for arbitrary `J3,J6>=0`. A simple corollary is

`J4,J5 <= 1/600000  =>  lambda2(M4)<=sigma_*`.

The small face box excluded by the uniform-gap cover is handed to the previously certified R7P-068 local cone after controlling parity drift.

### FR8 / FR9 — diagonal local staircase

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR8_FR9_local_staircase.md`
- `results/FR8_FR9_local_staircase.json`
- `tests/test_FR8_FR9_staircase.py`

FR8 (`r-v` arm):

- `|x| <= 1/6500`
- `u <= 1/8192`
- `v <= 1/4600`
- `e <= 1/100000`

FR9 (`r-u` arm):

- `|x| <= 1/6500`
- `u <= 1/2432`
- `v <= 1/8192`
- `e <= 1/100000`

The FR9 theorem deliberately uses `u<=1/2432`, although exploratory checking continued to pass to `1/2406` and failed at `1/2404`.

### FR10 — simultaneous `u-v` enlargement

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR10_diagonal_uv_local_box.md`
- `results/FR10_diagonal_uv_local_box.json`
- `tests/test_FR10_diagonal_uv.py`

Certified box:

- `|x| <= 1/8192`
- `u <= 1/4800`
- `v <= 1/4800`
- `e <= 1/3072`

This simultaneously enlarges both off-face directions and odd-parity mass relative to the original R7P-068 common-radius cube.

Negative controls retained by the same conservative checker:

- `e<=1/3052` fails;
- `u,v<=1/4700` at the certified `e` bound fails.

These are method failures, not ceiling violations.

### FR11 — low-odd-mass `r` arm

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR11_low_e_r_arm.md`
- `results/FR11_low_e_r_arm.json`
- `tests/test_FR11_low_e_r_arm.py`

Certified box:

- `|x| <= 1/7600`
- `u,v <= 1/4800`
- `e <= 1/81920`

Checker negative controls:

- `e<=1/75000` fails at these geometric radii;
- `|x|<=1/7500` fails at the certified `e` bound.

### FR12 — strong pure-`r` arm

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR12_pure_r_arm.md`
- `results/FR12_pure_r_arm.json`
- `tests/test_FR12_pure_r_arm.py`

Certified box:

- `|x| <= 1/4096`
- `u,v <= 1/8192`
- `e <= 1/1024`

This doubles the original signed `r` radius while allowing eight times the original odd-parity mass.

### FR13 — diagonal bulk box

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR13_diagonal_bulk_box.md`
- `results/FR13_diagonal_bulk_box.json`
- `tests/test_FR13_diagonal_bulk.py`

Certified box:

- `|x| <= 1/6200`
- `u <= 1/6800`
- `v <= 1/6800`
- `e <= 1/400`

This is a genuinely diagonal enlargement: all three geometric local directions are enlarged together, while odd-parity mass is more than twenty times the original R7P-068 radius.

The identical checker fails at `e<=1/399`; again this is only a proof-method negative control.

### FR14 — strengthened `r-v` arm

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR14_strengthened_rv_arm.md`
- `results/FR14_strengthened_rv_arm.json`
- `tests/test_FR14_strengthened_rv.py`

Certified box:

- `|x| <= 1/6400`
- `u <= 1/8192`
- `v <= 1/4600`
- `e <= 1/100000`

The same checker fails at `|x|<=1/6300` with the other bounds unchanged.

## 4. Intermediate continuation results that are NOT durable standalone theorems yet

The research session also produced FR4--FR7 intermediate results. They were used to guide later searches, but no separate proof/result/test package was frozen for them in the current tree. They must be replayed before promotion in a future report.

- **FR4:** hierarchical face-gap atlas around the boundary equality point.
- **FR5:** provisional global large-`J6` tail `exp(-2J6)<=1/5,000,000`, obtained from a boundary-gap atlas plus a covariance perturbation estimate.
- **FR6:** provisional `J4` anisotropic arm approximately `|x|<=1/8192`, `u<=1/2048`, `v<=1/8192`, `e<=1/30000`.
- **FR7:** provisional `J5` anisotropic arm approximately `|x|<=1/8192`, `u<=1/8192`, `v<=1/4096`, `e<=1/45000`.

Treat these as **session-certified / replay-required**, not as durable standalone claims. The later FR8--FR14 boxes are durable and independently replayed in this handoff.

## 5. What the continuation changed scientifically

The post-handoff picture is materially stronger than the 2026-09-14 final report:

1. The unresolved positive-orthant problem no longer includes the large-`J3`, large-`J4`, or moderately large-`J5` tails covered by FR1.
2. The extreme-face determinant reserve is stronger (`0.19`) and its worst parity endpoint is identified on the face.
3. A genuine off-face tube around the entire extreme face is available (FR3).
4. The R7P-068 equality neighborhood is not naturally a tiny isotropic cube. Multiple certified anisotropic boxes show a much larger, strongly direction-dependent safe region.
5. The local proof geometry is demonstrably anisotropic: expanding one coordinate can be safe while simultaneous expansion of another can make the conservative checker fail.
6. So far, checker failures have behaved as proof-method boundaries. They are not certified physical violations.

This continuation still does **not** prove that the entire remaining physical positive orthant satisfies the ceiling.

## 6. Current unresolved core

The durable unresolved set is the physical positive-orthant domain after removing the union of:

- all baseline boundary-Ising / extreme-face / intraparity certified regions;
- R7P-068 local cone;
- FR1 global `J3/J4/J5` tails;
- FR3 entire-face tube;
- FR8, FR9, FR10, FR11, FR12, FR13 and FR14 local anisotropic boxes.

FR4--FR7 may further reduce this set once replayed, but they must not be silently assumed in a rigorous global cover until their standalone artifacts are rebuilt.

No complete dependency-aware cover of this durable residual has yet been produced.

## 7. Numerical search status

During the continuation, adversarial searches repeatedly placed the best candidates on or infinitesimally beyond the artificial faces of the current certified-box union rather than at a stable interior violating point. The observed gaps remained negative (`lambda2-sigma_*<0`).

However, the raw frontier-search trajectories/seeds were not frozen as durable result files in the present continuation tree. Therefore this handoff records that behavior only as a **research-navigation observation**, not as evidence of a theorem or of absence of a counterexample.

A future search must save full fields, probabilities, seeds, objective values and the exact union-mask definition if it is to be used as promoted evidence.

## 8. Verification status

Baseline verification in the current continuation tree:

`PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src python verify.py`

Expected/current result:

- `VERIFY_PASS`
- 128 tasks
- 62 claims
- one declared historical `NONREPLAYED_INPUT` source ZIP.

Frontier sharded replay:

- `test_FR1_residual_tails.py`: 5 PASS
- `test_FR8_FR9_staircase.py`: 3 PASS
- `test_FR10_diagonal_uv.py`: 4 PASS
- `test_FR11_low_e_r_arm.py`: 3 PASS
- `test_FR12_pure_r_arm.py`: 2 PASS
- `test_FR13_diagonal_bulk.py`: 2 PASS
- `test_FR14_strengthened_rv.py`: 2 PASS

Total: **21 PASS, 0 FAIL** in the durable frontier test layer.

FR2/FR3 need standalone replay scripts if they are to be elevated to the same portability tier as FR1/FR8--FR14.

## 9. Exact next recommended atoms

### Priority 1 — replay/freeze FR5

Rebuild the boundary-gap atlas and perturbation calculation behind the provisional global `J6` tail. Completion criterion:

- a standalone proof/result/test artifact proving or correcting `exp(-2J6)<=1/5,000,000` globally in the physical four-amplitude model.

This is the most valuable missing durable piece because it would make the residual compact in all four field directions.

### Priority 2 — save a new adversarial residual search

Run a fixed-seed search on the complement of **only durable certified regions** first. Save:

- fields `(J3,J4,J5,J6)`;
- compact coordinates `(r,u,v,e)`;
- `lambda2-sigma_*`;
- seed and solver settings;
- exact region-mask reason for every excluded candidate.

Completion criterion: either a stable admissible positive-gap candidate or a reproducible negative-gap maximizer that identifies the next cover cell.

### Priority 3 — shifted-box proof just outside the local union

If the new maximizer again lies on an artificial wall of FR8--FR14, stop enlarging the equality-centered Taylor cube blindly. Use direct shifted intervals for the R7P-044 disjunction (`P<0` or `P1>0`) on a box centered at that wall point.

Completion criterion: one new certified shifted box with no dependence on a Taylor expansion centered exactly at the double root.

### Priority 4 — dependency-aware cover of the durable compact residual

After FR5 is replayed, return to the denominator-free determinant/inertia cover, now only on the finite residual not handled by tails/tubes/local boxes. Use exact physical parity coupling; do not replace it by independent `q in [1/2,1]`.

Completion criterion: zero unresolved leaves or one certified physical counterexample.

## 10. Do-not-repeat / guardrails

- Do not infer global 4D validity from numerical absence of a violation.
- Do not treat checker FAIL as a physical counterexample.
- Do not convexify the union of anisotropic local boxes unless convexity is separately proved.
- Do not assume monotonicity of `lambda2` in `q`; it was numerically falsified earlier in the continuation.
- Do not assume transverse monotonicity of the determinant reserve in `J4,J5`; that shortcut was also falsified.
- Do not enlarge the domain by freeing the physical parity relation and then interpret relaxed violations as physical.
- Do not transfer any C4 ceiling to the full seven-coordinate Hessian.
- Do not reopen the historical 128-task ledger; use a separate continuation ledger.

## 11. File map

Start here:

- `HANDOFF.md` — this continuation handoff.
- `HANDOFF_BASELINE_20260914.md` — frozen previous final handoff.
- `CONTINUATION_LEDGER.json` — machine-readable FR status.
- `results/FRONTIER_REPLAY_20260915.json` — replay summary.
- `logs/frontier_replay_20260915.log` — human-readable replay log.
- `proofs/FR*.md` — durable continuation proof notes.
- `results/FR*.json` — durable continuation results where available.
- `tests/test_FR*.py` — durable replay tests.
- `src/frontier_local_boxes.py` — generalized anisotropic R7P-068 checker.
- `src/frontier_residual_tails.py` — FR1 projected-tail checker.

The baseline report, claim ledger, certificates, source code and 2026-09-14 reproducibility material remain included unchanged underneath this continuation layer.


## 12. Post-handoff continuation update — FR5 / FR15 / FR16 / FR17

This section records work completed after the original 2026-09-15 handoff text above. It supersedes the earlier provisional FR5 status and the earlier exact-next-atom ordering.

### FR5 — global large-`J6` tail replayed and frozen

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Artifacts:

- `proofs/FR5_global_large_J6_tail.md`
- `results/FR5_global_large_J6_tail.json`
- `src/frontier_j6_tail.py`
- `tests/test_FR5_global_J6_tail.py`

The shifted boundary cover at `theta=sigma_*-10^-6` closes with 508 leaves and zero unresolved leaves: 356 `SAFE_A_NONPOS`, 150 `SAFE_B_NONNEG`, 2 `LOCAL_FR13`. Exact physical parity gives `e<=y/(1+y)` for `y=exp(-2J6)`. Finite feature-diameter enumeration gives `D^2<5`, hence the covariance perturbation is at most `5e`. Therefore `y<=1/5,000,000` implies perturbation `<10^-6`, while the two local leaves are covered by FR13.

**Consequence:** together with FR1, the durable unresolved region is bounded away from all four large-field tails; the remaining positive-orthant problem is compact in all four field directions.

### FR15 — strengthened FR9 parity arm

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Certified box:

- `|r-r_*| <= 1/6500`
- `u <= 1/2432`
- `v <= 1/8192`
- `e <= 1/14000`

The same checker fails at `e<=1/13000`; this is a proof-method negative control only. FR15 strictly subsumes FR9 in the parity direction and removes the first residual-search wall.

### FR16 — wall-adapted `r-v` box

Status: **INTERVAL_CERTIFIED_REPLAYED**.

Certified box:

- `|r-r_*| <= 1/5000`
- `u <= 1/8192`
- `v <= 1/5600`
- `e <= 1/100000`

This box was chosen from the post-FR15 optimizer wall and contains that wall point. The otherwise identical box with `|r-r_*|<=1/4800` fails the conservative checker.

### FR17 — frozen compact-residual adversarial search

Status: **NUMERICAL_SEARCH_ONLY**.

Artifacts:

- `src/frontier_residual_search.py`
- `results/FR17_residual_adversarial_search.json`

Four fixed-seed differential-evolution runs on the compact residual (after global FR1/FR5 tails and local FR8--FR16 masks) found no positive physical gap. The best saved point has negative `lambda2-sigma_*` and lies numerically on the FR16 `v=1/5600` wall to about `10^-12`. This is navigation evidence, not a theorem.

### Verification of the new durable layer

- `test_FR5_global_J6_tail.py`: 2 PASS
- `test_FR15_strengthened_parity_arm.py`: 1 PASS
- `test_FR16_wall_adapted_rv.py`: 1 PASS

Combined with the earlier 21 frontier tests, the durable frontier layer now has **25 PASS, 0 FAIL**.

### Updated next atom

The next scientifically useful step is no longer another blind equality-centered radius increase. The residual search tracks the boundary of the certified union. The preferred next step is a dependency-aware four-dimensional cover on the now compact residual, or a direct shifted-box/Weyl enclosure just outside the FR16 wall. The exact physical parity relation must remain coupled; do not relax `e` or the odd conditional law independently.
