# WORKLOG

## 2026-09-13T21:07:02.957048+00:00 — Campaign intake
- R7P-001 DONE: environment/provenance frozen. Current runtime has ~5.8 GiB RAM, 5 CPUs, no mounted Git repository.
- R7P-002 BLOCKED_DEPENDENCY: mandatory `fin_handoff_audit/` package not supplied/found. No regeneration from seeds attempted.
- R7P-003 DONE: 25 explicit claim atoms created; imported optimistic global H7 index claim corrected to REFUTED per master plan.
- R7P-004 DONE: missing-proof matrix created for phase census, boundary Ising, dominant mass, fold, global resolvent and full-7D counterexample replay.

## 2026-09-13T21:08:11.585604+00:00 — A infrastructure complete
- R7P-005 DONE: schema/immutable-input policy and negative fixtures.
- R7P-006 DONE: read-only verifier separated from result generator; mutation test rejects altered source.
- R7P-007 DONE: deterministic checkpoint/resume test passes.
- R7P-008 BLOCKED only by missing R7P-002 audit package.

## 2026-09-13 — Audit restored; packages B/C and stationary falsification
- R7P-002 unblocked: supplied audit files match historical hashes; fresh 19/19 rank-seven intake tests PASS; regenerated exact layer and manifest exactly match supplied results. 56 inherited regressions remain historical PASS only in this runtime.
- R7P-009..016 completed. New exact R7P-014 stationary inertia theorem: H7 and primal tangent Hessian have identical negative/null counts at common interior stationary states.
- R7P-017..021 completed; unique corrected sech-resolvent minimum interval-certified. R7P-023 fixed parity shifts certified negative; reoptimized slope not promoted.
- R7P-034/035 exact two-harmonic reduction + numerical atlas. At g=5 found full stationary candidates with H7 indices 2 and 3.
- R7P-036 COUNTEREXAMPLE_CERTIFIED: exact g=5 root box around J=1.24494819516,K=0.77955562099 passes uniform interval Krawczyk; H7 decomposes into positive 3sin and 36 blocks plus two disjoint indefinite 45 blocks, proving exact index 2, nullity 0. Universal stationary-point index<=1 claim is false without gain restriction.
- Same certified root is H4 index 1 but H7 index 2; extra negative direction is sine 45 sector. R7P-039 partial mismatch row saved.

## R7P-041--045 boundary-Ising algebra

- Rebuilt the exact four-state even-parity boundary with multiplicities `(1,2,1,2)` and retained the mandatory `-log(2)/2` degeneracy shift.
- Proved the three polynomial inequalities are exact on the positive interior but only an outer relaxation on the zero-probability boundary.  The true exponential-family closure has only the `{1}`, `{1,2}`, and `{1,3}` infinite-parameter supports with explicit ratio restrictions.
- Independently derived the covariance characteristic invariants and the determinant factor `3 lambda3 lambda4 lambda5 p1p2p3p4/8`.
- Strengthened the threshold test: after proving `P''(sigma*)/2>0` globally via an exact four-point enclosing ball, `lambda2>sigma*` is equivalent to `P(sigma*)>0` and `P'(sigma*)<0`.
- Corrected an earlier symbolic enclosing-radius formula while preserving its prior numerical value: the exact `R^2` is `(16 l3 l4+9 l4^2+10 l4 l5+l5^2)/(96 l4)`.
- Proved the special double-root identities exactly, including the remaining eigenvalue `l4 l5/(24 l3)` and strict gap below `sigma*`.

## 2026-09-14 — R7P-046--056 boundary-Ising theorem
- Corrected archive inspection: the full `fin_rank7_followup/` source tree is present; one separately hashed immutable source ZIP expected by top-level `verify.py` is missing from the supplied outer archive and remains an explicit packaging defect.
- R7P-046 exact compactification: `p ∝ (1,2st^3,rt^4,2rst)` on `[0,1]^3`; all zero-probability strata classified and rank<=1.
- R7P-047: three rational one-constraint relaxation witnesses are interval-certified to have `lambda2>sigma_*`.
- R7P-048--050: frozen exact A/B proof specification plus interval Bernstein engine and reason-code classifier.
- R7P-051 pilot: no remote unresolved components; four equality-neighborhood boxes retained honestly.
- R7P-052: exact zero value/gradient at the double root plus strict interval tangent-cone Hessian Schur certificate on `[51/128,103/256]×[255/256,1]^2`.
- R7P-053 refined cover: 436 leaves, 0 unresolved.
- R7P-054 independent checker: CHECK_PASS; altered spectrum, split corruption, leaf deletion and invalid reason mutations rejected.
- R7P-055 INTERVAL_CERTIFIED global boundary theorem: `lambda2<=sigma_*` on the exact four-state boundary-Ising closure, equality only at R7P-045 double root. No off-face/full-7D promotion.
- R7P-056 dependency propagation: G paid; H is next.

## 2026-09-14 — R7P-046..056 global boundary-Ising closure
- R7P-046: exact compactification to `[0,1]^3` with weights `(1,2 s t^3,r t^4,2 r s t)`; only zero-probability supports `{1,2}`, `{1,3}`, `{1}` and each has covariance rank <=1.
- R7P-047: interval-certified rational negative controls show each physical interior inequality is individually essential; these are relaxed-domain warnings, not physical counterexamples.
- R7P-048..050: froze exact polynomial criterion and implemented exact rational interval Bernstein/de Casteljau cover rules with adversarial classifier tests.
- R7P-051: pilot cover left only the known equality neighborhood ambiguous.
- R7P-052: exact `A=grad A=0` at the double root plus strict interval tangent-cone Hessian Schur signs prove `A<0` off the root inside `[51/128,103/256] x [255/256,1]^2`.
- R7P-053: refined cover completed with 436 leaves and zero unresolved leaves: 285 SAFE_A, 150 SAFE_B, one local equality certificate.
- R7P-054: independent checker PASS. Mutation tests reject removed leaf, corrupt split, flipped safe reason/sign, and altered spectral specification. The checker was optimized only by deterministic preflight/memoization; the valid tree still receives full polynomial replay.
- R7P-055: **INTERVAL_CERTIFIED GLOBAL BOUNDARY THEOREM** — throughout the exact boundary-Ising exponential-family closure, `lambda2<=sigma_*`; equality only at the exact R7P-045 double root. Relative to the original four-amplitude chart, the boundary model is the `J6->+infinity` parity-saturation limit.
- R7P-056: propagated dependencies. H can now use the boundary theorem, but finite-J6/off-face four-amplitude curvature remains open.
- Reproducibility caveat preserved: the historical immutable zip named by `INPUT_HASHES.json` is missing from the supplied handoff bundle; expected hashes were not rewritten to hide this packaging defect.

## 2026-09-14 — R7P-057 exact intraparity weight theorem
- Enumerated even/odd trigonometric features directly from the shared fields.
- Derived exact even/odd partition sums and proved `q0>=1/2` by monotonicity in `J3` plus a Taylor series whose coefficients are zero at orders 2 and 4 and strictly positive thereafter.
- Equality: `J3=J5=0`, arbitrary nonnegative `J4`; extension to `J6>=0` follows from the logistic parity factor.
- 3 new regression tests PASS, including a negative-J6 domain counterexample.

## 2026-09-14 — R7P-058..064 intraparity closure
- R7P-058: exact physical `(u,d)` odd-sector domain and covariance reconstructed. Strict spectral signs prove `sup lambda1(C_-)=(3lambda4+lambda5)/32`, attained only at an infinite-field closure point.
- R7P-059: `lambda_min(C_-)<sigma_*` globally; dangerous `lambda1(C_-)>=sigma_*` reduced exactly to the stated `u_±` interval and `d^2>=dmin^2(u)` with denominator signs certified.
- R7P-060: same-field dominant even mass reduces exactly to `J3=0`, `d=dmin(u)`. Numerical 1D minimum reproduces `u≈0.5280356754`, `d≈0.4522658586`, `p_dom≈0.719071046876`; locator remains numerical.
- R7P-061: exact-rational Bernstein covers certify `d/u>=0.8515`, `Y^2>=11.45`; independently recomputed rational consequence `p_dom>711/1000` (old imported `0.7112098557` not used as proof input).
- R7P-062: all dominant-vertex branches of the tetrahedral `e2` envelope reduced and interval-covered; `p_dom>=711/1000` implies `lambda2(C_+)<511/2000`.
- **R7P-063 INTERVAL_CERTIFIED:** two-case Weyl argument closes `lambda2(W_par)<=sigma_*` for all shared nonnegative fields. Final certified `2sigma - Cminus_sup - 511/2000` lower margin is about `0.00134546`.
- R7P-064: regression guards record failure of `q>=1/2` for negative J6 and rejection of arbitrary relaxed Cplus simplex laws. No theorem is inferred from a bounded negative search over decoupled physical sectors.

## 2026-09-14 — Work package I continuation (R7P-065--068 partial)
- R7P-065 DONE / EXACT_PROVED: froze the full four-amplitude target and the stable Schur/inertia representation. Exact strict intervals give `eta>=1-lambda6/(12 sigma_*)>0`; no full resolvent inverse is used.
- R7P-066 DONE / NUMERICAL_REPRODUCED: three fixed-seed DE runs over the exact compactified closure plus structured face/finite-J6/local slices found no positive gap beyond ~2.8e-16 roundoff at the known equality point. This is not a proof.
- R7P-067 DONE / EXACT_PROVED: exact compactification of the entire unbounded nonnegative field orthant. Boundary-aligned chart `(r,s,t,y)` has even weights `(1,2st^3,rt^4,2rst)` and odd weights `2 sqrt(r) y (s t^(2+sqrt3), t^2, s t^(2-sqrt3))`; `y=0` exactly recovers the certified G boundary cube.
- Upgraded the previously numerical reoptimized R7P-023 envelope slope to an interval-certified range around `-0.1312828584000`.
- R7P-068 IN_PROGRESS / PARTIAL_ANALYTIC: no first-order physical tangent through the equality point creates an immediate violation; acceptance still requires an explicit finite radius with validated mixed higher-order remainder.
