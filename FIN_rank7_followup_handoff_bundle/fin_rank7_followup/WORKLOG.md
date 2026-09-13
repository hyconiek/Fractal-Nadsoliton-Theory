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
