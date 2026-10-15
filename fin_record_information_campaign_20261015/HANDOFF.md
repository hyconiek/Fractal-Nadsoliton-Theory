# FIN RI-001–RI-005 — new research only, architect intake handoff

**Date:** 2026-10-15. **Execution directory:** `fin_record_information_campaign_20261015/`. **Package scope:** only RI work; **no earlier archive or full predecessor report is included**. The sources are referenced by SHA-256 in `SOURCE_HASHES.json`. An original repository intake review file is not supplied here as a mounted source; its content in the user message is a declared dependency, not a reconstructed hash.

## Terminal disposition

| Card | Terminal result | Scientific status |
|---|---|---|
| RI-001 | Exact tagged/count quotient and five distinct information contracts | **EXACT FINITE IDENTITIES** |
| RI-002 | 27 frozen joint cells and finite-record Bayes risk converge; strict qualitative advantage over 1/2 | **CONDITIONAL THEOREM; NUMERICS UNCERTIFIED** |
| RI-003 | Limiting stationary nonlinear causal filter, innovation and unbounded-observation moment payment | **LIMIT-MODEL ANALYTIC DERIVATION** |
| RI-004 | Exact finite parity-history point-process likelihood; compare full record to coarse limit | **LIKELIHOOD EXACT / FULL SUFFICIENCY OPEN** |
| RI-005 | Scope and single architect checkpoint | **HANDOFF ONLY** |

## Frozen nine-record experiment

`lambda=0`, `k=1`, `S=1`, two observations of the exact integrated tag history on scaled intervals `[0,1/4]`, `[1/4,1/2]`, quantizer LOW<=-1/2, MIDDLE<=+1/2, HIGH>+1/2; predict three-outcome sign of the collective mode at scaled time1 (MINUS/ZERO/PLUS). All 27 finite-N cells retained. The incoming asymptotic OBS result implies joint-TV convergence **of this fixed finite alphabet**, and 0–1 Bayes risk convergence. A Gaussian-association argument proves strict limiting risk below 1/2, without fitting bins or depending on numeric grid estimates.

`JOINT_TABLE_OR_BOUNDS.json` holds exact integral descriptions of all 27 cells alongside explicitly un-certified deterministic approximations. Grid risk values 0.41355, 0.41508 and 0.41556 are **not certified true decimals** or experimental claims; grid-difference is not a valid error enclosure. No useful finite-N N0 is given.

## Limiting causal filter

`d pi(f) = pi(Gf) ds + (k/2) Cov_pi(f,z) [dY-pi(z)ds]`, `Gf=(chi*lambda*z-z^3/3)f'+f''`. For finite horizons, stationary quartic tails satisfy Novikov's quadratic exponential moment via Jensen; the filter is the normalized Kallianpur–Striebel conditional law and its normalized innovation is Brownian in its observation filtration. It is **not** the finite-N microscopic filter and no Kalman closure is claimed.

## Exact open gate

One-tag exact microscopic integral path determines its parity jump history by its right derivative, whereas its limiting Brownian record has no such derivative. This prevents full-record sufficiency from following from weak convergence. `RI-004_FULL_RECORD_GATE.md` defines the actual finite hidden tagged/count process likelihood for observed switch times and the missing `RI-LIKELIHOOD-UI` comparison theorem. No claim of asymptotic full-record equivalence or independent-information conservation.

## Integrity and audit

`independent_checks.py` reconstructs small exact probability and labelled-generator algebra separately from `ri_exact.py`; **5626/5626** finite rational assertions passed in the producing session. These tests are not formalized stochastic theorems or physical validation. See `REPLAY.md`, `NUMERICAL_SCOPE.md`, individual task folders, the standalone `MANIFEST.sha256` and `verify_package.py`.

**Architect return:** do not automatically launch more campaigns; decide only whether `RI-LIKELIHOOD-UI` merits a bounded proof attempt. No modifications to AGENTS, past sources, physical inputs or the empirical lane.
