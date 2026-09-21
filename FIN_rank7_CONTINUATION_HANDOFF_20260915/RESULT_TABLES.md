# FIN rank-seven campaign — final result tables

Date: 2026-09-14. Values below are copied from saved campaign evidence; the proof-status column controls their permitted interpretation.

## Local radial/stationary events

| Result | Value / enclosure | Proof status | Evidence |
|---|---:|---|---|
| Four-amplitude equal-energy event | g around 3.71834489812038 | interval/Krawczyk local event | `certificates/R7P-026_equal_energy.json` |
| Simple fold | g in the saved R7P-031 certified interval around 3.51564471684 | interval-certified local fold | `certificates/R7P-031_fold.json` |
| Stationary index-2 counterexample | g=5 | certified stationary full-H7 index >=2 | `proofs/R7P-036_stationary_index2_counterexample.md` |
| Full-7D discovery atlas | orbit counts 3,3,3,15 at g=3.7, 3.7183449, 4, 5 | numerical saturation, not exhaustion | `results/R7P-037_full7_stationary_atlas.json` |

## Curvature / covariance domains

| Result | Domain | Proof status |
|---|---|---|
| Boundary-Ising ceiling | exact attainable boundary-Ising closure | exact/interval global boundary theorem; equality only at the certified double-root boundary point |
| Intraparity ceiling | shared J3,J4,J5,J6 >= 0 | certified `lambda2(W_par) <= sigma_*` |
| Extreme face | J4=J5=0, J3,J6>=0 | certified `lambda2(M4)<=sigma_*`; sharp only at compactified equality point |
| Local off-face cone | tangent box radius rho=1/8192 | interval-certified; checker rejects rho=1/4096 |
| Large-J5 tail | t=exp(-J5/2) <= 2^-11 | interval-certified by trace bound |
| Remaining positive-orthant 4D core | complement of the certified subdomains above | unresolved; no global ceiling claimed |

## Phase results

| Result | Saved result | Proof status |
|---|---:|---|
| Quartic candidate roots | 60 | all 60 locally isolated; complement not exhausted |
| Quartic Morse counts | 6 maxima, 42 saddles, 12 minima | locally certified at root boxes |
| Quartic complement cover | 1272 unresolved leaves after bounded pass | explicit resource-stop / unresolved map |
| Full roots continued/isolated near quartic roots | 60 | local one-to-one boxes; no theorem excluding extra full roots elsewhere |
| First pure-k6 angular instability | r in [0.41421132290, 0.41421132293] | interval-certified, k3 sector first |
| Full angular fold | r ~= 0.3463027188406942 | reconstructed; local proof level recorded in R7P-096 |
| Full angular coexistence | r ~= 0.3645555701282058 | locally interval-certified event |
| Angular saddle barrier | ~=1.2902908116e-5 | reconstructed full model |
| Quartic fold | ~=0.346613544874 | reconstructed; disagrees with imported 0.34660000115 and is preserved as a provenance discrepancy |

## Global rank-seven results

| Result | Value | Proof status / limitation |
|---|---:|---|
| First energetic-change bracket | `g_global in [2.8934, 3.71835]` | rigorous lower transfer via A7<=A_full + rational upper witness; attaining orbit not identified |
| g=4 best discovered stationary orbit | localized, Phi ~= -0.1388476630 | numerical atlas; global complement unresolved |
| g=4 saddle | Phi ~= 0.0225772250 | numerical atlas |
| g=4 declared Euclidean gradient flow | saddle branches numerically approach uniform / localized orbit | numerical trajectory geometry only; no validated heteroclinic tube or physical law |

## Passive finite-N / quantum comparison

| Result | Value | Status |
|---|---:|---|
| Hodge ranks | 11 tree + 55 cycle | structural theorem/checker |
| Distinct positive cycle eigenvalue clusters | 28 | numerical distinct-count only |
| Equilibrium mean-zero mode variance | 1/(12N) | exact |
| Alternating-mode excess kurtosis | -2/N | exact finite-N diagnostic |
| Visible passive-memory degree | 5 (residue ranks 1+2+2+0) | structural plus numerical pole-residue decomposition |
| Quantum loading robustness | at least (49/50)^2 times accepted mixed-family discord floor for gamma in [0.049,0.05] | conditional scoped theorem; historical gap certificate provenance retained |
| Quantum-to-localization causal bridge | none established | explicit nonidentifiability / missing-map result |
