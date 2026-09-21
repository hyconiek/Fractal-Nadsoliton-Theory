# R7N result tables

## Primary claims

| Claim | Final status | Certified result | Remaining obligation |
|---|---|---|---|
| Target P: `lambda2(M4)<=67/250` globally | **UNRESOLVED** | Current proof tree certifies `63.68950527593205%` of the compact hull by volume, plus accepted unbounded tails; 5,432 compact residual cells remain | New physical-coupled matrix enclosure; no third identical global `t` split |
| Target S: `lambda2(M4)<=sigma` globally | **OPEN / not re-entered** | Accepted baseline sigma-safe regions retained | A sharp-target strategy after/without Target-P closure |
| Quartic fixed-fixture phase census | **INTERVAL CERTIFIED — EXACT 60** | Full torus exhausted: 27,272 gradient-exclusion leaves + 640 root-collar leaves, zero unresolved | Optional amplitude stability only |
| Full log-mgf fixed-fixture phase census | **INTERVAL CERTIFIED — EXACT 60, independently replayed** | Full torus exhausted by K16/K20: 79,633 formula-replayed gradient leaves + 864 root-collar leaves, zero unresolved | Optional amplitude stability / homotopy only |

## Phase counts

| Fixture | Critical points | negative index 0 | index 1 | index 2 | index 3 | Exhaustive? |
|---|---:|---:|---:|---:|---:|---|
| Quartic K4 | 60 | 12 | 24 | 18 | 6 | Yes |
| Full log-mgf | 60 | 12 | 24 | 18 | 6 | Yes |

The two catalogs have a labelled one-to-one endpoint correspondence. No global `K_alpha` continuation theorem is claimed.

## Quartic complement proof

| Quantity | Value |
|---|---:|
| Processed cells | 55,760 |
| Gradient-exclusion terminal leaves | 27,272 |
| Certified root-collar leaves | 640 |
| Unresolved leaves | 0 |
| Certified roots | 60 |

## Full complement proof

| Layer | Safe gradient leaves | Root-collar leaves | Residual passed onward |
|---|---:|---:|---:|
| K16 baseline | 54,341 | 0 | 5,382 |
| K20 direct reclassification | 2,887 | 0 | 2,495 |
| K20 adaptive closure | 22,405 | 864 | 0 |
| **Total** | **79,633** | **864** | **0** |

K16 uses 25 retained resonances with uniform full-gradient error bound about `2.383381984359391e-8`. K20 uses 45 retained resonances and reduces the rigorous bound to about `1.9896999978556874e-10`.

## Full root collars

| Quantity | Value |
|---|---:|
| Locally certified full roots | 60 |
| Certified collar radii | 0.0003–0.0015 rad |
| Minimum distinct-center torus `L_inf` separation | ~0.7131676442373065 rad |
| Minimum separation after subtracting both collar radii | ~0.7123676442373066 rad |

## Target P bounded campaign

| Stage | Residual fraction of compact hull |
|---|---:|
| Before directed `t` refinement | ~54.4593% |
| After first directed `t` pass | ~41.8110% |
| After second directed `t` pass | **~36.3105%** |

The largest numerical `lambda2` at the centers of the final residual cells was ~`0.25392497687445353`, below `0.268`. This is a navigation diagnostic only, not a certified bound on those cells.

## Verification

| Check | Result |
|---|---|
| Fresh baseline replay | 119/119 tests, 29/29 files |
| Integration controls | 7/7 PASS |
| K16 formula-level replay | 54,341/54,341 PASS |
| K20 formula-level replay | 25,292/25,292 PASS |
| Phase geometry/mutation audit | PASS |
| Campaign promotion/firewall mutation audit | PASS |
| Consolidated intake source hashes | 23/23 match |

## Secondary full-X7 / energy lane

No expensive secondary campaign was launched. The successful K16/K20 tools are fixed-amplitude phase-torus tools and do not provide a licensed seven-dimensional exclusion method. No new full-X7 stationary exhaustion, global minimizer result, or improved energy lower bound is claimed.
