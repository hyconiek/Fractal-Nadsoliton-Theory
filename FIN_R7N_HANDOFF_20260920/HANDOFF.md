# FIN R7N campaign handoff

## 0. Identity and state

- Date: 2026-09-20. Campaign: R7N-001--R7N-064.
- Baseline: FR223 continuation package plus the September 19 consolidated intake review. Imported source archives were kept read-only; R7N work lives in `fin_rank7_next_campaign/`.
- Environment: Python 3.13.5; NumPy 2.3.5; SciPy 1.17.0; SymPy 1.14.0; mpmath 1.3.0.
- The current package contains source manifests, exact/interval proof objects, cover trees, replay checkpoints, logs and result ledgers.
- No background worker is intentionally left running.

## 1. Executive scientific change

The major new result is complete phase exhaustion at the exact declared fixed fixture.

- Quartic phase function: **exactly 60 critical points** on the full phase 3-torus, negative-index histogram `12,24,18,6`.
- Full log-mgf phase function: **exactly 60 critical points** on the full phase 3-torus, same histogram; the complete proof is independently replayed at geometry, mutation and formula levels.
- Target P at `tau0=67/250` remains globally open. The current proof tree certifies about `63.68950527593205%` of the compact hull by volume and leaves `36.31049472406795%` unresolved.
- Target S at `sigma` remains open and was not re-entered after the Target-P stop.
- No new full-X7/global-energy theorem was produced.

Baseline reproduction is not counted as new science.

## 2. Target P: rational practical ceiling

Threshold: exact `tau0=67/250`. Domain: all four nonnegative shared fields, decomposed into accepted unbounded tails plus the compact `(r,s,t,y)` hull.

Status: **UNRESOLVED PARTIAL CERTIFICATE**. The final compact residual occupies `0.3631049472406795` of compact-hull volume and contains 5,432 terminal residual cells. No admissible counterexample was certified. The largest sampled residual-center `lambda2` is about `0.25392497687445353`, but center values are navigation only and are not used as proof.

The first directed `t` pass reduced the incoming residual by about 23.23%; the second by about 13.16%. A third identical global split is prohibited without a new enclosure primitive.

Conditional gain consequence: if Target P is later proved globally, `H4=I4/g-M4` has at most one strictly negative eigenvalue for supplied `0<g<=250/67`, with the endpoint caveat stated in the plan. This campaign does not satisfy the antecedent globally.

## 3. Target S: sharp sigma ceiling

Status: **OPEN**. No tau0-only leaf is promoted to sigma. No sharp equality classification, physical gap or global sharp cover is claimed. The sharp lane was stopped after Target P retained a large residual and the primary phase lane offered higher proof value.

## 4. Phase results

### Exact fixture and coordinates

Fixed amplitudes: `r3=0.1131879146`, `r4=0.1698528641`, `r5=0.2269339093`, `z6=-0.3380663037`. Phase coordinates are normalized by `z=phi/(2*pi)` on the full periodic `[0,1]^3` torus.

### Quartic

- 60 local roots recertified.
- Certified uniqueness collars of radius `0.05` rad.
- Corrected complement cover: 27,272 gradient-exclusion leaves + 640 root-collar leaves; zero unresolved.
- Fixed-sign symmetry subgroup has 12 elements: even translations and reflections. Odd translations flip the alternating sign and are not licensed within this fixture.
- Exact total: **60**, with 9 subgroup orbits and index histogram `12,24,18,6`.

### Full log-mgf

- 60 local roots recertified; certified collar radii `0.0003`--`0.0015` rad.
- Minimum pairwise root-center torus `L_inf` distance ~`0.7131676442373065` rad; minimum collar-disjointness margin ~`0.7123676442373066` rad.
- K16: 25 resonances, rigorous gradient error ~`2.383381984359391e-8`; 54,341 safe leaves, 5,382 passed to K20.
- K20: 45 resonances, rigorous gradient error ~`1.9896999978556874e-10`; 2,887 direct safe leaves plus adaptive closure with 22,405 safe and 864 root-collar leaves; zero unresolved.
- Formula replay: K16 54,341/54,341, K20 25,292/25,292, zero failures.
- Exact total: **60**, index histogram `12,24,18,6`.

Quartic/full correspondence is only a labelled endpoint bijection. No global `K_alpha` homotopy theorem is claimed. Amplitude robustness was optional and not attempted.

## 5. Optional full-X7 / global-energy results

No secondary theorem was launched. The new phase tools do not justify restricting generic X7 states to the fixed phase/amplitude fixture, and Target P remains partial. The known g=5 index-two stationary witness remains a regression guard. No improved energy lower endpoint is claimed.

## 6. Complete R7N-001--064 task ledger

| ID | Execution status | Scientific status | Evidence | Remaining atom |
|---|---|---|---|---|
| R7N-001 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-001_intake_manifest.json, INPUT_MANIFEST.sha256, environment.json | — |
| R7N-002 | DONE | NUMERICAL_REPRODUCED | results/R7N-002_fresh_119_replay.json, logs/R7N-002_integration_7_compat.log | — |
| R7N-003 | DONE | INTERVAL_CERTIFIED | results/safe_union_v2_audited.json, results/R7N-003_canonical_vs_alternate_repairs.json | — |
| R7N-004 | DONE | INTERVAL_CERTIFIED | results/R7N-004_coverage_tests.json | — |
| R7N-005 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-005_readonly_verification.json | — |
| R7N-006 | DONE | INTERVAL_CERTIFIED | results/R7N-006_arithmetic_benchmark.json, results/R7N-006_arithmetic_policy.json | — |
| R7N-007 | DONE | NOT_A_SCIENTIFIC_CLAIM | TASKS.json, CLAIMS.json, STATE_MAP.md | — |
| R7N-008 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-008_first_batch_decision.json | — |
| R7N-009 | DONE | NUMERICAL_REPRODUCED | results/R7N-009_compact_formula_diagnostics.json | — |
| R7N-010 | DONE | CONDITIONAL_LEMMA | proofs/R7N-010_residual_hull.md | — |
| R7N-011 | DONE | CONDITIONAL_LEMMA | results/R7N-011_physical_coupling.json | — |
| R7N-012 | DONE | INTERVAL_CERTIFIED | results/R7N-012_schur_identity_diagnostics.json, results/R7N-012_014_generic_threshold_tests.json | — |
| R7N-013 | DONE | INTERVAL_CERTIFIED | candidates/R7N-013_compression_candidate.json, results/R7N-013_compression_probe100_v2.json, src/compression_combined.py | — |
| R7N-014 | DONE | INTERVAL_CERTIFIED | results/R7N-012_014_generic_threshold_tests.json, src/target_p_e2_bound.py, src/compression_centered_moment.py | — |
| R7N-015 | DONE | INTERVAL_CERTIFIED | results/R7N-015_016_cover_infrastructure.json, src/target_p_trace_cover.py, results/safe_union_v2_audited.json | — |
| R7N-016 | DONE | NOT_A_SCIENTIFIC_CLAIM | proofs/global_cover_format.json, results/R7N-015_016_cover_infrastructure.json | — |
| R7N-017 | DONE | INTERVAL_CERTIFIED | results/R7N-017_target_p_implication.json | — |
| R7N-018 | DONE | NUMERICAL_NEW | results/R7N-018_target_p_navigation.json | — |
| R7N-019 | DONE | UNRESOLVED | checkpoints/R7N-019_trace_cover.json, checkpoints/R7N-020_tight_rounded_trace_only.json, results/R7N-019_024_practical_lane_final.json | Residual compact-hull cells remain; no global PASS. |
| R7N-020 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-021_t_refine_analysis.json, results/R7N-021_second_t_analysis.json, results/R7N-019_024_practical_lane_final.json | — |
| R7N-021 | DONE | UNRESOLVED | checkpoints/R7N-021_t_refine_once_v2.json, checkpoints/R7N-021_t_refine_second_cheap_v1.json, results/R7N-021_t_refine_analysis.json, results/R7N-021_second_t_analysis.json | 36.31049472406795% of compact hull volume remains uncertified; new physical-coupled enclosure needed. |
| R7N-022 | BLOCKED_DEPENDENCY | NOT_A_SCIENTIFIC_CLAIM | — | No complete Target-P cover exists; a full global-PASS checker is therefore inapplicable. Accepted local leaves are checked by their component checkers. |
| R7N-023 | DONE | INTERVAL_CERTIFIED | results/R7N-019_024_practical_lane_final.json, results/R7N-017_target_p_implication.json | Target P remains globally open. |
| R7N-024 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-019_024_practical_lane_final.json | — |
| R7N-025 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Sharp-sigma lane was not re-entered after Target P retained a large certified residual; no tau0-only evidence is transferred to sigma. |
| R7N-026 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Sharp-sigma lane was not re-entered after Target P retained a large certified residual; no tau0-only evidence is transferred to sigma. |
| R7N-027 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Sharp-sigma lane was not re-entered after Target P retained a large certified residual; no tau0-only evidence is transferred to sigma. |
| R7N-028 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Sharp-sigma lane was not re-entered after Target P retained a large certified residual; no tau0-only evidence is transferred to sigma. |
| R7N-029 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Sharp-sigma lane was not re-entered after Target P retained a large certified residual; no tau0-only evidence is transferred to sigma. |
| R7N-030 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Sharp-sigma lane was not re-entered after Target P retained a large certified residual; no tau0-only evidence is transferred to sigma. |
| R7N-031 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Sharp-sigma lane was not re-entered after Target P retained a large certified residual; no tau0-only evidence is transferred to sigma. |
| R7N-032 | DONE | UNRESOLVED | results/R7N-032_four_amplitude_report.md, results/R7N-019_024_practical_lane_final.json | Target P and Target S both remain open globally. |
| R7N-033 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-033_phase_cover_audit.json | — |
| R7N-034 | DONE | CONDITIONAL_LEMMA | proofs/R7N-034_normalized_phase_contract.json | — |
| R7N-035 | DONE | INTERVAL_CERTIFIED | results/R7N-035_quartic_uniqueness_collars.json, results/R7N-035_041_phase_local_recertification.json | — |
| R7N-036 | DONE | INTERVAL_CERTIFIED | results/R7N-036_revalidated_old_safe_terms.json | — |
| R7N-037 | DONE | INTERVAL_CERTIFIED | checkpoints/R7N-037_quartic_cover_normalized.json | — |
| R7N-038 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | — |
| R7N-039 | DONE | INTERVAL_CERTIFIED | results/R7N-039_symmetry_and_count.json | — |
| R7N-040 | DONE | INTERVAL_CERTIFIED | results/R7N-040_quartic_phase_theorem.json, results/R7N-040_quartic_phase_theorem.md | — |
| R7N-041 | DONE | INTERVAL_CERTIFIED | results/R7N-041_full_uniqueness_collars.json | — |
| R7N-042 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-043_K16_surrogate.json, results/R7N-043_K20_surrogate.json | — |
| R7N-043 | DONE | INTERVAL_CERTIFIED | results/R7N-043_K16_surrogate.json, results/R7N-043_K20_surrogate.json, results/R7N-043_K16_numeric_validation.json | — |
| R7N-044 | DONE | INTERVAL_CERTIFIED | results/R7N-044_full_phase_theorem.md, results/R7N-044_full_census_candidate.json, checkpoints/R7N-044_full_cover_k16.json, checkpoints/R7N-044_K20_residual.json, checkpoints/R7N-044_K20_adaptive.json | — |
| R7N-045 | DONE | NUMERICAL_REPRODUCED | results/R7N-045_endpoint_matching.json | — |
| R7N-046 | DONE | INTERVAL_CERTIFIED | results/R7N-046_phase_exhaustion_final.json, results/R7N-046_phase_exhaustion_audit.json, checkpoints/R7N-046_K16_formula_replay.json, checkpoints/R7N-046_K20_formula_replay.json | — |
| R7N-047 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | — | Optional amplitude-stability theorem not attempted; fixed-fixture result is complete and no uniform amplitude-box margins were paid. |
| R7N-048 | DONE | INTERVAL_CERTIFIED | results/R7N-048_phase_lane_summary.json, results/R7N-046_phase_exhaustion_final.json | — |
| R7N-049 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-049_056_secondary_stop.json | — |
| R7N-050 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-049_056_secondary_stop.json | — |
| R7N-051 | BLOCKED_DEPENDENCY | NOT_A_SCIENTIFIC_CLAIM | results/R7N-049_056_secondary_stop.json | — |
| R7N-052 | BLOCKED_DEPENDENCY | NOT_A_SCIENTIFIC_CLAIM | results/R7N-049_056_secondary_stop.json | — |
| R7N-053 | BLOCKED_DEPENDENCY | UNRESOLVED | results/R7N-049_056_secondary_stop.json | — |
| R7N-054 | OUT_OF_SCOPE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-049_056_secondary_stop.json | — |
| R7N-055 | BLOCKED_DEPENDENCY | NOT_A_SCIENTIFIC_CLAIM | results/R7N-049_056_secondary_stop.json | — |
| R7N-056 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-049_056_secondary_stop.json | — |
| R7N-057 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-057_traceability.json, THEOREM_REGISTER.json | — |
| R7N-058 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-058_mutation_audit.json, results/R7N-046_phase_exhaustion_audit.json | — |
| R7N-059 | DONE | NUMERICAL_REPRODUCED | verification.json, results/R7N-059_verification.json, logs/R7N-059_integration_fresh.log | — |
| R7N-060 | DONE | NOT_A_SCIENTIFIC_CLAIM | RESULT_TABLES.md | — |
| R7N-061 | DONE | NOT_A_SCIENTIFIC_CLAIM | REPORT.md, THEOREM_REGISTER.json | — |
| R7N-062 | DONE | NOT_A_SCIENTIFIC_CLAIM | AGENTS_PATCH.md, ACCEPTED_RESULTS_R7N.md | — |
| R7N-063 | DONE | NOT_A_SCIENTIFIC_CLAIM | results/R7N-063_bundle_replay.json, portable_verify.py, logs/R7N-063_clean_replay.log | — |
| R7N-064 | DONE | NOT_A_SCIENTIFIC_CLAIM | HANDOFF.md, verification.json | — |

## 7. Theorem register

The machine-readable register is `THEOREM_REGISTER.json`. The accepted theorem-level entries are:

1. `QUARTIC_PHASE_EXACT_60` — interval-certified exact fixed-fixture root census.
2. `FULL_PHASE_EXACT_60` — interval-certified exact fixed-fixture root census with independent geometry/mutation/formula replay.
3. `TARGET_P_PARTIAL` — explicitly partial-domain result only.
4. `TARGET_P_GAIN_IMPLICATION_CONDITIONAL` — exact/interval-certified conditional implication; antecedent not globally proved.

Each entry states assumptions, proof files, replay commands and nonconclusions.

## 8. Counterexamples and failed methods

- No Target-P counterexample was certified. Conservative interval-checker failure is never treated as one.
- A single fixed global hyperplane/compression candidate fails as a global Target-P proof; adversarial navigation found a restricted value around `0.285996 > 0.268` near the large-J6 boundary.
- Direct full-gradient interval evaluation was too dependency-inflated on microscopic annuli around phase roots.
- K16 alone hit the 200,000-cell budget without exhausting the thin annuli. The campaign pivoted to K20, whose ~119-fold smaller uniform gradient error closed the residual.
- The physical-coupled Schur prototype remained too slow and interval-inflated on wide cells for a global Target-P run.
- Independent weight-box/Bernstein relaxation lost too much shared `(r,s,t,y)` dependence and did not close the Target-P residual.

## 9. Verification and portability

`verification.json` records:

- 119/119 baseline tests in 29/29 files,
- fresh 7/7 integration controls,
- 79,633/79,633 formula-level K16/K20 leaf replays,
- independent geometry/mutation audit PASS,
- campaign firewall/mutation audit PASS,
- 23/23 consolidated source hashes matching.

`portable_verify.py` was also executed from a separate clean directory containing only the copied campaign package and declared Python dependencies; it passed. The final archive hash is stored in the sidecar `FIN_R7N_HANDOFF_20260920.zip.sha256`. Full all-leaf formula replay scripts are included for deeper independent recomputation.

## 10. Ranked next atoms

1. **Physical-coupled Target-P enclosure:** preserve common `(r,s,t,y)` dependencies, particularly `t`, `t^3`, `t^4` and `t^(2+/-sqrt(3))`, to attack the remaining 36.31% compact residual.
2. **Optional amplitude-stability phase theorem:** only if robustness around the exact fixed fixture is scientifically needed; pay root, index and complement margins uniformly.
3. **Validated full-X7 exclusion primitive:** develop a genuinely seven-dimensional method before launching another stationary-complement campaign.

No next campaign is started automatically by this handoff.
