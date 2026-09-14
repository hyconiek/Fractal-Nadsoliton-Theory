# FIN rank-seven follow-up research handoff

## 0. Identity and execution state
- Campaign version and date: post-handoff 128-task campaign, final assembly 2026-09-14.
- Repository baseline: supplied portable handoff bundle; no Git repository was available in this execution environment.
- Environment: see `environment.json`; final campaign replay recorded on Python 3.13.5.
- Executed task range: R7P-001 through R7P-128. Terminal task count: **128/128**.
- Package manifest: `FINAL_MANIFEST.sha256` (generated after clean-copy smoke replay).
- Running/stopped process inventory: no campaign worker is intended to remain running after archive creation.

## 1. Executive result
- **Most important genuinely new theorem:** the campaign rebuilt the exact boundary-Ising ceiling and extended rigorous four-amplitude control through the shared-field intraparity theorem, an off-face cone, and a large-J5 tail; separately it established a rigorous global energetic bracket `g_global in [2.8934,3.71835]`.
- **Most important certified counterexample:** a stationary full-seven-coordinate root at g=5 has H7 index at least 2, refuting the stationary-only universal index-one conjecture.
- **Most important reproduced numerical finding:** the phase program recovers 60 quartic critical candidates and locally isolates all 60 plus 60 nearby full roots, while preserving the unresolved complement rather than claiming exhaustion.
- **Most important unresolved gap:** the compact residual core of the positive-orthant four-amplitude ceiling; global phase/stationary exhaustion and unique global minimizer orbit also remain open.
- **Already known before this campaign:** the strict spectral provider, weighted Hodge identities at intake scope, Gibbs duality, local crossing candidate, corrected face formulas, separate quantum package results, and the nonstationary full-7D curvature counterexample.

## 2. Source and coordinate contract
- Input hashes: `INPUT_HASHES.json`; one absent historical ZIP is explicitly declared in `NONREPLAYED_INPUTS.json`.
- Baseline audit: `references/starting_audit/`.
- `X7` is the full seven-column real Fourier feature matrix; `C4` is only its reflection-fixed cosine/alternating restriction.
- `p` is always a full 12-label probability distribution. Gibbs completion does not impose `p-u in Range(A7)`.
- Exact strict spectral intervals are distinct from rounded display values and from historical rational discord certificates.

## 3. Complete task ledger

| Task ID | Execution status | Claim status | Main output | Remaining atom |
|---|---|---|---|---|
| R7P-001 | DONE | NOT_A_SCIENTIFIC_CLAIM | `environment.json` | — |
| R7P-002 | DONE | NOT_A_SCIENTIFIC_CLAIM | `logs/verify_initial.log` | 56 inherited regressions remain historical PASS only because their full source packages were not supplied locally; this does not block rank7 exact-layer dependencies. |
| R7P-003 | DONE | NOT_A_SCIENTIFIC_CLAIM | `CLAIMS.json` | — |
| R7P-004 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-004_missing_evidence_matrix.json` | — |
| R7P-005 | DONE | NOT_A_SCIENTIFIC_CLAIM | `src/schema.py` | — |
| R7P-006 | DONE | NOT_A_SCIENTIFIC_CLAIM | `verify.py` | — |
| R7P-007 | DONE | NOT_A_SCIENTIFIC_CLAIM | `src/resumable_cover.py` | — |
| R7P-008 | DONE | NOT_A_SCIENTIFIC_CLAIM | `STATE_MAP.md` | — |
| R7P-009 | DONE | NUMERICAL_REPRODUCED | `src/model.py` | — |
| R7P-010 | DONE | NUMERICAL_REPRODUCED | `src/model.py` | — |
| R7P-011 | DONE | EXACT_PROVED | `proofs/R7P-011_primal_dual_theorem.md` | — |
| R7P-012 | DONE | NOT_A_SCIENTIFIC_CLAIM | `src/derivatives.py` | — |
| R7P-013 | DONE | NOT_A_SCIENTIFIC_CLAIM | `src/intervals.py` | — |
| R7P-014 | DONE | EXACT_PROVED | `proofs/R7P-014_stationary_inertia_transfer.md` | — |
| R7P-015 | DONE | EXACT_PROVED | `proofs/R7P-015_compact_domains.md` | — |
| R7P-016 | DONE | NUMERICAL_NEW | `proofs/R7P-016_domain_charts.md` | — |
| R7P-017 | DONE | INTERVAL_CERTIFIED | `src/face_certificate.py` | — |
| R7P-018 | DONE | EXACT_PROVED | `src/face_certificate.py` | — |
| R7P-019 | DONE | INTERVAL_CERTIFIED | `certificates/R7P-017_021_face_certificate.json` | — |
| R7P-020 | DONE | INTERVAL_CERTIFIED | `certificates/R7P-017_021_face_certificate.json` | — |
| R7P-021 | DONE | INTERVAL_CERTIFIED | `certificates/R7P-017_021_face_certificate.json` | — |
| R7P-022 | DONE | INTERVAL_CERTIFIED | `src/face_gap.py` | — |
| R7P-023 | DONE | INTERVAL_CERTIFIED | `src/parity_asymptotics.py` | — |
| R7P-024 | DONE | INTERVAL_CERTIFIED | `src/face_api.py` | — |
| R7P-025 | DONE | NUMERICAL_REPRODUCED | `—` | — |
| R7P-026 | DONE | INTERVAL_CERTIFIED | `—` | — |
| R7P-027 | DONE | INTERVAL_CERTIFIED | `—` | — |
| R7P-028 | DONE | INTERVAL_CERTIFIED | `—` | — |
| R7P-029 | DONE | INTERVAL_CERTIFIED | `—` | — |
| R7P-030 | DONE | NUMERICAL_REPRODUCED | `—` | — |
| R7P-031 | DONE | INTERVAL_CERTIFIED | `—` | — |
| R7P-032 | DONE | NUMERICAL_NEW | `results/R7P-032_local_branch_diagram.json` | Overlapping validated continuation tubes joining the isolated fold and equal-energy neighborhoods. |
| R7P-033 | DONE | EXACT_PROVED | `proofs/R7P-033_full7_everywhere_counterexample.md` | — |
| R7P-034 | DONE | EXACT_PROVED | `src/two_harmonic.py` | — |
| R7P-035 | DONE | NUMERICAL_NEW | `results/R7P-034_035_two_harmonic_stationary_atlas.json` | — |
| R7P-036 | DONE | COUNTEREXAMPLE_CERTIFIED | `src/two_harmonic_certificate.py` | — |
| R7P-037 | DONE | NUMERICAL_REPRODUCED | `results/R7P-037_full7_stationary_atlas.json` | numerical saturation only; complement not covered |
| R7P-038 | DONE | INTERVAL_CERTIFIED | `—` | — |
| R7P-039 | DONE | INTERVAL_CERTIFIED | `results/R7P-039_restricted_full_mismatch_table.json` | Complete mismatch table for remaining certified/numerical reflection-fixed roots as they are added to atlas. |
| R7P-040 | DONE | COUNTEREXAMPLE_CERTIFIED | `—` | — |
| R7P-041 | DONE | EXACT_PROVED | `results/R7P-041_045_boundary_ising.json` | — |
| R7P-042 | DONE | EXACT_PROVED | `results/R7P-041_045_boundary_ising.json` | — |
| R7P-043 | DONE | EXACT_PROVED | `results/R7P-041_045_boundary_ising.json` | — |
| R7P-044 | DONE | INTERVAL_CERTIFIED | `results/R7P-041_045_boundary_ising.json` | — |
| R7P-045 | DONE | INTERVAL_CERTIFIED | `results/R7P-041_045_boundary_ising.json` | — |
| R7P-046 | DONE | EXACT_PROVED | `results/R7P-046_050_boundary_cover_foundation.json` | — |
| R7P-047 | DONE | COUNTEREXAMPLE_CERTIFIED | `results/R7P-046_050_boundary_cover_foundation.json` | — |
| R7P-048 | DONE | NOT_A_SCIENTIFIC_CLAIM | `certificates/R7P-048_boundary_proof_spec.json` | — |
| R7P-049 | DONE | EXACT_PROVED | `results/R7P-046_050_boundary_cover_foundation.json` | — |
| R7P-050 | DONE | EXACT_PROVED | `src/boundary_cover.py` | — |
| R7P-051 | DONE | INTERVAL_CERTIFIED | `results/R7P-051_boundary_cover_pilot.json` | — |
| R7P-052 | DONE | INTERVAL_CERTIFIED | `certificates/R7P-052_local_equality.json` | — |
| R7P-053 | DONE | INTERVAL_CERTIFIED | `results/R7P-053_refined_boundary_cover.json` | — |
| R7P-054 | DONE | INTERVAL_CERTIFIED | `results/R7P-054_checker_result.json` | — |
| R7P-055 | DONE | INTERVAL_CERTIFIED | `certificates/R7P-055_boundary_theorem.json` | — |
| R7P-056 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-056_dependency_impact.json` | — |
| R7P-057 | DONE | EXACT_PROVED | `results/R7P-057_intraparity_shared_fields.json` | — |
| R7P-058 | DONE | INTERVAL_CERTIFIED | `results/R7P-057_059_intraparity.json` | — |
| R7P-059 | DONE | INTERVAL_CERTIFIED | `results/R7P-057_059_intraparity.json` | — |
| R7P-060 | DONE | NUMERICAL_REPRODUCED | `results/R7P-060_063_intraparity_closure.json` | — |
| R7P-061 | DONE | INTERVAL_CERTIFIED | `results/R7P-060_063_intraparity_closure.json` | — |
| R7P-062 | DONE | INTERVAL_CERTIFIED | `results/R7P-060_063_intraparity_closure.json` | — |
| R7P-063 | DONE | INTERVAL_CERTIFIED | `results/R7P-060_063_intraparity_closure.json` | — |
| R7P-064 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-064_relaxation_regressions.json` | — |
| R7P-065 | DONE | EXACT_PROVED | `certificates/R7P-065_off_face_target.json` | — |
| R7P-066 | DONE | NUMERICAL_REPRODUCED | `results/R7P-065_067_off_face.json` | — |
| R7P-067 | DONE | EXACT_PROVED | `certificates/R7P-067_compactification.json` | — |
| R7P-068 | DONE | INTERVAL_CERTIFIED | `certificates/R7P-068_first_order_partial.json` | — |
| R7P-069 | DONE | INTERVAL_CERTIFIED | `results/R7P-069_072_off_face_partial.json` | global compact residual core remains unresolved |
| R7P-070 | DONE | UNRESOLVED | `results/R7P-065_067_off_face.json` | no certified admissible violation found; search cannot prove absence |
| R7P-071 | DONE | CONDITIONAL_LEMMA | `proofs/R7P-069_072_off_face_final.md` | only union of certified 4D subdomains is proved |
| R7P-072 | DONE | NOT_A_SCIENTIFIC_CLAIM | `proofs/R7P-069_072_off_face_final.md` | global positive-orthant ceiling and full-7D transfer remain prohibited |
| R7P-073 | DONE | EXACT_PROVED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-074 | DONE | NUMERICAL_REPRODUCED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-075 | DONE | EXACT_PROVED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-076 | DONE | EXACT_PROVED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-077 | DONE | INTERVAL_CERTIFIED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-078 | DONE | EXACT_PROVED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-079 | DONE | INTERVAL_CERTIFIED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-080 | DONE | EXACT_PROVED | `proofs/R7P-073_080_cooperativity.md` | — |
| R7P-081 | DONE | NUMERICAL_REPRODUCED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-082 | DONE | EXACT_PROVED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-083 | DONE | EXACT_PROVED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-084 | DONE | EXACT_PROVED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-085 | DONE | NUMERICAL_REPRODUCED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-086 | DONE | UNRESOLVED | `proofs/R7P-081_090_phase_foundation.md` | Sharp uniform C2 remainder bound preserving phase/resonance cancellations. |
| R7P-087 | DONE | EXACT_PROVED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-088 | DONE | NUMERICAL_NEW | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-089 | DONE | NUMERICAL_REPRODUCED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-090 | DONE | INTERVAL_CERTIFIED | `proofs/R7P-081_090_phase_foundation.md` | — |
| R7P-091 | DONE | UNRESOLVED | `proofs/R7P-091_095_phase_continuation.md` | Certified exclusion of the 1272 retained complement leaves. |
| R7P-092 | DONE | INTERVAL_CERTIFIED | `proofs/R7P-091_095_phase_continuation.md` | — |
| R7P-093 | DONE | UNRESOLVED | `proofs/R7P-091_095_phase_continuation.md` | Direct full-function complement cover or quartic gradient gap exceeding a rigorous uniform gradient remainder. |
| R7P-094 | DONE | CONDITIONAL_LEMMA | `proofs/R7P-091_095_phase_continuation.md` | — |
| R7P-095 | DONE | INTERVAL_CERTIFIED | `proofs/R7P-091_095_phase_continuation.md` | — |
| R7P-096 | DONE | NUMERICAL_NEW | `results/R7P-096_angular_diagram.json` | — |
| R7P-097 | DONE | EXACT_PROVED | `results/R7P-097_104_global_frontier.json` | attainment/uniqueness of ratio infimum unresolved |
| R7P-098 | DONE | CONDITIONAL_LEMMA | `proofs/R7P-097_104_global_frontier.md` | lower range uses accepted upstream ST448 plus exact PSD order; no new full-domain rank7 cover |
| R7P-099 | DONE | INTERVAL_CERTIFIED | `results/R7P-097_104_global_frontier.json` | bracket remains wide and does not identify attaining orbit |
| R7P-100 | DONE | NUMERICAL_REPRODUCED | `results/R7P-037_full7_stationary_atlas.json` | 7D complement not exhausted |
| R7P-101 | DONE | UNRESOLVED | `results/R7P-097_104_global_frontier.json` | g=4 unique global D12 minimizer orbit not certified |
| R7P-102 | DONE | EXACT_PROVED | `results/R7P-097_104_global_frontier.json` | no theorem forces generic global minimizer into C4/two-harmonic/phase subspaces |
| R7P-103 | DONE | NUMERICAL_NEW | `results/R7P-097_104_global_frontier.json` | no validated trajectory tube; declared Euclidean gradient flow only |
| R7P-104 | DONE | NOT_A_SCIENTIFIC_CLAIM | `proofs/R7P-097_104_global_frontier.md` | global uniqueness and exhaustion remain open |
| R7P-105 | DONE | EXACT_PROVED | `results/R7P-105_112_passive_finiteN.json` | none within weighted Hodge assumptions |
| R7P-106 | DONE | NUMERICAL_REPRODUCED | `results/R7P-105_112_passive_finiteN.json` | 28 distinct positive eigenvalues not upgraded to exact multiplicity theorem |
| R7P-107 | DONE | EXACT_PROVED | `results/R7P-105_112_passive_finiteN.json` | none for the declared independent-walker model |
| R7P-108 | DONE | EXACT_PROVED | `results/R7P-105_112_passive_finiteN.json` | no universal finite-N Gaussian theorem; OU is limiting/equilibrium linear-noise |
| R7P-109 | DONE | CONDITIONAL_LEMMA | `results/R7P-105_112_passive_finiteN.json` | visible McMillan degree five under retained symmetric residue structure |
| R7P-110 | DONE | EXACT_PROVED | `proofs/R7P-105_112_passive_finiteN.md` | only symmetric passive class; driven/nonnormal models separate |
| R7P-111 | DONE | CONDITIONAL_LEMMA | `results/R7P-105_112_passive_finiteN.json` | beta, J, g, mediator rank, pump and clock remain supplied |
| R7P-112 | DONE | NOT_A_SCIENTIFIC_CLAIM | `proofs/R7P-105_112_passive_finiteN.md` | gain/source problem remains open |
| R7P-113 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-113_120_quantum_bridge.json` | scope matrix frozen |
| R7P-114 | DONE | EXISTING_ACCEPTED_SCOPED | `results/R7P-113_120_quantum_bridge.json` | 5e-16 replay provenance gap for historical delta_L |
| R7P-115 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-113_120_quantum_bridge.json` | no supplied causal coupling map |
| R7P-116 | DONE | EXISTING_ACCEPTED_SCOPED | `proofs/R7P-113_120_quantum_bridge.md` | existing theorem replay only |
| R7P-117 | DONE | CONDITIONAL_LEMMA | `results/R7P-113_120_quantum_bridge.json` | loading robustness conditional on accepted historical gap certificates |
| R7P-118 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-113_120_quantum_bridge.json` | resource table complete for discussed examples |
| R7P-119 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-113_120_quantum_bridge.json` | no new legacy comparison needed |
| R7P-120 | DONE | NOT_A_SCIENTIFIC_CLAIM | `proofs/R7P-113_120_quantum_bridge.md` | no causal spectral bridge; one narrow loading robustness result only |
| R7P-121 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-121_claim_traceability.json` | — |
| R7P-122 | DONE | NOT_A_SCIENTIFIC_CLAIM | `results/R7P-122_mutation_review.json` | — |
| R7P-123 | DONE | NOT_A_SCIENTIFIC_CLAIM | `verification_final.json` | — |
| R7P-124 | DONE | NOT_A_SCIENTIFIC_CLAIM | `RESULT_TABLES.md` | — |
| R7P-125 | DONE | NOT_A_SCIENTIFIC_CLAIM | `REPORT_FINAL.md` | — |
| R7P-126 | DONE | NOT_A_SCIENTIFIC_CLAIM | `AGENTS_PATCH.md` | — |
| R7P-127 | DONE | NOT_A_SCIENTIFIC_CLAIM | `REPLAY.md` | — |
| R7P-128 | DONE | NOT_A_SCIENTIFIC_CLAIM | `HANDOFF.md` | — |

## 4. Accepted new theorems

### 4.1 Boundary-Ising global ceiling
The exact attainable boundary domain, its strata, the shifted-characteristic/inertia criterion, complete cover, local double-root treatment and independent checker establish the stated second-eigenvalue ceiling. Equality is confined to the certified compactified double-root case. See `proofs/R7P-046_056_boundary_global_theorem.md` and the G certificates/checker.

### 4.2 Shared-field intraparity ceiling
For shared J3,J4,J5,J6>=0, `lambda2(W_par)<=sigma_*`. The proof pays the conditional-sector coupling constraints and dominant-mass/envelope bounds; it does not optimize parity classes independently. See `proofs/R7P-057_064_intraparity_closure.md`.

### 4.3 Partial positive-orthant four-amplitude theorem
The exact extreme face, the certified tangent cone of radius 1/8192, and the full tail `exp(-J5/2)<=2^-11` satisfy the ceiling. The residual compact core is not covered. See `proofs/R7P-069_072_off_face_final.md`.

### 4.4 Stationary full-7D counterexample theorem
At g=5 an interval-isolated stationary root has at least two negative full-H7 directions. This refutes the stationary-only universal index-one conjecture but says nothing about all stationary roots at all gains.

### 4.5 Global energetic bracket
With `g_global=inf 2D/Q`, PSD ordering transfers the certified full-model uniform region only in the permitted lower-bound direction, while a rational rank-seven competitor is strictly negative at exact g=3.71835. Thus `g_global in [2.8934,3.71835]`. The attaining orbit and uniqueness remain open.

### 4.6 Other scoped theorems
The campaign also certifies the pure-k6 first angular instability, local phase root boxes/indexes, passive finite-N generator identities, passive Schur/positive-real memory in the symmetric passive class, physical cooperativity on its paid finite-group hypotheses, and the narrow gamma-loading quantum robustness statement recorded in O. Each retains its own domain/nontransfer restrictions.

## 5. Counterexamples and rejected conjectures
- The everywhere full-seven-coordinate Hessian index<=1 claim is refuted by the intake nonstationary witness.
- The stationary-only universal index<=1 claim is separately refuted at g=5.
- Relaxing any one of the key physical boundary-Ising inequalities produces negative-control violations; these do not refute the physical domain.
- A simple norm/Weyl off-face bound is grossly too loose and is retained as a method counterexample.
- Fixed Perron-compression directions do not establish the global 4D ceiling.
- Numerical saturation of phase/stationary searches is never used as an exhaustion theorem.
- The imported quartic-fold decimal is not reproduced by the reconstructed generator; the discrepancy is preserved.

## 6. Numerical findings not yet certified
- Full-7D stationary atlas: 3/3/3/15 D12 orbit candidates at g=3.7/3.7183449/4/5 under two recorded batches; no complement cover.
- g=4 localized stationary orbit is the lowest found candidate, not a globally unique minimizer theorem.
- Declared Euclidean gradient-flow branches from the g=4 saddle numerically approach uniform/localized endpoints; no validated trajectory tube.
- Distinct positive cycle spectrum count 28 is numerical only.
- Sampled C2 phase remainder diagnostics are not substituted for uniform bounds.

## 7. Root and continuation atlas
Certified local equal-energy/fold/root boxes are in `certificates/` and corresponding proofs. R7P-032 retains certified nodes plus numerical continuation segments without claiming validated joins. R7P-037 is discovery saturation only. D12 orbit/stabilizer data are saved with the numerical atlas.

## 8. Cover certificates and unresolved domains
- Boundary-Ising: global cover complete; independent checker and mutations pass.
- Four-amplitude: extreme face + local cone rho=1/8192 + t<=2^-11 tail certified. Residual compact core explicitly unresolved.
- Quartic phase torus: 60 root neighborhoods locally certified; bounded complement pass leaves 1272 unresolved leaves.
- Full phase torus: 60 local corresponding full roots; extra roots elsewhere are not excluded.

## 9. Phase results
The fixed-amplitude fixture was reconstructed from the angular-radius convention. Cubic locks, quartic cumulants, full derivatives and root candidates were independently rebuilt. All 60 quartic roots and 60 nearby full roots have local boxes; global count remains open. The pure-k6 constrained-sphere instability is certified. Full angular fold/coexistence/barrier are reconstructed with proof levels recorded in R7P-096.

## 10. Full-seven-coordinate/global results
The nonstationary index-two witness is preserved; a stationary index-two witness is newly certified. Full stationary exhaustion is open. The rigorous energetic bracket is `[2.8934,3.71835]`; first-transition orbit and g=4 unique global minimizer remain unresolved. The chosen gradient flow is an analysis convention, not a physical law.

## 11. Conditional interpretation and resource assumptions
Passive walkers, Hodge cycle space and passive memory do not supply active gain. The conditional Gibbs model supplies a statistical realization only after g=beta J is given. Quantum preparation/access assumptions remain separate from the classical rank-seven landscape. No selector, pump, clock, apparatus or causal quantum-to-localization map is produced by this campaign.

## 12. Failures, resource stops, and do-not-repeat list
- Do not rerun the old boundary proof on a relaxed domain and treat relaxation violations as physical counterexamples.
- Do not use the simple Weyl norm bound to claim the off-face ceiling.
- Do not infer phase or stationary exhaustion from repeated random/Sobol saturation.
- Do not convert the 1272 unresolved quartic-cover leaves into “numerical noise.”
- Do not treat the current rounded strict provider as independently replaying the historical `delta_L`; its lower endpoint is 5e-16 lower.
- Do not include stale `__pycache__` in portable archives; it was observed to retain old absolute source paths.
- Do not regenerate the absent historical source ZIP merely to silence the input ledger.

## 13. Verification and portability
- Canonical verifier: `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.:..:src python verify.py`.
- Claim traceability: `python src/final_audit.py`.
- Campaign suite: **94 passed, 0 failed, 0 skipped**, run as fresh-process shards; see `verification_final.json` and `replay_sharded.py`.
- Mutation review: `results/R7P-122_mutation_review.json`.
- Clean-copy smoke replay: manifest PASS, verifier/traceability PASS, and 27 representative tests PASS; see `results/R7P-127_clean_replay.json`.
- Historical baseline: 19 intake tests + 56 inherited regressions under the recorded Python 3.12.3 environment.
- Known nonreplayed material: historical source ZIP listed in `NONREPLAYED_INPUTS.json`; extracted constituent evidence is included.

## 14. Proposed AGENTS.md update
See `AGENTS_PATCH.md`. It is a reviewable patch; it is not silently substituted for the earlier intake/audit history.

## 15. Ranked next frontier
1. **Residual 4D core:** produce either a complete dependency-aware inertia cover of the compact residual core or a certified admissible counterexample. Completion = zero unresolved cells or one certified physical witness.
2. **Phase complement exhaustion:** close the 1272 quartic unresolved leaves and then transfer/redo the complement for the full function. Completion = exact zero remaining complement cells.
3. **Full-7D stationary exhaustion:** replace numerical atlas saturation by a complete symmetry-aware complement exclusion at one fixed gain (g=4 is the clearest target). Completion = certified all-root orbit list or explicit bounded unresolved cells.
4. **Narrow global energetic bracket:** develop a rigorous full-domain lower bound above 2.8934 and/or certified competitor below 3.71835 without assuming a symmetry-reduced minimizer. Completion = strictly narrower interval.
5. **Physical source law:** only if new external scientific premises are supplied, state and test an explicitly typed active law that determines gain/clock/selector. Completion = sourced law + mathematical coupling map; otherwise this remains outside the current model.
