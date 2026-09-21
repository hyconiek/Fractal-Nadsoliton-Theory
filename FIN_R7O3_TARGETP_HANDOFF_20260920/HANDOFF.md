# FIN R7O3 Target-P completion handoff

## 0. Identity and provenance

- Date: 2026-09-20.
- Predecessor R7N ZIP required identity: `ba56cc0ecf69ac970e1d6fe719e0770753437213a4ef345949ffcf1411941d9c`.
- R7O2 continuation identity: `b3d04f6b25c451ddeac0adcf22adfa175b7e3950db1bff22bffa9e2a092e7731`.
- Frozen original residual: 5,432 parents, recorded in `inherited/frozen_r7n_residual_5432.json`.
- No source archive or authoritative AGENTS.md was modified.

## 1. Executive result

**Target P is a fully verified theorem candidate for supervisory review.**

Within the declared shared nonnegative-field C4 family, the proof set establishes `lambda2(M4)<=67/250`, conditional only on the already accepted R7N tail prerequisites outside the compact hull. R7O3 removes the complete 5,432-parent compact residual: 5,432/5,432 parents are certified closed and no terminal residual remains.

The strongest new statement is deliberately limited to Target P. Target S at sigma and full X7 remain open.

## 2. Primitive validation

The proof is based on the centered-moment Loewner inequality and min-max argument in `proofs/R7O3-006_centered_moment.md`. The transformed shared-parameter jets and `t^sqrt(3)` enclosure are documented in `proofs/R7O3-007_transformed_interval_jets.md`. `certify_fixed` accepts saved rational `B,c`, verifies exact rank/Gram evidence, and performs no eigensolver/optimizer call. Arithmetic comparison and negative controls are in `results/R7O3-010_arithmetic_validation.json` and the mutation suite.

## 3. Canonical parent registry

| State | Parent count | Whole-parent fraction of compact hull | Terminal residual fraction |
|---|---:|---:|---:|
| CERTIFIED_CLOSED | 5,432 | `23991168798917221399653867639993952899253335329138479599153946323/66072271890625000000000000000000000000000000000000000000000000000` ≈ 0.363104947241 | 0 |
| PARTIAL_REPAIR_REQUIRED | 0 | 0 | 0 |
| UNPROCESSED | 0 | 0 | 0 |
| INVALID | 0 | 0 | 0 |

All original IDs are accounted for exactly once. Replacement lineage is normalized into one active proof version per parent. The prior exploratory R7O2 2.26% number is not used in this canonical accounting.

## 4. First pass

The original R7O3 queue contained 1,972 genuinely unprocessed parent IDs, including the noncontiguous holes 3454 and 3458. The fixed first-pass policy used whole-cell certification followed by the bounded `t -> r -> s -> t -> r` tree. Work was published atomically per parent; interrupted batches did not mark unfinished parents safe.

All originals eventually received a produced state. The first pass left 81 partial parents / 255 unresolved terminals, which were frozen before repair.

## 5. Repairs

Repair-only work preserved each verified safe prefix and replaced only unresolved subtrees. After bounded local fallback, all 81 partial parents became closed. The active normalized proof set contains 913 R7O3 repair SAFE leaves and zero unresolved terminals. No surviving-cell analysis or counterexample certification was needed.

## 6. Independent replay

- Geometry: 5,432/5,432 parent covers PASS; no gaps or interior overlaps.
- Producer formula replay: 12,425/12,425 PASS.
- Clean-directory formula replay: 12,425/12,425 PASS with exact rational-bound match.
- Compact global join: 18,663 exact R7N leaves reconstruct the declared root hull.
- Parent identity: all 5,432 former residual cells match exactly.
- Accepted tail prerequisite hashes: PASS.
- Hostile mutations: 12/12 rejected.
- Final read-only verifier: `global_pass=true`, zero errors.

## 7. Global join and theorem statement

The 13,231 previously accepted compact SAFE leaves plus the now-certified 5,432 residual parents form the entire compact hull. Together with the accepted R7N tails this yields the Target-P theorem candidate `lambda2(M4)<=67/250` on the declared shared nonnegative-field family.

Therefore `index_negative(I4/g-M4)<=1` for supplied `0<g<=250/67`. At `g=250/67` a zero eigenvalue is not excluded by this non-strict bound.

No statement at sigma or in full X7 is exported.

## 8. Counterexamples / bounded failures

No admissible Target-P counterexample survives. There is no final resource-stop or mathematical residual. Historical failed certificates were method failures handled by later physical-coupled trees; they are not counterexamples.

## 9. Full R7O3-001--040 task ledger

| ID | Execution state | Scientific state | Artifacts | Remaining atom |
|---|---|---|---|---|
| R7O3-001 | DONE | VERIFIED | INPUTS.json, STATE_MAP.md | — |
| R7O3-002 | DONE | VERIFIED | INPUTS.json | — |
| R7O3-003 | DONE | VERIFIED | src/fixed_witness_checker.py, verify.py | — |
| R7O3-004 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-005 | DONE | VERIFIED | results/R7O3-036_final_accounting.json | — |
| R7O3-006 | DONE | VERIFIED | proofs/R7O3-006_centered_moment.md | — |
| R7O3-007 | DONE | VERIFIED | proofs/R7O3-007_transformed_interval_jets.md | — |
| R7O3-008 | DONE | VERIFIED | src/fixed_witness_checker.py | — |
| R7O3-009 | DONE | VERIFIED | certificates/active_leaf_certificates.jsonl | — |
| R7O3-010 | DONE | VERIFIED | results/R7O3-010_arithmetic_validation.json | — |
| R7O3-011 | DONE | VERIFIED | certificates/reconstructed_parent_trees.jsonl, src/verification_core.py | — |
| R7O3-012 | DONE | VERIFIED | parent_registry.json, certificates/reconstructed_parent_trees.jsonl | — |
| R7O3-013 | DONE | VERIFIED | results/R7O3-010_arithmetic_validation.json | — |
| R7O3-014 | DONE | VERIFIED | results/R7O3-032_full_formula_replay.json | — |
| R7O3-015 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-016 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-017 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-018 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-019 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-020 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-021 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-022 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-023 | DONE | VERIFIED | results/R7O3-032_full_formula_replay.json, certificates/reconstructed_parent_trees.jsonl | — |
| R7O3-024 | OUT_OF_SCOPE | NOT_APPLICABLE | — | No surviving leaves after R7O3 repairs; arithmetic validation was already paid by R7O3-010. |
| R7O3-025 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-026 | OUT_OF_SCOPE | NOT_APPLICABLE | — | No surviving terminals after R7O3-025. |
| R7O3-027 | OUT_OF_SCOPE | NOT_APPLICABLE | — | No surviving terminals after R7O3-025. |
| R7O3-028 | OUT_OF_SCOPE | NOT_APPLICABLE | — | No surviving terminals after R7O3-025. |
| R7O3-029 | OUT_OF_SCOPE | NOT_APPLICABLE | — | No admissible counterexample candidate survived; all parents certified closed. |
| R7O3-030 | DONE | VERIFIED | parent_registry.json | — |
| R7O3-031 | DONE | VERIFIED | parent_registry.json, certificates/reconstructed_parent_trees.jsonl | — |
| R7O3-032 | DONE | VERIFIED | results/R7O3-032_full_formula_replay.json, results/R7O3-037_clean_directory_math_replay.json | — |
| R7O3-033 | DONE | VERIFIED | results/R7O3-033_independent_global_join.json | — |
| R7O3-034 | DONE | VERIFIED | results/R7O3-034_hostile_mutations.json | — |
| R7O3-035 | DONE | VERIFIED | THEOREM_TARGET_P.md | — |
| R7O3-036 | DONE | VERIFIED | results/R7O3-036_final_accounting.json | — |
| R7O3-037 | DONE | VERIFIED | results/R7O3-037_clean_directory_math_replay.json, REPLAY.md | — |
| R7O3-038 | DONE | ARTIFACT_COMPLETE | REPORT.md, AGENTS_PROPOSED_PATCH.md | — |
| R7O3-039 | DONE | ARTIFACT_COMPLETE | MANIFEST.sha256 | — |
| R7O3-040 | DONE | ARTIFACT_COMPLETE | HANDOFF.md | — |

## 10. Portability

The package includes all active parent trees, 12,425 fixed-witness certificates, frozen proof inputs, clean checker code, complete replay shards, the exact R7N compact partition used for the join, accepted tail dependencies and source/provenance hashes. Executable proof code contains no required producer-specific `/mnt/data` path. `MANIFEST.sha256` covers package contents except itself.

## 11. Integration proposal only

`AGENTS_PROPOSED_PATCH.md` and `ACCEPTED_CLAIM_CANDIDATE.md` are proposals. They have **not** been applied to the authoritative repository. Supervisory review should decide whether to promote the theorem candidate.

## 12. Next atoms

1. Independent supervisory review/integration of R7O3 Target P using the portable replay package.
2. Only after acceptance, decide whether a new campaign should attack the sharper Target S (`sigma`) ceiling; Target P must not be silently relabeled as Target S.
3. Full-X7/global-energy questions remain separate and should start only with an explicit new plan.

No background worker is intentionally left running by this handoff.
