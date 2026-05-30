# Strict alpha full-Aut forbidden-shell lattice audit

Status: `full-aut-forbidden-shell-lattice-classified-previous-d1-d5-d6-closure-is-minimal-but-not-unique-unsat-blocker`

## Finite model

- Ring: `Z_12`
- Aut units: `[1, 5, 7, 11]`
- Shell orbits: `[[1, 5], [2], [3], [4], [6]]`
- Target active count: `5`
- Previous full-Aut closure: `[1, 5, 6]`

## Lattice summary

- Total full-Aut masks: `32`
- k=5 UNSAT masks: `25`
- k=5 SAT masks: `7`
- Minimal k=5 UNSAT masks: `5`
- Minimal blockers: `[[4], [2, 3], [3, 6], [1, 2, 5], [1, 5, 6]]`
- Previous closure is minimal blocker: `True`
- Previous closure alpha: `3`
- Previous closure k=5 count: `0`

## Proof certificate

- `full_aut_lattice`: Aut(Z_12) has shell orbits [[1,5],[2],[3],[4],[6]], so exactly 32 full-Aut-invariant forbidden-shell masks are audited.
- `classification_counts`: Of 32 masks, 25 make k=5 UNSAT and 7 leave at least one 5-support.
- `minimal_blockers`: There are 5 inclusion-minimal full-Aut k=5 blockers: [[4], [2, 3], [3, 6], [1, 2, 5], [1, 5, 6]].
- `previous_closure_position`: The previous closure [1, 5, 6] is one minimal blocker with alpha=3 and k=5 count=0.
- `non_uniqueness_warning`: Because other minimal full-Aut blockers exist, the full-Aut UNSAT fact alone does not derive the d1+d6 clause origin or chi_11 source.
- `conditional_selector_reading`: The successful A5/d5 exact-cover selector still requires the non-full-Aut d1-vs-d5 shell-label premise before the full-Aut closure is avoided.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, or exact-cover clauses from strict nadsoliton geometry.
- The result is a finite full-Aut forbidden-shell lattice audit, not a selector-origin theorem.
- No QW-2191 discharge.
- No ToE closure.
