# Strict alpha chi_11 d1+d6 local necessity certificate

Status: `d1-d6-is-locally-necessary-and-upward-isolated-inside-chi11-shell-exclusion-grammar`

## Finite model

- Ring: `Z_12`
- Active count: `5`
- chi_11 kernel units: `[1, 11]`
- Center forbidden shells: `[1, 6]`

## Local boundary summary

- Center selects A5/d5 uniquely: `True`
- Deletion neighbors: `2`
- Addition neighbors: `4`
- All deletions overselect: `True`
- All additions UNSAT: `True`
- delete d1 solution count: `192`
- delete d6 solution count: `36`
- addition solution counts: `{'add_d2': 0, 'add_d3': 0, 'add_d4': 0, 'add_d5': 0}`

## Proof certificate

- `finite_domain`: Each local system enumerates all C(12,5)=792 supports exactly under the chi_11-stable shell-exclusion grammar.
- `center_certificate`: The center {d1,d6} has 12 solutions, one dihedral orbit, and only histogram [0,3,2,1,4,0], so it selects A5/d5.
- `d1_necessity`: Deleting d1 leaves only d6 and gives 192 solutions over 9 histogram classes, including A1/contiguous; hence d1 is necessary for specificity.
- `d6_necessity`: Deleting d6 leaves only d1 and gives 36 solutions over 3 gap necklaces; hence d6 is necessary for exact cover uniqueness.
- `maximality_boundary`: Adding any one missing shell to {d1,d6} gives UNSAT counts [0, 0, 0, 0], so the selector is upward-isolated.
- `conditional_scope`: Necessity is local and conditional on the chi_11/shell-label premise plus the pure shell-exclusion exact-cover grammar; no strict source theorem is claimed.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.
- The result is a conditional finite local necessity certificate inside an axiom-augmented grammar, not a strict selector-origin theorem.
- No QW-2191 discharge.
- No ToE closure.
