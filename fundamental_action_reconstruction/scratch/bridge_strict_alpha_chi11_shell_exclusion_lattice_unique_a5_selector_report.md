# Strict alpha chi_11 shell-exclusion lattice unique A5 selector audit

Status: `chi11-shell-exclusion-lattice-has-unique-a5-d5-selector-d1-d6-conditional-on-shell-label-premise`

## Finite model

- Ring: `Z_12`
- Active count: `5`
- chi_11 kernel units: `[1, 11]`
- Folded shells: `[1, 2, 3, 4, 5, 6]`
- Target A5/d5 histogram: `[0, 3, 2, 1, 4, 0]`

## Lattice summary

- Total chi_11 shell masks: `64`
- A5/d5 selecting mask count: `1`
- A5/d5 selecting forbidden shells: `[[1, 6]]`
- A1/contiguous selecting mask count: `1`
- A1/contiguous selecting forbidden shells: `[[5, 6]]`
- UNSAT mask count: `51`
- SAT mask count: `13`
- Minimal UNSAT forbidden shells: `[[4], [1, 2], [2, 3], [2, 5], [3, 6], [1, 5, 6]]`

## Proof certificate

- `finite_domain`: For every one of the 64 chi_11-stable shell-exclusion masks, all C(12,5)=792 supports are enumerated exactly.
- `chi11_shell_stability`: Units {1,11} fix every folded shell label d1..d6; unlike full Aut, they do not force d1 and d5 to share a clause orbit.
- `unique_a5_selector`: Exactly 1 mask selects the A5/d5 orbit; its forbidden shells are [1, 6].
- `conjugate_boundary`: Exactly 1 mask selects the A1/contiguous orbit; this records the d5+d6 conjugate boundary separately from the A5 selector.
- `unsat_boundary`: There are 51 masks with no 5-support at all, so uniqueness of A5/d5 is not just generic overconstraint.
- `conditional_scope`: This is uniqueness only after adjoining the chi_11/shell-label premise and restricting the grammar to pure shell-exclusion exact-cover clauses.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, or exact-cover clauses from strict nadsoliton geometry.
- The result is a conditional finite uniqueness certificate inside an axiom-augmented grammar, not a strict selector-origin theorem.
- No QW-2191 discharge.
- No ToE closure.
