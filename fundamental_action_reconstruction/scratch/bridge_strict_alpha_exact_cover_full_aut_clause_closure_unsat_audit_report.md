# Strict alpha exact-cover full-Aut clause-closure UNSAT audit

Status: `full-aut-closure-of-d1-d6-exact-cover-is-unsat-chi_11-kernel-certificate-only`

## Finite model

- Ring: `Z_12`
- Active Boolean variables: `5` of 12
- Cardinality-compatible assignments checked: `792`
- Aut units: `[1, 5, 7, 11]`
- Base forbidden shells: `[1, 6]`

## Clause closure computation

- Orbit of `d1` under full Aut: `[1, 5]`
- Orbit of `d6` under full Aut: `[6]`
- Full-Aut closed forbidden shells: `[1, 5, 6]`
- chi_11-kernel closed forbidden shells: `[1, 6]`

## Exact enumeration results

- d1+d6 / chi_11-kernel solution count: `12`
- d1+d6 / chi_11-kernel dihedral orbit count: `1`
- d1+d6 selects A5/d5 histogram: `True`
- Full-Aut closure d1+d5+d6 solution count: `0`
- Full-Aut closure is UNSAT: `True`
- Aut-conjugate d5+d6 solution count: `12`
- Aut-conjugate d5+d6 histograms: `[[4, 3, 2, 1, 0, 0]]`

## Proof certificate

- `finite_domain`: All C(12,5)=792 Boolean cardinality-compatible supports are enumerated exactly.
- `base_certificate`: Forbidding d1 and d6 gives 12 solutions, one dihedral orbit, with histogram [0,3,2,1,4,0] (A5/d5).
- `unit_closure_step`: Under Aut(Z_12)={1,5,7,11}, the orbit of d1 is {d1,d5}; d6 is fixed. Full-Aut closure of d1+d6 is therefore d1+d5+d6.
- `unsat_result`: Forbidding d1,d5,d6 has 0 solutions among all 792 five-supports, so full-Aut closure destroys the exact-cover certificate.
- `conditional_positive`: The d1+d6 certificate survives exactly at the chi_11-kernel level {1,11}; it is not a full-Aut-invariant strict-core selector theorem.
- `source_obstruction`: A strict derivation would still need an internal source for the chi_11/shell-label reduction before using the successful exact-cover selector.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is used or claimed.
- No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.
- No theorem derives the chi_11-kernel, shell-label d1 vs d5, or unit-axis bit from strict nadsoliton geometry.
- The result is a finite full-Aut clause-closure UNSAT audit, not a selector-origin theorem.
- No QW-2191 discharge.
- No ToE closure.
