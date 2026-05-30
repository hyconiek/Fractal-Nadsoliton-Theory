# D12 chi_11 character module certificate

Status: `exact-d12-quotient-chi11-module-rank-13-but-requires-unit-axis-premise`

## Finite model

- Ring: `Z_12`
- Enumerated supports: `792`
- D12 orbit count: `38`
- Residual unit-5 permutation size: `38`

## Module summary

- Tau two-cycles: `13`
- Tau fixed D12 orbits: `12`
- Integer chi_11 module rank: `13`
- Full-Aut invariant intersection rank: `0`
- Global nonzero ±1 character count: `0`
- Ternary chi_11-covariant count: `1594323`
- Branch-normalized ternary count: `531441`

## Branch generator

- Cycle index: `0`
- Orbit pair: `[0, 37]`
- Basis values low/high: `[-1, 1]`
- Low/high gap necklaces: `[1, 1, 1, 1, 8]` / `[2, 2, 3, 2, 3]`
- Requires premise: unit-axis orientation selecting the D12 quotient component pair; not full-Aut strict-core provenance

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are quotient-enumerated by D12.
- `residual_action`: Modulo D12, units 1 and 11 are trivial while units 5 and 7 induce the same involution on D12 orbits.
- `module_equations`: A D12-invariant chi_11-covariant support function is exactly a function f on D12 orbits satisfying f(tau(O))=-f(O), tau=unit5.
- `rank_certificate`: The tau action has 13 two-cycles and 12 fixed points; hence the integer chi_11 module has rank 13 and fixed-point values are forced to 0.
- `full_aut_intersection`: Full-Aut invariance requires f(tau(O))=f(O); together with chi_11 covariance this gives f(O)=0 on every orbit over characteristic zero.
- `branch_generator`: The A1/A11 orbit and A5/A7 orbit form one tau two-cycle, so there is a branch-local generator with values -1 and +1 on that pair.
- `not_strict_source`: Selecting that generator requires the D12 quotient and unit-axis orientation; it is not exported by full-Aut strict support data alone.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
- Result is a complete D12-quotient character-module certificate, not a complete strict-source provenance theorem.
