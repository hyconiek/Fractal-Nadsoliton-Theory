# Cyclic quotient character projection chi_11 certificate

Status: `translation-quotient-character-projection-rank-13-for-chi11-but-not-strict-source`

## Finite model

- Ring: `Z_12`
- Enumerated supports: `792`
- Translation orbit count: `66`
- Residual unit orbit count: `25`
- Residual unit orbit size counts: `{1: 4, 2: 11, 4: 10}`

## Trace / character summary

- Fixed translation-orbit counts by unit: `{'1': 66, '5': 10, '7': 14, '11': 10}`
- Trivial/full-Aut invariant dimension: `25`
- chi_11 dimension: `13`
- Full-Aut trivial intersection with chi_11 rank: `0`

## Branch chi_11 basis certificate

- Basis index: `0`
- Translation orbit indices: `[0, 65]`
- Stabilizer units: `[1, 11]`
- Values by translation orbit: `{'0': 1, '65': -1}`
- Normalized branch values: `{'A1_k1': -1, 'A5_k5': 1, 'A7_k7': 1, 'A11_k11': -1}`
- Requires premise: choice of chi_11 character / unit-axis; not a full-Aut invariant strict source

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are quotiented by C12 translations into 66 cyclic support orbits.
- `character_projection`: For each U12 character, dim V_chi = (1/4) * sum_u chi(u) * Fix(u) is computed exactly from the residual unit permutation traces.
- `chi11_rank`: For chi_11 with kernel {1,11}, the projection rank is 13.
- `branch_basis`: The A1/A11 and A5/A7 translation orbits form chi_11 basis block 0; unit 11 stabilizes sign pairs while unit 5 and unit 7 flip them.
- `orthogonality_boundary`: Full-Aut invariant data is the trivial character sector; its intersection with the nontrivial chi_11 sector is zero over characteristic zero.
- `strict_limit`: This computes the representation space available after translation quotienting; it does not derive the unit-axis premise needed to choose chi_11 as a strict source.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
- Result is a character-projection certificate on the translation quotient; it does not supply an internal strict source for selecting chi_11.
