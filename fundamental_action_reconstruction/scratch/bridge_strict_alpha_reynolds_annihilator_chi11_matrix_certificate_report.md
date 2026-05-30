# Reynolds annihilator chi_11 matrix certificate

Status: `exact-integer-matrix-annihilator-computed-on-translation-quotient`

## Finite model

- Ring: `Z_12`
- Enumerated supports: `792`
- Translation orbit count / matrix size: `66` / `66`
- Residual unit orbit count: `25`
- Residual unit orbit size counts: `{1: 4, 2: 11, 4: 10}`

## Exact matrix certificate

- Reynolds numerator rank: `25`
- chi_11 numerator rank: `13`
- R_num*C_chi11_num zero: `True`
- C_chi11_num*R_num zero: `True`
- Zero entries in annihilator product: `4356/4356`

## Branch generator witness

- Quotient pair indices: `[0, 65]`
- Normalized branch values: `{'A1_k1': -1, 'A11_k11': -1, 'A5_k5': 1, 'A7_k7': 1}`
- A1 representative support: `[0, 1, 2, 3, 4]`
- A5 representative support: `[0, 2, 4, 7, 9]`
- Reynolds numerator kills branch: `True`
- chi_11 numerator returns 4*branch: `True`

## chi_11 basis annihilation

- chi_11 basis rows: `13`
- Every row has zero Reynolds sum: `True`

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are first quotiented by translations into a 66-point exact finite representation.
- `matrix_construction`: For each unit u in U12, P_u is the exact 66x66 permutation matrix on translation orbits; R_num=sum_u P_u and C_chi11_num=sum_u chi_11(u)P_u.
- `annihilator_identity`: The computed integer products R_num*C_chi11_num and C_chi11_num*R_num are zero matrices, so full-Aut Reynolds data annihilates the chi_11 sector exactly.
- `branch_witness`: The compact A1/A11=-1 and A5/A7=+1 translation vector is a chi_11 eigenvector with C_chi11_num*v=4v and R_num*v=0.
- `logical_boundary`: This strengthens the obstruction to full-Aut export of chi_11 polarity, but the chi_11/unit-axis choice remains an extra premise, not a strict-core theorem.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No theorem derives a full-Aut internal chi_11 polarity source.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
