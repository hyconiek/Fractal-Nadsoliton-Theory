# D12 chi_11 sparsest-extension selector audit

Status: `sparsest-branch-normalized-d12-chi11-extension-is-unique-but-sparsity-is-extra`

## Finite model

- Ring: `Z_12`
- Enumerated supports: `792`
- D12 orbit count: `38`
- chi_11 two-cycle count: `13`

## Sparsity summary

- Branch cycle orbit pair: `[0, 37]`
- Free nonbranch cycle count: `12`
- Branch-normalized ternary functions: `531441`
- Minimum D12-component support size: `2`
- Minimum-support function count: `1`
- Unique minimum is branch generator: `True`

## Support-size distribution

- support `2`: `1` functions (nonzero free coefficients `0`)
- support `4`: `24` functions (nonzero free coefficients `1`)
- support `6`: `264` functions (nonzero free coefficients `2`)
- support `8`: `1760` functions (nonzero free coefficients `3`)
- support `10`: `7920` functions (nonzero free coefficients `4`)
- support `12`: `25344` functions (nonzero free coefficients `5`)
- support `14`: `59136` functions (nonzero free coefficients `6`)
- support `16`: `101376` functions (nonzero free coefficients `7`)
- support `18`: `126720` functions (nonzero free coefficients `8`)
- support `20`: `112640` functions (nonzero free coefficients `9`)
- support `22`: `67584` functions (nonzero free coefficients `10`)
- support `24`: `24576` functions (nonzero free coefficients `11`)
- support `26`: `4096` functions (nonzero free coefficients `12`)

## Proof certificate

- `finite_domain`: All C(12,5)=792 supports are quotiented by D12; the residual unit-5 action gives the chi_11 two-cycle basis.
- `module_reuse_boundary`: The D12 chi_11 module has 13 two-cycle basis generators; this probe only studies branch-normalized ternary extensions inside that already computed module.
- `branch_constraint`: A1/A11=-1 and A5/A7=+1 fixes the branch two-cycle coefficient to +1.
- `free_coefficients`: The other 12 two-cycle coefficients are independently free in (-1, 0, 1), giving 3^12=531441 branch-normalized ternary extensions.
- `support_formula`: If k non-branch coefficients are nonzero, the D12-component support size is 2+2k.
- `unique_minimum`: The unique minimum-support extension has k=0, support size 2, and equals the branch generator only.
- `strict_limit`: The sparsest-extension rule is an added selector criterion, not a strict-core derivation of chi_11 or QW-2191 discharge.

## Hard limits

- No identity K_legacy_ont == K_strict_gate is asserted.
- No legacy physical-role transfer onto K_strict_gate is used.
- No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.
- No theorem derives a sparsest-extension selector from strict core.
- No QW-2191 discharge is claimed.
- No ToE closure is claimed.
