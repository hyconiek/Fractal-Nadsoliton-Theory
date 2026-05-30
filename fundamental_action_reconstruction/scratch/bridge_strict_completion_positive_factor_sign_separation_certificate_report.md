# Strict completion positive-factor sign-separation certificate

Status: `positive-alpha-and-damping-factors-carry-zero-z2-sign-so-completion-sign-is-phase-only`

## Result

This audit proves a finite sign-separation statement for the selected
completion factorization: positive alpha-normalization and damping factors
carry zero Z2 sign bits, so the full completion sign is phase-only.

## Summary

- All alpha-normalization factors positive: `True`
- All damping/compression factors positive: `True`
- All positive-factor Z2 bits zero: `True`
- All factor signs equal phase signs: `True`
- All completion flips equal phase flips: `True`
- Derived sign pattern: `[1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]`
- Derived flip edges: `['1->2', '5->6', '7->8', '9->10']`

## Edge sign separation rows

| edge | positive-factor bit | phase bit | completion bit | phase flip |
|---|---:|---:|---:|---:|
| 0->1 | 0 | 0 | 0 | False |
| 1->2 | 0 | 1 | 1 | True |
| 2->3 | 0 | 0 | 0 | False |
| 3->4 | 0 | 0 | 0 | False |
| 4->5 | 0 | 0 | 0 | False |
| 5->6 | 0 | 1 | 1 | True |
| 6->7 | 0 | 0 | 0 | False |
| 7->8 | 0 | 1 | 1 | True |
| 8->9 | 0 | 0 | 0 | False |
| 9->10 | 0 | 1 | 1 | True |
| 10->11 | 0 | 0 | 0 | False |

## Proof certificate

- `factorization_step`: Use the completion factorization T(d)=A(d)*P(d)*D(d) from the necessity certificate.
- `positivity_step`: The alpha-normalization and damping/compression factors are positive at every audited node; exact damping positivity is separately certified.
- `node_sign_step`: Multiplication by positive factors does not change node sign, so sign(T(d))=sign(P(d)).
- `edge_z2_step`: Positive factors contribute zero Z2 edge bits, so the completion sign coboundary equals the phase sign coboundary.
- `nonduplication`: This is a positive-factor sign-separation certificate, not another phase-zero placement, damping calculus, or real-valued transport fit.
- `theoretical_limit`: The certificate proves finite sign bookkeeping for the selected completion ansatz; it does not derive A(d), P(d), D(d), omega/phi, beta/eta, or transport from strict nadsoliton dynamics.

## Hard limits

- K_strict_gate remains the current live/full operational kernel.
- No unqualified identity K_legacy_ont == K_strict_gate is claimed.
- No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or transport from strict nadsoliton dynamics.
- No beta_tors -> chi_11 theorem is claimed.
- No legacy physical-role transfer to K_strict_gate is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.
