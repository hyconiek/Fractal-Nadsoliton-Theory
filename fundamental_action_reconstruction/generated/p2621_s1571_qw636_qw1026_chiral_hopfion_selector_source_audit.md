# P2621/S1571 QW-636/QW-1026 chiral-Hopfion selector source audit

Status: `P2621_CONDITIONAL_CHIRAL_HOPFION_SELECTOR_SOURCE_AUDIT_NO_FULL_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE`

## Professorial verdict on the proposed next step

Scientifically usable only after narrowing: QW-636 and QW-1026 are heuristic/numerical prior-art scripts, not strict proofs by themselves. They can be formalized into a conditional theorem exporting an orientation-odd source premise, but they do not unconditionally discharge QW-2191 or repair the full P2620 bridge.

## Anti-duplication grep audit

- `new_packet`: 24 hits; samples retained in JSON certificate.
- `qw636_prior_art`: 67 hits; samples retained in JSON certificate.
- `qw1026_prior_art`: 41 hits; samples retained in JSON certificate.
- `selector_chirality_prior_art`: 7967 hits; samples retained in JSON certificate.
- `bridge_guards`: 6280 hits; samples retained in JSON certificate.

## Theorem export

**Claim.** If the strict kernel is supplied with a typed nonzero chiral Hopfion/anomaly source, then sigma=sign(Re Tr(gamma5 K^3)) (equivalently the normalized QW-636 parity-odd Hopfion energy delta) is a C2-odd topological phase selector.  This repairs the P2619/P2620 orientation atom only under that explicit source premise.

Proof skeleton:
- For QW-636-style directional Hermitian hopping t_sigma(r)=a_r exp(i sigma theta), E_sigma(k)-E_sigma(-k)=-4 sigma sin(theta) sum_r a_r sin(k r).
- When sin(theta) and the declared sine probe are nonzero, the normalized parity-odd energy delta recovers sigma exactly and flips under orientation reversal.
- For QW-1026-style anomaly A_sigma=Tr(gamma5 K_sigma^3), a typed chiral kernel satisfying Re A_-sigma=-Re A_sigma and Re A_sigma != 0 yields sigma=sign(Re A_sigma) after an orientation convention.
- The C2 enumeration verifies that a freely acting chiral two-point torsor admits equivariant sign maps, unlike the invariant scalar inputs blocked by P2619.

## Computational certificates

- Hopfion energy selector recovers sigma on all rows: `True`.
- Chiral anomaly selector recovers sigma on all rows: `True`.
- C2 equivariant sign maps from the chiral torsor: `2`.
- P2620 conditional accepting rows: `1` of 4.

## Guarded conclusion

Orientation atom status: `CONDITIONALLY_REPAIRED_IF_TYPED_NONZERO_CHIRAL_SOURCE_IS_ADMITTED`.

This does **not** repair the nonlinear damping completion atom and therefore does **not** by itself reopen role-bearing `L_total`, role transfer, QW-2191 discharge, or ToE closure.

Not licensed:
- unconditional strict-core selector closure
- nonlinear damping completion source
- full P2620 bridge-source cut repair without the damping atom
- GF(2) bridge revalidation
- beta_tors -> chi11 route reopening
- role-bearing L_total re-enable
- legacy physical-role transfer
- QW-2191 discharge
- ToE closure

Certificate hash: `51750ca6e6faf463719a5eed63ad851544b710a6e007a36ee787164faf800063`
