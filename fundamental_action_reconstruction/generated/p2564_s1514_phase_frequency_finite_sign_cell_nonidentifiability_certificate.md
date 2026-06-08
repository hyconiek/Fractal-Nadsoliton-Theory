# P2564/S1514 phase/frequency finite sign-cell nonidentifiability certificate

Status: `P2564_PHASE_FREQUENCY_FINITE_SIGN_CELL_NONIDENTIFIABILITY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_phase_frequency_source`.
- Source class audited: `finite d=0..11 phase-sign/topological-bit pattern as selector for exact omega and phi`.
- P2563 rational-winding obstruction inherited: `True`.
- P2561 phase/frequency residual atom inherited: `True`.
- P2415 affine nonautomorphism inherited: `True`.
- Strict sign pattern d=0..11: `[1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1]`.
- Certified epsilon omega: `0.0032376530502126466`.
- Certified epsilon phi: `0.0032376530502126466`.
- Open sign-cell box area: `4.1929589094205015e-05`.
- Grid witnesses preserving signs: `25/25`.
- Nontrivial same-sign witnesses: `24`.
- Finite sign pattern selects unique tuple: `False`.

## Proof interpretation

The strict finite sign pattern has a positive clearance from every audited cosine zero.  A small open box around `(omega, phi)=(743/4000,13/80)` therefore preserves all audited signs, so finite GF(2)/phase-sign data cannot by itself select the exact strict phase/frequency tuple.

## Recommended next honest step

Do not treat GF(2) or finite phase-sign reconstruction as a source for the exact strict phase/frequency tuple. The next honest step is to add a stronger selector candidate with a computable objective inside the open sign cell, then test whether it has a unique minimizer at omega=743/4000 and phi=13/80 without discharging QW-2191 by assumption.

## Negative controls

No finite-sign phase/frequency source, strict phase/frequency source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`a1028e734ffb5f6180f753fffc90e10a88da63c80f4d3f872e8daf87bbe0bda2`
