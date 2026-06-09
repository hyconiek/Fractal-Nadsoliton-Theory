# P2622/S1572 QW-636/QW-1026 physical-rigor nonpromotion audit

Status: `P2622_QW636_QW1026_PRIOR_ART_NONPROMOTION_NO_SELECTOR_SOURCE_NO_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE`

## Anti-duplication grep audit

- `new_packet`: 24 hits; samples retained in JSON certificate.
- `qw636_qw637_prior_art`: 260 hits; samples retained in JSON certificate.
- `qw1026_prior_art`: 90 hits; samples retained in JSON certificate.
- `selector_guard_prior_art`: 8170 hits; samples retained in JSON certificate.
- `methodology_risk_terms`: 1013 hits; samples retained in JSON certificate.

## Theorem export

**Claim.** QW-636 and QW-1026 contain useful parity/chirality diagnostics, but as written they do not export the missing strict orientation_odd_selector_source.  QW-636's phase asymmetry is gauge/covariance sensitive and needs an already-directed momentum or loop orientation.  QW-1026's sign is tied to an unsourced alternating gamma5 convention.  Therefore P2621 remains only a conditional schema, not an unconditional source theorem.

Positive retained content:
- QW-636 can be retained as a diagnostic for a typed chiral/flux source once a gauge-invariant orientation source is supplied.
- QW-1026 can be retained as a diagnostic for a typed chirality operator once gamma5 is physically sourced and basis-safe.
- P2621 remains a conditional implication: typed nonzero chiral source implies an odd selector.

Obstructions:
- QW-636: E_sigma(k)=E_-sigma(-k), sorted spectra for sigma=+/- match, and open-chain Peierls phases are pure gauge.
- QW-1026: gamma5=diag((-1)^i) changes sign under the alternating-origin choice, flipping Re Tr(gamma5 K^3).
- Neither script supplies the independent nonlinear damping completion atom required by P2620.

## Computational certificates

- QW-636 parity covariance max defect: `0.0`.
- QW-636 sorted sigma-spectrum defect: `1.7763568394002505e-15`.
- QW-636 open-chain gauge-equivalence defect: `1.1102230246251565e-16`.
- QW-1026 opposite gamma5 sign sum real/imag: `0.0`, `0.0`.
- Accepted orientation source from QW-636/QW-1026 prior art alone: `False`.
- P2620 accepting rows from prior art alone: `0`.

## Next admissible targets
- derive a gauge-invariant Wilson-loop/flux orientation source on a closed typed cycle with fixed orientation convention sourced internally
- derive gamma5/chirality as a basis-independent operator from strict nadsoliton field content rather than alternating labels
- continue nonlinear damping completion separately; selector progress alone cannot repair P2620

Not licensed:
- promotion of QW-636/QW-1026 prior art alone to strict selector source
- unconditional P2621 source closure
- P2620 bridge-source cut repair
- role-bearing L_total
- QW-2191 discharge
- ToE closure

Certificate hash: `39c34b0f12e9fc61a12aef9bb86f2b865c41e66e3116c94fccb330976320c253`
