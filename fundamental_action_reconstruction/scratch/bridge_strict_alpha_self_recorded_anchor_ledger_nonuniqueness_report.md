# Scratch strict-alpha self-recorded anchor ledger nonuniqueness probe

Status: endpoint self-record nonuniqueness audit; no strict selector discharge.

- Positive ordered ledgers: `35`; endpoint-distinct self-recording: `22`; endpoint-equal ambiguous: `13`.
- Canonical partition counts: `{'4,1,1,1,1': 5, '3,2,1,1,1': 20, '2,2,2,1,1': 10}`.
- Honest read: endpoint self-recording is a decoder/consistency property, not a unique ledger selector.
- Target replay: `q^5=256/243`, eta residual `0.000e+00`.
- No false pass: no strict ledger theorem, no QW-2191 discharge, no ToE closure.
