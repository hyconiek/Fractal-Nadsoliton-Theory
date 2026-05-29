# Scratch strict-alpha phase Fourier selector discriminator probe

Status: Fourier selector discriminator for eta=9/5; no strict selector theorem.

- Exact Parseval discriminator: `total_non_dc_power = N*sum(e_i^2) - (sum e_i)^2`.
- Minimum-ripple winner: `[2, 2, 2, 1, 1]`; this matches the balanced ledger.
- Maximum total-resonance-power winner: `[4, 1, 1, 1, 1]`; this is not the balanced ledger.
- Bounded Z12 placement scan: every canonical ledger can reach single-mode power `64`, so naive highest single-mode resonance is placement-dominated.
- Honest read: phase/Fourier language supports `(2,2,2,1,1)` only as a minimum-ripple/coherence-compression selector, not as an automatic highest-resonance selector.
- No false pass: no phase/Fourier strict selector theorem, no QW-2191 discharge, no ToE closure.
