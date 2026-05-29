# Scratch strict-alpha Hebbian prospective horizon value-iteration no-go probe

Status: finite-horizon full-Aut-equivariant lookahead cannot create singleton d5 from invariant future value.

- Supports scanned: `792`; histogram classes: `35`.
- Equivariant kernels audited: `13` with `T_a = [[a, 1-a], [1-a, a]]`.
- Future reward grid: `169` rewards; invariant rewards: `13`; invariant+d5-dominant rewards: `0`.
- Horizons checked: `0..8`; invariant reward tie violations: `0`.
- Controlled invariant Bellman trace remains tied: `True`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: highest-future-resonance selects d5 only if the future value already names d5.
- No false pass: no future-value source theorem, no QW-2191 discharge, no ToE closure.
