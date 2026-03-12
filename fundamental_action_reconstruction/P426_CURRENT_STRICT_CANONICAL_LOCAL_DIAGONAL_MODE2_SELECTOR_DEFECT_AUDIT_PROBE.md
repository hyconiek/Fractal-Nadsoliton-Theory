# P426 Current Strict Canonical Local-Diagonal Mode‑2 Selector‑Defect Audit Probe

Status: `P426_EXECUTED_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_MODE2_SELECTOR_DEFECT_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `N465/N466/P424/P425`, the strict picture around “real physical accelerator of choice” is:

1. the certified host operator `A = K_total + m0^2 I` is **isotropic** on `pair1 = span{c1,s1}` and therefore cannot
   cut the `O(2)` family (`N465`),
2. a diagonal/local sector can cut `O(2)` on `pair1` **iff** its diagonal profile has a nonzero mode‑2 Fourier
   coefficient `F2(d)` (`N466`).

`P426` asks the narrowest constructive question, with no false pass:

```text
Can we at least export the exact canonical-diagonal mode‑2 defect F2(d),
in a form reduced to the actually exported six opposite-pair sums from R18,
so the strict-core O(2)-cut question becomes a single checkable expression?
```

This probe does **not** claim that the defect is nonzero.

## Inputs reused

1. `R15`
   - exports the canonical local diagonal sector and its residual decomposition
     `D_canonical = m0^2 I + D_local_residual`.
2. `R18`
   - exports the six opposite-pair residual sums on the 12-slot carrier:
     `Sigma_psi0_psi6, ..., Sigma_psi5_psi11`.
3. `N466`
   - mathematical criterion: diagonal `pair1` `O(2)`-cut iff `F2(d) ≠ 0`.

## Result

`P426` exports a persisted artifact computing:

1. the reduced six-sum form of `F2(d)` for `n=12`:
   `F2(d) = Σ_{k=0..5} S_k exp(i*pi*k/3)` with `S_k := d_k + d_{k+6}`,
2. the explicit `Re(F2)` and `Im(F2)` linear combinations of the six exported `R18` residual sums,
3. the induced `pair1` diagonal `O(2)`-cut signature:
   `(a1-d1, b1) = ((1/6)Re(F2), (1/12)Im(F2))`.

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p426_current_strict_canonical_local_diagonal_mode2_selector_defect_audit_probe.json`
- `fundamental_action_reconstruction/generated/p426_current_strict_canonical_local_diagonal_mode2_selector_defect_audit_probe_summary.json`

## Honest verdict (no false pass)

On the current repo exports, the defect is **symbolic** (in terms of `g4_psi*`, `g6_psi*`, `gY*`, `m2_psi*`, `vpsi*`,
`vphi`), and no strict-derived coefficient instantiation exists. Therefore:

```text
P426 exports an exact reduced defect expression,
but it does not yield a strict-core O(2)-cut witness or any QW-2191 discharge.
```

## Recommended next move (strict)

To turn this into a strict-core `O(2)`-cut witness, the repo must export a strict-derived
non-translation-invariant diagonal/local profile (or equivalently strict-derived values sufficient to decide the
nonvanishing of the exported `F2` defect).

