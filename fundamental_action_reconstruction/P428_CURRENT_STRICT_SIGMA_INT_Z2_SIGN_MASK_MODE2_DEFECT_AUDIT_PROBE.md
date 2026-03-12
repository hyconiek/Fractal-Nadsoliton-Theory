# P428 Current Strict Sigma-Int `Z2` Sign‑Mask Mode‑2 Defect Audit Probe

Status: `P428_EXECUTED_CURRENT_STRICT_SIGMA_INT_Z2_SIGN_MASK_MODE2_DEFECT_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

The strict diagonal `pair1` `O(2)`‑cut criterion (`N466`) makes one mode‑2 defect explicit:

```text
diagonal/local sector cuts O(2) on pair1  <=>  F2(d) ≠ 0.
```

Independently, the strict sigma‑int lane exports a strict‑provenance (premise‑based) FR‑derived `Z2` sign mask:

$$
b_{1,k}=(-1)^k,\qquad b_{2,k}=(-1)^{k+1}.
$$

`P428` checks the narrow question (no false pass):

```text
does the exported Z2 parity sign mask itself have a nonzero mode‑2 defect (F2 ≠ 0),
hence could it even in principle act as a diagonal pair1 O(2)-cut ingredient?
```

## Inputs reused

1. exported mask value object:
   - `fundamental_action_reconstruction/generated/b_sigma_int_E_pair_sign_mask_strict_provenance_v1.json`
2. strict diagonal criterion:
   - `N466`

## Result

For both exported masks (`pair1` and `pair2`), the audited mode‑2 coefficient is exactly zero (within tolerance),
while the mode‑6 coefficient is nonzero (the pattern is pure parity / Nyquist mode).

So the parity mask does **not** cut `O(2)` on `pair1` via the strict diagonal criterion.

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p428_current_strict_sigma_int_z2_sign_mask_mode2_defect_audit_probe.json`
- `fundamental_action_reconstruction/generated/p428_current_strict_sigma_int_z2_sign_mask_mode2_defect_audit_probe_summary.json`

## Honest verdict (no false pass)

On the current repo exports:

```text
the strict sigma-int FR-derived Z2 parity mask has F2=0,
therefore it cannot supply a diagonal/local O(2) cut on pair1.
```

Formal theorem closure: `N469`.

