# N497 Current First Strict Shannon vs Diagonal Mode‑Index Assignment Axis Alignment Up to Residual `Z2` (Value‑Instantiation) Theorem (No False‑PASS)

Status: `N497_DISCHARGED_CURRENT_FIRST_STRICT_SHANNON_VS_DIAGONAL_MODE_INDEX_ASSIGNMENT_AXIS_ALIGNMENT_UP_TO_RESIDUAL_Z2_VALUE_INSTANTIATION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Two independent strict exports now provide explicit mode‑index assignment basis objects on the same strict `n=12` Fourier scaffold:

1. diagonal/local strict-derived assignment:
   `ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1` (`F453`), and
2. Shannon element‑order reference strict-core assignment:
   `ModeIndexAssignment_shannon_element_order_reference_strict_core_v1` (`F454`).

This theorem packages one narrow **value‑instantiation** conclusion:

```text
on the current exported n=12 instantiations, the two assignments select the same axes
inside each pair_m plane (m=1..5), up to the unavoidable residual Z2 sign.
```

This does **not** claim any global discharge of `QW-2191`, does **not** claim strict‑core selector closure, and does **not** claim ToE closure.

## Strict-admissible inputs reused

1. `F453`
   - exported diagonal/local mode-index assignment basis object:
     `fundamental_action_reconstruction/generated/mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json`.
2. `F454`
   - exported Shannon element‑order reference mode-index assignment basis object:
     `fundamental_action_reconstruction/generated/mode_index_assignment_shannon_element_order_reference_strict_core_v1.json`.
3. `P455`
   - mechanical alignment audit:
     `fundamental_action_reconstruction/generated/p455_current_strict_mode_index_assignment_shannon_vs_diagonal_alignment_audit_probe_summary.json`.

## Theorem (value‑instantiated axis alignment on current exported objects)

On the current repo state, `P455` reports:

1. `aligned_all_pairs_up_to_residual_z2 = true`,
2. `max_theta_star_mod_pi_abs_diff ≈ 4.44e-16`,
3. `min_primary_abs_dot_across_pairs ≈ 1.0`,
4. `max_cross_abs_dot_across_pairs ≈ 7.44e-16`.

Therefore, on the exported strict `n=12` instantiations, the diagonal/local and Shannon element‑order reference mode‑index assignments agree on the **axis choice** on every Fourier-degenerate pair plane `pair_m (m=1..5)` up to residual `Z2` sign (basis-vector sign flips).

Equivalently:

```text
for each m=1..5, span{u_{m,+}^diag, u_{m,-}^diag} = span{u_{m,+}^sha, u_{m,-}^sha}
and the ordered axes match up to sign flips.
```

## Meaning (no false‑PASS)

This theorem means only:

1. the axis selection is not an artifact of a single lane: two independent strict constructions converge on the same axes (up to sign) on `n=12`,
2. residual `Z2` sign remains a convention unless separately proven gauge‑irrelevant for a given downstream observable.

## What N497 does not claim

`N497` does not claim:

1. a theorem identifying the diagonal/local residual profile with `ord_Z12(x)` (or any stronger bridge),
2. strict-core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.

