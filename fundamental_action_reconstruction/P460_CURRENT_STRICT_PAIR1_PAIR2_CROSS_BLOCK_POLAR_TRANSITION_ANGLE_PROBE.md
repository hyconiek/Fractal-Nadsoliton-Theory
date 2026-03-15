# P460 Current Strict Pair1/Pair2 Cross-Block Polar Transition Angle Probe

Status: `P460_EXECUTED_CROSS_BLOCK_POLAR_TRANSITION_ANGLE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Given the exported **value-instantiated** declared residual local-diagonal control pullback on the control basis
`(c1,s1,c2,s2)` (`P459`),
extract the `pair1 -> pair2` cross-block and compute a canonical **polar orthogonal factor**:

```text
B = Q S   (polar decomposition),
Q ∈ O(2).
```

When `det(Q)=+1`, interpret `Q` as a rotation and export the corresponding lane-scoped transition angle:

```text
alpha_cross := atan2(Q[1,0], Q[0,0])  (mod 2π).
```

This produces a concrete, testable candidate for the missing “transition object” data in the control-lane sense, while
keeping explicit that this is:

- **not** a global selector atlas/gluing object,
- **not** a strict-core selector-closure discharge,
- and depends on the conditional `N477` rewrite used to instantiate `P459`.

## Inputs

1. `M_control_residual_diag_declared_value_instantiated_v1` (from `P459`):
   - 4×4 symmetric matrix on `(c1,s1,c2,s2)`,
   - contains the explicit 2×2 cross-block `B`.
2. `Alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1` (from `F457`):
   - lane-scoped `alpha_12 := (theta_2-theta_1) mod 2π`,
   - used only for comparison/hygiene, not for derivation.

## Output artifacts

- `fundamental_action_reconstruction/generated/p460_current_strict_pair1_pair2_cross_block_polar_transition_angle_probe.json`
- `fundamental_action_reconstruction/generated/p460_current_strict_pair1_pair2_cross_block_polar_transition_angle_probe_summary.json`

## Hard limits (no false pass)

This probe does **not** claim:

1. that `alpha_cross` is a physical mixing angle,
2. that any global selector transition/gluing object exists,
3. that `QW-2191` is discharged,
4. that host matching/cancellation is achieved,
5. that any strict-core selector closure / admissible `S_sel_int` is exported.

