# F461 Current Strict Sigma-int Pair1/Pair2 Selector Chart Transition Object Export Packet

Status: `F461_DRAFT_EXPECTED_EXECUTED_BY_f461_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Export one **explicit, lane-scoped** chart-transport / gluing **operator-level** object that acts **between**
the two strict sigma-int corridor local selector charts (`pair1` and `pair2`) and transports the exported
pair representatives in a **sign-gauge-safe** (projective) manner.

This packet is strictly scoped:

- carrier: `n=12` real Fourier scaffold (`QW-2190`),
- charts: `pair1=(c1,s1)` and `pair2=(c2,s2)` only,
- source: strict sigma-int slot-free theta-pair supply (`F451/N489`) and its derived transition angle export (`F457`),
- output: an explicit orthogonal operator `O_12` transporting `pair1` ↔ `pair2` while leaving the orthogonal complement fixed.

It does **not** claim a global selector atlas/gluing structure (`H41`) nor a global selector transition object
discharging `QW-2191` (`H40` remains globally open).

## Strict inputs

1. `ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1` (`F451`)
   - provides `u_1 ∈ span{c1,s1}` and `u_2 ∈ span{c2,s2}` (declared scope) and their angles `theta_1, theta_2`.
2. `Alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1` (`F457`)
   - exports `alpha_12 := (theta_2 - theta_1) mod 2π` and the 2×2 rotation `G(alpha_12)`.
3. `QW-2190` real Fourier scaffold (as a declared carrier structure)
   - provides canonical `c_m,s_m` pair planes for `n=12` used to build the operator embedding.

## Construction (definition)

Let:

- `V1 := span{c1,s1}` and `V2 := span{c2,s2}` in the `n=12` carrier,
- `C1 := [c1 s1]` and `C2 := [c2 s2]` (12×2 matrices with orthonormal columns),
- `Π1 := C1 C1^T`, `Π2 := C2 C2^T`, `Π_rest := I - Π1 - Π2`,
- `G(α) ∈ SO(2)` the standard 2×2 rotation by `α := alpha_12 (mod 2π)`.

Define the explicit 12×12 orthogonal chart-transport operator:

```text
O_12 := C2 G(α) C1^T  +  C1 G(-α) C2^T  +  Π_rest.
```

Then, in the declared sigma-int corridor scope, `O_12` transports the exported pair representatives:

```text
O_12 u_1 ≈ u_2
P(u_2) ≈ O_12 P(u_1) O_12^T,  where P(u) := u u^T.
```

Crucially, the projector transport is **sign-gauge-invariant**: it is unchanged under `u -> -u`.

## Exported artifact

This packet exports:

- `O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1`  
  in `fundamental_action_reconstruction/generated/o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json`
  with a `_summary.json` audit.

## Hard limits (no false pass)

This packet does **not** claim:

1. any global selector atlas / overlap-domain declaration / cocycle data (`H41` remains open),
2. any global selector transition/gluing object discharging `QW-2191` (`H40` remains open globally),
3. any sign-sensitive physical orientation datum (residual `Z2` remains a convention/gauge layer),
4. any physical interpretation of `alpha_12` as a fundamental mixing angle,
5. strict-core selector closure / admissible `S_sel_int`,
6. ToE closure.

