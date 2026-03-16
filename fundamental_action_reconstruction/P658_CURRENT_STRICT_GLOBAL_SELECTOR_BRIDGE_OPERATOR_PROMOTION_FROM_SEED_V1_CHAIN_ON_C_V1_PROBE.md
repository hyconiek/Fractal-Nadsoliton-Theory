# P658 Current Strict Global Selector Bridge Operator Promotion From Seed‑v1 Chain On `C_v1` Probe (No False‑PASS)

Status: `P658_CURRENT_STRICT_GLOBAL_SELECTOR_BRIDGE_OPERATOR_PROMOTION_FROM_SEED_V1_CHAIN_ON_C_V1_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Audit that `F658` really exported a **global** `C_v1`‑typed selector bridge operator family promoted from the seed‑v1 local `B_sel` on `pair1`,
and that it is glued consistently (projector/section level) by the already exported global transition data, while keeping strict-core closure claims explicitly out.

## Inputs

1. `F658`
   - exported global object:
     `SelectorBridgeOperator_global_C_v1_seed_v1_promoted_strict_v1`.
2. `F469/N515`
   - global selector transition/gluing object on `C_v1` (`T170` discharge).
3. `F470/N516`
   - global projective selector state object on `C_v1` (projector/span semantics).
4. `F655/N547`
   - seed‑v1 local `B_sel` on `pair1`.
5. `N512`
   - boundary: no operator-level transition groupoid promotion on the full carrier.

## Checks (acceptance)

The probe verifies:

1. each chartwise `B_sel(pair_m)` is symmetric and involutive (`B^2 = I`) on the declared pair plane,
2. overlap transport consistency: for every exported transition edge `pair_i -> pair_j` with SO(2) transport `G_ij`,
   the exported matrices satisfy:
   `B_j = G_ij B_i G_ij^T` within tolerance,
3. pair1 alignment: the global chartwise `B_1` matches the seed‑v1 local export (`F655`) within tolerance,
4. projective consistency: the induced `P_plus = (I+B)/2` matches the already exported local projectors `A_m` of the
   global projective selector state (within tolerance),
5. strict discipline flags remain explicit: no selector-closure / no `QW-2191` discharge claims.

## Output

Executed by:

```text
python3 fundamental_action_reconstruction/p658_current_strict_global_selector_bridge_operator_promotion_from_seed_v1_chain_on_c_v1_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p658_current_strict_global_selector_bridge_operator_promotion_from_seed_v1_chain_on_c_v1_probe_summary.json`

