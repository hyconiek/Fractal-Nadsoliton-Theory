# F00 Sandbox Strict-Core Ingredient Candidate Packet

Status: `F00_SANDBOX_STRICT_CORE_INGREDIENT_CANDIDATE_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Build one explicit sandbox-only candidate object that is more ambitious than
another extension lift, but still remains honestly below strict-core
discharge.

## Sandbox candidate

Define one removable sandbox-only candidate object:

```text
G_strict_core_selector_source_sandbox_candidate_v0 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v0,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

with the following intended meaning:

1. `S_sel_int_candidate_seed_v0`
   - reuses the existing strict-side seed carrier from `F36`,
2. `rho_int_orientation_request_slot_v0`
   - one sandbox-only request slot for a future strict-core internal
     orientation-bearing datum,
3. `beta_strict_selector_bridge_request_slot_v0`
   - one sandbox-only request slot for a future strict-core bridge from the
     seed carrier into selector-bearing structure,
4. `lambda_pair1_reachability_request_slot_v0`
   - one sandbox-only request slot for a future strict-core reachability path
     to `A_1(pair1)` or an equivalent selector-facing downstream target.

## Why this is stricter than another extension lift

This object attempts to move in a different direction than `N279..N281`.

Instead of lifting one more admissibility clause in
`strict_extension_only` scope, it introduces one candidate strict-core
construction scaffold that explicitly names the missing positive strict-core
ingredients inside a single object.

## Why this is still not a false pass

`F00` does not claim:

1. that `rho_int_orientation_request_slot_v0` is actual `E_orient`,
2. that `beta_strict_selector_bridge_request_slot_v0` is an actual strict-core
   bridge,
3. that `lambda_pair1_reachability_request_slot_v0` reaches `A_1(pair1)`,
4. that the whole candidate is admissible `S_sel_int`,
5. that strict-core selector closure follows,
6. that ToE closure follows.

## Kernel-split safety

This candidate remains kernel-split-safe because:

1. it does not identify `K_legacy_ont` with `K_strict_gate`,
2. it does not transfer legacy physical-role claims onto `K_strict_gate`,
3. it does not treat `tau_src_candidate_v1` as a hidden strict-core selector
   source.

## Immediate use

The candidate is useful only as one sandbox object to be checked against the
`F29` admission contract axis-by-axis.
