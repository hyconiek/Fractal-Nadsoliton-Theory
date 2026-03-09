# F01 Sandbox Rho-Int Orientation Slot Attack Packet

Status: `F01_SANDBOX_RHO_INT_ORIENTATION_SLOT_ATTACK_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Refine the sandbox orientation slot into something more informative than a
generic placeholder, while still remaining below any real orientation export.

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v1 :=
(
  rho_int_orientation_request_slot_v0,
  residual_orientation_datum_target_slot_export_present,
  future_E_orient_admission_contract_present,
  current_internal_orientation_datum_absent,
  current_sigma_int_residual_bridge_absent,
  current_orientation_export_branch_undischarged
)
```

where:

1. `rho_int_orientation_request_slot_v0`
   - is the original sandbox request slot from `F00`,
2. `residual_orientation_datum_target_slot_export_present`
   - reuses the actual target-slot object from `R1`,
3. `future_E_orient_admission_contract_present`
   - reuses the strict admissibility contract from `F32`,
4. `current_internal_orientation_datum_absent`
   - reuses the strict audit verdict from `B2`,
5. `current_sigma_int_residual_bridge_absent`
   - reuses the route-negative theorem from `N7`,
6. `current_orientation_export_branch_undischarged`
   - reuses the branch-level negative probe from `P135`.

## Updated sandbox candidate

Define the sandbox candidate refinement:

```text
G_strict_core_selector_source_sandbox_candidate_v1 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v1,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This attack does not supply a new strict-core internal orientation datum.

It supplies something narrower:

1. one request slot now aligned to the real strict-core residual-orientation
   target-slot language,
2. one request slot now constrained by the real future `E_orient` admission
   contract,
3. one request slot now annotated by the exact current negative boundaries
   that still block discharge.

So the rho slot is no longer a free-floating placeholder. It is now a
repo-consistent target-slot-bound request scaffold.

## What this still does not claim

`F01` does not claim:

1. actual internal orientation datum,
2. actual `E_orient`,
3. actual residual-datum bridge discharge,
4. actual selector reduction discharge,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. ToE closure.
