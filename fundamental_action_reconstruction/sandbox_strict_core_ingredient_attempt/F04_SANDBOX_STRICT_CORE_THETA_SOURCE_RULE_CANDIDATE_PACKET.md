# F04 Sandbox Strict-Core Theta-Source Rule Candidate Packet

Status: `F04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_PACKET_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

Upgrade the sandbox from one non-placeholder strict-core theta-source
skeleton attempt to one explicit conditional rule candidate, while preserving
the `C50` blocker.

## Conditional strict-core theta-source rule candidate

Define:

```text
theta_source_rule_candidate_v0 :=
(
  populated_basis_pair_instance_if_supplied_then_phase_serialization_defined,
  no_actual_populated_basis_pair_instance_exported,
  no_actual_theta_values_exported
)
```

with the following meanings:

1. `populated_basis_pair_instance_if_supplied_then_phase_serialization_defined`
   - reuse `C33` and `C49`:

   ```text
   if populated_instance(theta_1,theta_2) contains actual u_1,u_2,
   then
     theta_1 = atan2(<s_1,u_1>, <c_1,u_1>),
     theta_2 = atan2(<s_2,u_2>, <c_2,u_2>)
   ```

   In words:
   once an actual populated pair instance is present, phase serialization is
   rule-defined by the strict-core local phase formulas.

2. `no_actual_populated_basis_pair_instance_exported`
   - preserve `C49/C50`:
     the conditional schema is packet-ready, but no actual populated
     `u_1,u_2` instance is exported from strict core.

3. `no_actual_theta_values_exported`
   - preserve `C50`:
     no actual strict-core `theta_1`, `theta_2` values are exported.

## Refined rho slot

Define:

```text
rho_int_orientation_request_slot_v4 :=
(
  rho_int_orientation_request_slot_v3,
  theta_source_rule_candidate_v0
)
```

## Updated sandbox candidate

Define:

```text
G_strict_core_selector_source_sandbox_candidate_v4 :=
(
  S_sel_int_candidate_seed_v0,
  rho_int_orientation_request_slot_v4,
  beta_strict_selector_bridge_request_slot_v0,
  lambda_pair1_reachability_request_slot_v0
)
```

## Meaning

This attack does not create a packet-ready strict-core minimal source
skeleton.

It creates something narrower:

1. one explicit conditional rule candidate showing how actual populated
   `u_1,u_2` would serialize phases if strict core ever supplied such an
   instance,
2. one preserved record that strict core still does **not** supply that
   populated instance,
3. one preserved record that strict core still does **not** export actual
   `theta_1`, `theta_2`.

So the sandbox no longer stops at a skeleton attempt.
It now contains one strict-core-only conditional theta-source rule candidate.

## What this still does not claim

`F04` does not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. a packet-ready strict-core minimal source skeleton in the strong sense
   denied by `C50`,
4. actual internal orientation datum,
5. actual `E_orient`,
6. admissible `S_sel_int`,
7. strict-core selector closure,
8. ToE closure.
