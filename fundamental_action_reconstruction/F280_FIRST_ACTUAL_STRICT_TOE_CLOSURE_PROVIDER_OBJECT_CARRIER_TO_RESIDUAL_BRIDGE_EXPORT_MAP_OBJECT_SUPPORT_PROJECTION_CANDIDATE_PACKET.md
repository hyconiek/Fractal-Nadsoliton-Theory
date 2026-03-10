# F280 First Actual Strict ToE Closure Provider-Object Carrier to Residual Bridge/Export-Map Object-Support Projection Candidate Packet

Status: `F280_EXECUTED_FIRST_ACTUAL_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_TO_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PROJECTION_CANDIDATE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

After `N302`, the residual bridge/export-map route is blocked below **actual**
object support.

After `N385`, the repo exports one projection layer into the object-support
frontier, but still remains below actual object support.

`F280` does **not** claim `N302` is discharged.

`F280` executes a weaker strict move:

```text
export one explicit bridge-facing projection *candidate*
from a provider-object carrier candidate (T126/F279)
into the residual bridge/export-map object-support-facing record,
using only noncyclic inputs and no K_obs as primary source.
```

## Inputs reused

1. `T126/F279/N391`
   - one explicit provider-object carrier orbit-quotient candidate is exported,
2. `T127`
   - projection candidate spec from provider-carrier into residual bridge
     frontier,
3. `T115/F268/N380`
   - strict phasor reduction pair-map-rule candidate form exists,
4. `N302`
   - boundary below actual object support remains in force.

## Candidate projection result

`F280` exports one actual packaged projection candidate:

```text
Pi_strict_provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_v1
```

with the intended scoped meaning:

```text
ProviderObjectCarrier_pair^{cand,orbit}
  -> E_pair^{cand,prov} (finite, noncyclic carrier)
  -> (theta_1^cand, theta_2^cand) via T115
  -> residual bridge/export-map object-support-facing projection record (candidate)
```

and with explicit guardrail fields:

```text
theta_inputs_used = false
populated_instance_inputs_used = false
K_obs_used_as_primary_source = false
actual_bridge_export_map_object_support_present = false
admissible_S_sel_int_present = false
strict_core_selector_closure_present = false
QW_2191_discharged = false
ToE_closed = false
```

## Persisted candidate artifact instance

`F280` also records one persisted candidate artifact instance:

```text
fundamental_action_reconstruction/generated/provider_object_carrier_to_residual_bridge_export_map_object_support_projection_candidate_instance.json
```

This file is a carrier artifact only. It is not an export-spec and not a
discharge.

## What F280 does not claim

`F280` does not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. actual export-map object export,
4. actual theta export,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. `QW-2191` discharge,
8. ToE closure.

