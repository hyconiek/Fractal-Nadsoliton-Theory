# F298 First Actual Strict ToE Closure Provider-Object Carrier Contraction-Parameter Source Map Packet

Status: `F298_EXECUTED_FIRST_ACTUAL_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_CONTRACTION_PARAMETER_SOURCE_MAP_PACKET_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

After `T142/N409`, the provider-object carrier lane has one explicit upstream
missing ingredient class named sharply:

```text
a strict source-side derivation/source-upgrade for the contraction parameters (a,b)
```

The next honest move is to attempt a narrow discharge that does not touch:

1. provider-object realization,
2. residual bridge/export-map object support (`N302`),
3. selector closure.

`F298` executes exactly that narrow move:

```text
export one explicit source-side map tau_src_candidate_v1 -> (a,b)
derived only from already exported strict source topology control data.
```

## Inputs reused

1. `F127/N235`
   - `tau_src_candidate_v1` exists as an exported strict source packet,
2. `F74`
   - the strict-kernel core flow amplitude is carried as:
     `T_flow^(0) := cos(phi) * e_topo`,
3. `F141/N249`
   - barrier-protected sign and margin witnesses ensure `0 < cos(phi) < 1` on
     the declared core branch,
4. `T143`
   - contraction-parameter source-map spec.

## Packet result (actual map object)

`F298` exports one actual packaged source-side map:

```text
A_strict_provider_object_contraction_parameter_source_map_v1 :
  tau_src_candidate_v1 -> (a,b)
```

defined by:

```text
A_strict_provider_object_contraction_parameter_source_map_v1(tau_src_candidate_v1)
  := (cos(phi), cos(phi)).
```

## Target discharge statement

This map satisfies the acceptance tests of:

```text
Delta_strict_provider_object_carrier_contraction_parameter_strict_derivation_source_upgrade_target_v1
```

from `T142`, on the declared strict core branch.

## Status discipline

This packet does **not** claim:

1. discharge of `Epsilon_strict_provider_object_carrier_layer_target_v1`,
2. discharge of `N302` (actual residual bridge/export-map object support),
3. any export-map object export (`N300`),
4. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
5. ToE closure.

It claims only:

1. the provider-object carrier contraction parameters `(a,b)` are now exported
   as source-side strict-derived outputs of `tau_src_candidate_v1` via one
   explicit map object,
2. this removes one upstream “free knob” from the provider-object carrier lane,
   without touching the downstream bridge/export-map object-support frontier.

