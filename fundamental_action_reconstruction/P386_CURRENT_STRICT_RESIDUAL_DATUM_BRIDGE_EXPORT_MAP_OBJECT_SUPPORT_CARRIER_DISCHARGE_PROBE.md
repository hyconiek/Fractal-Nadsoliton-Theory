# P386 Current Strict Residual Datum Bridge/Export-Map Object-Support Carrier Discharge Probe

Status: `P386_EXECUTED_CURRENT_STRICT_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_DISCHARGE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Re-test the strict residual-datum bridge/export-map lane at the post-witness
carrier layer after:

1. `N410` (source-derived contraction parameters), and
2. `N412` (actual provider-object carrier layer),

to determine whether the repo can now export an **actual** post-witness
object-support carrier discharging:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_target_v1 (T140/N407),
```

at least on the provider-object witness track (`N405`), while remaining below
actual bridge/export-map object support above the map object (`N395`) and below
selector closure.

On the updated repo state, the strict-core bridge/export-map object is already
exported (`F311/N422`), so `N300` is historical only and may not be used as a
current-state “below export” clause.

## Probe table (T140 acceptance tests)

| Acceptance test (T140) | Verdict | Evidence |
|---|---|---|
| typed carrier object above witness | YES | `T146/F301` exports `Omicron_residual_datum_bridge_export_map_object_support_carrier_v1` above `N405` |
| route reference discipline | YES | carrier object records `upgrades_witness = Kappa_*_provider_object_*` |
| noncyclic contract (no theta/populated inputs) | YES | finite nad12-depth orbit summary; no theta inputs; no populated-instance inputs |
| observer-free contract | YES | explicit internal `U(1)` gauge; no `K_obs` |
| map-object compatibility | YES | carrier is explicitly scoped below actual object support (`N395`) and does not assert any new export-map object beyond `F311/N422` |
| sigma-int discipline if used | N/A | this carrier is scoped to provider-object witness track only |
| selector neutrality | YES | no `S_sel_int` or `QW-2191` claim |

## Exact verdict

Under the construction specified by `T146`, the repo can export one actual
post-witness object-support carrier object:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_v1,
```

thereby discharging:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_target_v1.
```

This discharge remains explicitly below actual bridge/export-map object support
above the map object (`N395`) and does not imply selector closure or ToE
closure.
