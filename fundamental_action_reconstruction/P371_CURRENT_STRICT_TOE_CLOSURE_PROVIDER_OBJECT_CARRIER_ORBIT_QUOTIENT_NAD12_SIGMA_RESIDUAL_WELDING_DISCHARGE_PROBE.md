# P371 Current Strict ToE Closure Provider-Object Carrier Orbit-Quotient ↔ Nad12-Sigma Residual Welding Discharge Probe

Status: `P371_EXECUTED_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_ORBIT_QUOTIENT_NAD12_SIGMA_RESIDUAL_WELDING_DISCHARGE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

Confirm whether the current repo exports a theorem-facing discharge witness for:

```text
Zeta_strict_provider_object_carrier_orbit_quotient_nad12_sigma_residual_welding_target_v1
```

as specified by `T128`.

`P371` must remain:

1. below any discharge of `N302`,
2. below provider-object realization,
3. below any selector closure / `QW-2191` discharge.

## What P371 checks

`P371` checks only:

1. `F285` exports a discharge witness packet
   `Eta_strict_provider_object_carrier_orbit_quotient_nad12_sigma_residual_welding_discharge_witness_packet_v1`,
2. the welded dictionary alias
   `W_weld_v1 := Pi_strict_provider_object_carrier_orbit_quotient_nad12_sigma_residual_welding_dictionary_candidate_v1`
   is explicit,
3. the `T128` acceptance tests are satisfied **in the declared interface
   sense**:
   - explicit weld data as data,
   - explicit 12-slot identification,
   - explicit `E_pair` carrier-field builder compatible with `T116` on
     `Nad12Index_v1 = {0..11}` (no physical canonicalization implied),
   - residual-frontier record-class compatibility (via `T127/T118` alignment),
   - noncyclic + observer-free preserved,
   - selector neutrality preserved.

## Probe table

| Check | Verdict | Evidence |
|---|---|---|
| weld target exists | YES | `T128/F281/N393` |
| weld dictionary candidate exists | YES | `T129/F282/N394` |
| discharge witness packet exists | YES | `F285` |
| weld target discharged (T128 satisfied) | YES | `F285` + `T127` + `T129` |
| `N302` discharged | NO | boundary remains |

## Exact verdict

The strongest honest current verdict is:

```text
welding target (T128) : discharged at declared interface level
N302                  : still in force
ToE closure           : not proved
```

