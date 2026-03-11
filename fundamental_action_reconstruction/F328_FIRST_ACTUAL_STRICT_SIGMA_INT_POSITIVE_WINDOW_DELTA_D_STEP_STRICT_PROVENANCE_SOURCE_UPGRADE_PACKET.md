# F328 First Actual Strict Sigma-Int Positive-Window Delta_d Step Strict-Provenance Source-Upgrade Packet

Status: `F328_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_POSITIVE_WINDOW_DELTA_D_STEP_STRICT_PROVENANCE_SOURCE_UPGRADE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After the positive-window corridor spec (`T119`) and the delta_d sensitivity
audit (`P403/N437`), the strict sigma-int → theta candidate lane contains a real
selector slot:

```text
delta_d ∈ (0, delta_max]   where   delta_max := d_local/11.
```

After `T158/F327/N439`, that missing ingredient is now sharply named as one
future-only target object with explicit acceptance tests.

`F328` performs the next honest move demanded by those acceptance tests:

```text
export one dedicated delta_d value object with explicit strict provenance
(strict-source-upgraded via an explicit premise),
so downstream strict sigma-int → theta candidate work cannot silently smuggle
in a delta_d choice as if it were strict-derived.
```

## Inputs reused (strict-admissible)

1. `T119`
   - corridor definitions `delta_barrier`, `d_local`, `delta_max := d_local/11`,
2. `P403/N437`
   - delta_d nonuniqueness/sensitivity is explicit (no false pass),
3. `T158`
   - strict-provenance delta_d target spec and acceptance tests.

## Exported strict delta_d value object (source upgrade by explicit premise)

`F328` exports one dedicated delta_d step value object:

```text
delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max.
```

Provenance classification:

- `strict_source_upgraded` (explicit strict-side premise: corridor saturation),
- observer-free (no `K_obs` primary selection),
- noncyclic (no `theta_{1,2}` inputs and no populated basis-pair inputs),
- scope-limited to the `T119` positive-window sigma-int lane.

Persisted artifact:

```text
fundamental_action_reconstruction/generated/delta_d_sigma_int_positive_window_step_strict_provenance_v1.json
```

## Implementation

Executed by:

```text
python3 fundamental_action_reconstruction/f328_first_actual_strict_sigma_int_positive_window_delta_d_step_strict_provenance_source_upgrade_packet.py
```

Outputs:

```text
fundamental_action_reconstruction/generated/delta_d_sigma_int_positive_window_step_strict_provenance_v1.json
fundamental_action_reconstruction/generated/f328_first_actual_strict_sigma_int_positive_window_delta_d_step_strict_provenance_source_upgrade_packet_summary.json
```

## Status discipline

This packet does **not** claim:

1. strict derivation or uniqueness of delta_d,
2. actual strict-core `theta_1`, `theta_2` export,
3. discharge of object-support above the exported map object (`N395/T130`),
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

