# F317 First Actual Strict Sigma-Int `E_pair` Eps-Parameter Strict-Provenance Source-Upgrade Packet

Status: `F317_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_E_PAIR_EPS_PARAMETER_STRICT_PROVENANCE_SOURCE_UPGRADE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `F315/N426` and `T150/F316/N427`, the strict sigma-int lane has already
exported:

1. a strict sigma-int source upgrade (`F307/N418`),
2. an actual strict-core sigma-int → residual target-slot export-map object
   (`F311/N422`),
3. strict-input candidate theta-pair instantiations (`F312/N423`, `F314/N425`),

but the sigma-int driven `E_pair` generator still depended on a free amplitude
parameter:

```text
eps ∈ [0,1].
```

This packet performs the next strict anti-false-pass move demanded by `T150`:

```text
export one dedicated eps value object with explicit strict provenance
(strict-source-upgraded via an explicit premise),
so downstream strict sigma-int instantiations no longer rely on a silent
parameter choice such as eps = 1/2.
```

## Inputs reused (strict-admissible)

1. `T117/F270/N382`
   - the sigma-int-driven `E_pair` generator candidate and its parameter
     contract `eps ∈ [0,1]`,
2. `F315/N426`
   - the explicit audit showing eps had no strict provenance object and was
     used only as a fixed numeric choice in instances,
3. `T150/F316/N427`
   - the strict eps-provenance target spec, target-name record, and theorem-level
     statement keeping the missing ingredient explicit.

## Exported strict eps value object (source upgrade by explicit premise)

`F317` exports one strict eps value object:

```text
eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2.
```

Provenance classification:

- `strict_source_upgraded` (explicit strict-side premise; not strict-derived),
- observer-free (no `K_obs` primary selection),
- noncyclic (no `theta_{1,2}` inputs and no populated basis-pair inputs),
- semantically separated from unrelated `epsilon_*` radius/stability symbols
  (no bridge identifying those symbols with this generator amplitude parameter).

Persisted artifact:

```text
fundamental_action_reconstruction/generated/eps_sigma_int_E_pair_amplitude_strict_provenance_v1.json
```

## Selector discipline

This discharge keeps the global discipline unchanged:

- no implied selector closure,
- `QW-2191` remains open.

`eps` is fixed only as a strict-side source-upgrade ingredient for removing a
silent free-parameter from the sigma-int-driven generator lane. It is not a
complete selector or theta-source export.

## Status discipline

This packet does **not** claim:

1. strict derivation or uniqueness of the `E_pair` generator (`T117` remains candidate),
2. actual strict-core `theta_1`, `theta_2` export,
3. actual populated basis-pair instance,
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

