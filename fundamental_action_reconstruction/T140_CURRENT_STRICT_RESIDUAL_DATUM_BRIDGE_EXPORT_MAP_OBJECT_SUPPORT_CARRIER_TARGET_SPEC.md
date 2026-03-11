# T140 Current Strict Residual Datum Bridge/Export-Map Actual Object-Support Carrier Target Spec

Status: `T140_CURRENT_STRICT_RESIDUAL_DATUM_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_CARRIER_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `N387` and `N405`, the strict residual bridge/export-map lane now exports
two separate **witness-level** arrivals at the bridge/export-map object-support
frontier:

1. sigma-int third-provider witness (`N387`),
2. provider-object carrier residual witness (`N405`),

while still remaining explicitly below **actual bridge/export-map object
support** (`N302`).

`T130` and `T139` already name the missing **post-witness layer** as a
future-only target object (`Lambda_*_object_support_target_v1`) and use the
placeholder requirement:

```text
"export one explicit typed object-support carrier O_support_v1 above witness"
```

The next honest audit-safe move (still below closure) is to remove vagueness:

```text
name the missing post-witness object-support carrier itself
as one explicit future-only target object with explicit acceptance tests,
so that a future “object support” claim cannot be satisfied by relabeling
projection/witness packets.
```

`T140` does **not** claim any such carrier is currently exported.

## Target object

If the current repo cannot export an actual object-support carrier but can name
the missing layer sharply, export one explicit future-only target object:

```text
Omicron_residual_datum_bridge_export_map_object_support_carrier_target_v1
```

with the intended meaning:

```text
one explicit future-only target object for the next missing "actual object
support carrier" layer above the current witness layer(s)
(`N387` sigma-int witness and/or `N405` provider-object witness),
and above the exported strict-core export-map object (`F311/N422`).
```

This target is scoped to the strict residual-datum bridge/export-map lane and
is intended to be compatible with both route-local witnesses:

1. `Kappa_residual_datum_sigma_int_bridge_export_map_object_support_witness_v1`
   (`N387`),
2. `Kappa_residual_datum_provider_object_carrier_bridge_export_map_object_support_witness_v1`
   (`N405`).

On the current repo state, this target is discharged by the exported carrier
object `T146/F301/N413`.

## Acceptance tests (what would count as discharge)

An **actual** inhabitant of
`Omicron_residual_datum_bridge_export_map_object_support_carrier_target_v1`
must at minimum provide:

1. **Typed carrier object (new layer above witness):** one explicit object
   `O_support_v1` that is declared as the route-local bridge/export-map
   object-support carrier sitting strictly above:
   - the sigma-int witness layer (`N387`), and/or
   - the provider-object carrier witness layer (`N405`),
   and strictly below:
   - actual bridge/export-map object support above the map object (`N395`).
2. **Route reference discipline:** `O_support_v1` must carry, as explicit data,
   which strict witness layer it upgrades (sigma-int, provider-object, or both),
   by referencing the exact witness object name(s) above.
3. **Noncyclic contract:** `O_support_v1` must be constructed without taking as
   input:
   - any `theta_{1,2}`,
   - any populated basis-pair instance,
   (respects sandbox `N18`).
4. **Observer-free contract:** `O_support_v1` must not use `K_obs`-indexed
   selection as a primary source of uniqueness.
5. **Map-object compatibility:** the construction must not silently claim that
   the bridge/export map is absent. The strict-core export-map object is
   already exported (`F311/N422`).
6. **Sigma-int discipline (if used):** if the construction uses
   `sigma_int_candidate` as a strict-core source datum, it must keep the two
   missing prerequisites explicit:
   - gauge-quotient safety (`N388`) must be discharged or explicitly kept open
     (no silent gauge fixing),
   - strict derivation/source upgrade (`N389`) must be discharged or explicitly
     kept open (no silent candidate→source upgrade).
7. **Selector neutrality:** `O_support_v1` must not imply:
   - admissible `S_sel_int`,
   - strict-core selector closure,
   - `QW-2191` discharge,
   - ToE closure.

## Hard limits

`T140` must not claim:

1. discharge of `N302`,
2. actual bridge/export-map object support,
3. any new export-map object export beyond `F311/N422`,
4. actual theta export / pair population,
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.
