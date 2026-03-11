# T130 Current Strict Sigma-Int Residual Bridge/Export-Map Actual Object Support Target Spec

Status: `T130_CURRENT_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`N302` freezes one current-state boundary on the residual-datum /
`sigma_int_candidate` third-provider route:

```text
projection + witness layers exist,
but actual bridge/export-map object support is still absent.
```

As of the current repo state, the sigma-int residual branch exports:

1. an object-to-map support projection layer into the object-support frontier
   (`N385`),
2. an object-support witness candidate and witness (`N386`, `N387`),

while remaining explicitly below:

1. actual bridge/export-map object support,
2. strict-core theta export / pair population (`N18` loop not broken),
3. admissible `S_sel_int` / selector closure.

On the updated repo state, the strict-core bridge/export-map object is already
exported (`F311/N422`), so `N300` is superseded and `N301` is discharged.

The next honest question is therefore narrower than closure:

```text
can the repo at least sharply name the missing "actual object support" layer
as one explicit future-only target object with explicit acceptance tests?
```

`T130` does **not** claim actual object support exists.

`T130` does something narrower and audit-safe:

```text
name the missing actual object support layer as one explicit future-only target
object with explicit acceptance tests
```

so that the post-`N387` frontier is no longer “missing something above witness”
but a typed, checkable target.

## Target object

If the repo cannot yet export actual object support but can name it sharply,
export one future-only target object:

```text
Lambda_residual_datum_sigma_int_bridge_export_map_object_support_target_v1
```

with the intended meaning:

```text
one explicit future-only target object for the next missing layer
above the current object-support witness (N387) and above the exported strict-core
export-map object (F311/N422)
on the sigma-int residual third-provider route
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Lambda_residual_datum_sigma_int_bridge_export_map_object_support_target_v1`
must at minimum provide:

1. **Typed object-support carrier:** one explicit object `O_support_v1`
   belonging to the residual bridge/export-map object-support lane and sitting
   above the witness layer (`N387`).
2. **Noncyclic contract:** no `theta_{1,2}` inputs and no populated basis-pair
   instance as input (respects `N18`).
3. **Observer-free contract:** no `K_obs`-indexed selection as primary source.
4. **Sigma-int discipline:** if the construction uses `sigma_int_candidate` as
   a strict-core source datum, it must keep the missing prerequisites explicit:
   - gauge-quotient safety (`N388`) must be discharged or explicitly kept open
     (no silent gauge fixing),
   - strict derivation/source upgrade (`N389`) must be discharged or explicitly
     kept open (no silent candidate→source upgrade).
5. **Map-object compatibility:** no claim that the bridge/export map is
   absent; the strict-core bridge/export-map object is already exported
   (`F311/N422`), and this target concerns only the missing object-support
   layer above that map object.
6. **Selector neutrality:** no implied admissible `S_sel_int`, no implied
   selector closure, and no implied `QW-2191` discharge.

## Hard limits

`T130` must not claim:

1. actual object support is already present,
2. discharge of `Lambda_residual_datum_sigma_int_bridge_export_map_object_support_target_v1`,
3. actual theta export / pair population,
4. admissible `S_sel_int`,
5. selector closure,
6. `QW-2191` discharge,
7. ToE closure.
