# T158 Current Strict Sigma-Int Positive-Window Delta_d Step Strict-Provenance Target Spec

Status: `T158_CURRENT_STRICT_SIGMA_INT_POSITIVE_WINDOW_DELTA_D_STEP_STRICT_PROVENANCE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After the positive-window corridor spec (`T119`) and the strict-input instantiations
(`F314/N425`, `F325/N436`), the strict sigma-int → theta candidate lane uses a
typed corridor step parameter:

```text
delta_d ∈ (0, delta_max]   where   delta_max := d_local/11.
```

After the delta_d sensitivity audit (`P403/N437`), the current repo state also
exports one explicit no-false-pass blocker fact:

```text
the computed theta-pair depends on the admissible corridor step delta_d;
therefore delta_d is a real selector slot and cannot be silently treated as “derived”.
```

The repo already accepts one **extension-scope** convention fixing this slot:

```text
AX17 (strict_extension_only): delta_d := delta_max.
```

Before `F328/N440`, strict core still lacked a dedicated delta_d value object
with explicit strict provenance/selection method.

`T158` makes that missing ingredient sharp, to prevent false pass by wording
such as “delta_d is derived from the kernel tuple” without an explicit
selection premise or an internal strict derivation.

Update (current repo state):

On the current repo state (`F328/N440`), the strict sigma-int lane now exports
one dedicated delta_d value object with explicit strict provenance
(strict-source-upgraded by explicit premise):

```text
delta_d_sigma_int_positive_window_step_strict_provenance_v1 := delta_max.
```

Therefore `T158` is no longer a “current missing-object naming” spec.
It is kept as:

1. a sharp historical target-spec for the delta_d provenance acceptance tests, and
2. an admissible target-name record (now superseded by the actual delta_d export).

## Scope

`T158` is scoped only to the **delta_d step selection** used by the strict
positive-window sigma-int → theta candidate constructions.

It does not decide:

1. strict-core `theta_1`, `theta_2` export,
2. discharge of object-support above the exported map object (`N395/T130`),
3. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
4. ToE closure.

## Target object

If the current repo cannot yet export a dedicated strict-provenance delta_d step
value object, export one future-only target object:

```text
Delta_sigma_int_positive_window_delta_d_step_strict_provenance_target_v1
```

with intended meaning:

```text
export one dedicated, observer-free, noncyclic delta_d step value object
with explicit provenance/selection method,
so downstream strict sigma-int → theta candidate work cannot silently smuggle
in a delta_d choice as if it were strict-derived.
```

This target name is now superseded by the actual exported delta_d value object:

```text
delta_d_sigma_int_positive_window_step_strict_provenance_v1
```

exported by `F328/N440`.

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Delta_sigma_int_positive_window_delta_d_step_strict_provenance_target_v1`
must at minimum provide:

1. **Dedicated delta_d value object:** an exported object/value
   `delta_d_sigma_int_positive_window_step_strict_provenance_v1` with contract:
   - `0 < delta_d_sigma_int_positive_window_step_strict_provenance_v1 <= delta_max`,
   - where `delta_max := d_local/11` is the corridor bound from `T119`.
2. **Explicit provenance classification:** the delta_d value object must be
   explicitly marked as one of:
   - `strict_derived` (with an explicit strict derivation chain on a declared strict domain), or
   - `strict_source_upgraded` (with an explicit strict-side premise selecting a value),
   and it must not be left as an implicit convention inside an instantiation artifact.
3. **Noncyclic contract:** no `theta_{1,2}` inputs and no populated basis-pair
   instance may be used as inputs (respects `N18`).
4. **Observer-free contract:** no `K_obs`-indexed selection may serve as the
   primary source of the delta_d choice.
5. **Selector discipline:** the discharge must explicitly state whether it:
   - keeps `QW-2191` open (default), or
   - adds a separated selector/symmetry-breaking premise in a non-strict scope.
In either case, no implied strict-core selector closure is permitted.

## Relation to existing extension convention

`AX17` already accepts `delta_d := delta_max` in the separated
`strict_extension_only` scope.

That acceptance is an explicit convention (premise) and does **not** by itself
count as a strict-core discharge of this target unless a dedicated delta_d value
object is exported with the required provenance classification.

## Hard limits

`T158` must not claim:

1. that delta_d selection is already strict-derived in strict core,
2. that `AX17` constitutes a strict-core discharge (it does not),
3. actual strict-core theta export / pair population,
4. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
5. ToE closure.
