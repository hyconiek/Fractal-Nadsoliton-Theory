# T142 Current Strict ToE Closure Provider-Object Carrier Contraction-Parameter Strict-Derivation/Source-Upgrade Target Spec

Status: `T142_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_CONTRACTION_PARAMETER_STRICT_DERIVATION_SOURCE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

The current strict provider-object carrier lane exports one explicit **carrier
candidate form** (`T126/F279/N391`), but it is parameterized by free strict
contraction parameters:

```text
0 < |a| < 1,
0 < |b| < 1.
```

On the current repo state, `a,b` are **not** exported as strict-derived
source-side outputs of `tau_src_candidate_v1`. They are only part of a
permitted candidate construction form.

This matters because:

1. `N327/T125/N390` require a genuinely **source-side** provider-object carrier
   layer (no extra free knobs),
2. downstream bridge-facing projections and carrier preobjects (`T127`, `T141`)
   reuse `|a|,|b|` to build finite noncyclic carriers and weights,
3. without a strict derivation/source-upgrade for `a,b`, any attempted
   promotion of the provider-object lane into an **actual** strict-core
   ingredient would risk an implicit “parameter fit” / hidden premise,
4. therefore the next honest strict-only move is to name the missing
   strict-derivation/source-upgrade ingredient for the contraction parameters
   sharply, as one future-only target with explicit acceptance tests.

This is analogous in *role* (not in content) to the sigma-int missing
strict-derivation/source-upgrade target (`T124/N389`).

`T142` does **not** claim such a strict derivation exists.

`T142` only names the missing ingredient as a target object.

## Scope

`T142` is scoped only to the strict source status of the provider-object
carrier contraction parameters `(a,b)` used in the provider-object carrier
candidate lane (`T126/F279`).

It does not decide:

1. discharge of `Epsilon_strict_provider_object_carrier_layer_target_v1`
   (`T125/N390`),
2. discharge of residual bridge/export-map object support (`N302`),
3. export of any bridge/export map object (`N300` remains in force),
4. discharge of sigma-int prerequisites (`N388/N389`),
5. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
6. ToE closure.

## Target object

If the answer is negative with respect to actual strict derivation but positive
with respect to sharp naming, export one future-only target object:

```text
Delta_strict_provider_object_carrier_contraction_parameter_strict_derivation_source_upgrade_target_v1
```

with the intended meaning:

```text
the current repo now names one exact future-only target object for the missing
strict derivation / source-object upgrade required before the provider-object
carrier contraction parameters (a,b) may be treated as strict-derived source-side
outputs rather than free candidate parameters.
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of this target must at minimum provide:

1. **Explicit source-side parameter map (no extra knobs):**
   one exported map on a declared strict source domain:

   ```text
   A_strict_provider_object_contraction_parameter_source_map_v1 :
     tau_src_candidate_v1 -> (a,b)
   ```

   with:

   ```text
   0 < |a| < 1,  0 < |b| < 1.
   ```

2. **Strict derivation / source-upgrade statement:**
   an explicit theorem-level statement that the output parameters are
   strict-derived on the declared domain (not injected as external knobs).
3. **Noncyclic contract (L5/L12, N18):**
   the construction must not take as input:
   - any `theta_{1,2}`,
   - any populated basis-pair instance.
4. **Observer-free contract:**
   no `K_obs`-indexed selection is used as a primary source of uniqueness.
5. **Gauge/symmetry discipline:**
   any gauge/symmetry action used in the carrier lane must remain explicit as an
   internal symmetry action (no silent gauge fixing counted as derivation).
6. **No silent promotion:**
   the discharge must explicitly state which statuses it upgrades and which it
   does **not** upgrade (e.g. it must not silently imply `N302` discharge).

## Hard limits

`T142` must not claim:

1. strict derivation/source upgrade already discharged,
2. discharge of `Epsilon_strict_provider_object_carrier_layer_target_v1`,
3. actual residual bridge/export-map object support (`N302`),
4. any bridge/export-map object export (`N300`),
5. admissible `S_sel_int`,
6. strict-core selector closure or `QW-2191` discharge,
7. ToE closure.

