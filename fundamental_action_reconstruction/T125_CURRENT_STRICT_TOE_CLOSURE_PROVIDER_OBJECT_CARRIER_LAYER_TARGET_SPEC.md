# T125 Current Strict ToE Closure Provider-Object Carrier Layer Target Spec

Status: `T125_CURRENT_STRICT_TOE_CLOSURE_PROVIDER_OBJECT_CARRIER_LAYER_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

`N327` sharpens the dominant missing strict-closure ingredient class as one:

```text
source-side, observer-free, pair-indexed, noncyclic
strict selector/provider object-carrier layer.
```

`N370-N376` further isolate that gap as the provider-object realization-side
arm, but they do not yet export the missing carrier layer itself.

`T125` does **not** claim any realization.

`T125` does something narrower and more audit-safe:

```text
name one explicit future-only target object
for the missing provider-object carrier layer
with an explicit minimal type-shape and explicit acceptance tests
```

so that “provider object” does not remain a purely narrative placeholder.

## Context (current strict closure lane)

On the current repo state:

1. `N327` exports the dominant missing ingredient class diagnosis.
2. `N370` exports the noncyclic realization split target
   `Delta_prov / Rho_orient`.
3. `N376` exports one provider-object arm support witness, but remains strictly
   below actual provider-object realization.
4. `N302` blocks the nearest residual-datum / `sigma_int_candidate` branch
   below actual bridge/export-map object support.
5. `N385` (above `N384`) exports an actual object-to-map support projection
   layer, but still remains below actual bridge/export-map object support.
6. `N388/N389` name two explicit sigma-int prerequisites (future-only targets)
   that must not be silently assumed discharged if `sigma_int_candidate` is
   used as an internal source datum.

Therefore the next honest strict-only move is still not a closure claim.
It is to sharply name the missing provider-object carrier layer as a typed
future-only target with explicit acceptance tests.

## Target object

If the repo cannot yet export an actual provider-object carrier layer but can
still name it sharply, export:

```text
Epsilon_strict_provider_object_carrier_layer_target_v1
```

with the intended minimal type-shape:

```text
C_prov_target_v1 :
  tau_src_candidate_v1 -> ProviderObjectCarrier_pair
```

where:

1. `tau_src_candidate_v1` is the strict source-side candidate domain used
   across the strict lane (`N284/T26`),
2. `ProviderObjectCarrier_pair` is a pair-indexed carrier type (at minimum on
   the designated pair family `{pair1, pair2}`),
3. the map is **source-side** and **observer-free** (no `K_obs` as primary
   selector source),
4. the map is **noncyclic** with respect to the sandbox loop (`N18`):
   - no `theta_1,theta_2` as inputs,
   - no populated basis-pair instance as input.

## Acceptance tests (what would count as discharge)

An **actual** inhabitant of
`Epsilon_strict_provider_object_carrier_layer_target_v1` must at minimum
provide:

1. **Explicit pair indexing:** an explicit pair index set and an explicit
   indexing map carried as data (not implied by narrative).
2. **Observer-free contract:** no use of `K_obs` or observer-indexed selection
   as a source of uniqueness.
3. **Noncyclic contract:** no `theta` inputs and no populated-instance inputs
   (respects `N18`).
4. **Bridge-facing projection interface:** an explicit projection interface
   from provider-object carrier data into the residual bridge/export-map
   object-support frontier, i.e. something of the form:

   ```text
   ProviderObjectCarrier_pair
     -> ResidualBridgeExportMapSupportCarrier_pair
   ```

   where the output is typed so it can be used to attack the `N302` object
  -support gap without introducing cycles. (It must remain below actual
   bridge/export-map object support.)
5. **Sigma-int discipline (if used):**
   - the construction must not silently assume strict derivation/source-upgrade
     of `sigma_int_candidate` (see `N389`),
   - the construction must not silently assume gauge-quotient safety of
     `sigma_int_candidate` (see `N388`).
6. **Selector neutrality:** no claim of admissible `S_sel_int`, no implicit
   strict-core selector closure, and no implied `QW-2191` discharge.

## Hard limits

`T125` must not claim:

1. actual provider-object realization,
2. actual bridge/export-map object support,
3. actual `E_orient`,
4. admissible `S_sel_int`,
5. strict-core selector closure,
6. `QW-2191` discharge,
7. ToE closure.

