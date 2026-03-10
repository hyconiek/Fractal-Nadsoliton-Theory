# T144 Current Strict Alpha-Geo Equipartition Witness (16 Microstates) Strict-Derivation/Source-Upgrade Target Spec

Status: `T144_CURRENT_STRICT_ALPHA_GEO_EQUIPARTITION_WITNESS_16_MICROSTATE_STRICT_DERIVATION_SOURCE_UPGRADE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-10`

## Goal

On the current repo state:

1. `F1` exports `alpha_geo := 4 ln 2` only as a
   `canonical-ontology-supported` parameter identity (not strict-derived).
2. `S2` treats `4 ln 2` as a strict-side strategic premise for candidate
   construction, while explicitly forbidding any silent promotion into actual
   discharge.

This creates one named strict-side tension:

```text
alpha_geo = 4 ln 2 is usable as a strategic constant,
but the repo still lacks any exported strict-core derivation/source object
that makes "16 equally weighted microstates" operational as an actual source
ingredient.
```

`T144` names that missing ingredient sharply as one future-only target object
with explicit acceptance tests, following the strict status vocabulary.

## Scope

`T144` is scoped only to the strict derivation / source-upgrade of:

```text
alpha_geo = 4 ln 2
```

It does not decide:

1. strict derivation/source-upgrade of `sigma_int_candidate` (`T124`),
2. theorem-level gauge-quotient safety for `sigma_int_candidate` (`T123`),
3. discharge of `N302` (actual residual bridge/export-map object support),
4. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
5. ToE closure.

## Target object

If the repo cannot yet export an actual strict equipartition witness but can
still name it sharply, export one future-only target object:

```text
Delta_alpha_geo_equipartition_witness_16_microstate_strict_derivation_source_upgrade_target_v1
```

with the intended meaning:

```text
export an observer-free, source-side, noncyclic strict object/packet that
forces an actual equipartition witness of exactly 16 microstates and upgrades
alpha_geo from canonical/premise status to strict-derived source status.
```

## Acceptance tests (what would count as discharge)

An **actual** discharge of
`Delta_alpha_geo_equipartition_witness_16_microstate_strict_derivation_source_upgrade_target_v1`
must at minimum provide:

1. **Strict microstate object:** an exported strict object `Omega_16_v1` with
   `|Omega_16_v1| = 16`.
2. **Equipartition measure:** an exported strict probability measure
   `mu_eq_v1` on `Omega_16_v1` with:
   - `mu_eq_v1({ω}) = 1/16` for every `ω ∈ Omega_16_v1`, and
   - an explicit proof obligation showing why equipartition is forced (not a
     free choice).
3. **Four-bit structure witness:** an exported strict witness that the
   microstate space carries exactly four independent binary degrees of freedom,
   e.g. an isomorphism:

   ```text
   Omega_16_v1 ≅ {0,1}^4
   ```

   or an equivalent invariant structure that forces `2^4` microstates.
4. **Observer-free / gauge discipline:** an explicit observer-free invariance
   contract (no `K_obs` as a primary selector source), ideally expressed as a
   gauge/symmetry group action `G ⟲ Omega_16_v1` with:
   - `mu_eq_v1` gauge-invariant, and
   - the four-bit witness gauge-invariant.
5. **Noncyclic contract:** no use of `theta_{1,2}` and no populated basis-pair
   instance as input (respects `N18`).
6. **Explicit alpha-geo derivation:** an exported strict derivation/witness
   that computes the Shannon entropy of `mu_eq_v1`:

   ```text
   H(mu_eq_v1) = ln(16) = 4 ln(2)
   ```

   and upgrades the repo status of `alpha_geo` by exporting an explicit object
   of the form:

   ```text
   alpha_geo_strict_derived_v1 := H(mu_eq_v1)
   ```

7. **Kernel-split and legacy discipline:** the construction must not rely on
   legacy-only operator decompositions (e.g. `K_geo/K_res/K_tors/K_topo`)
   unless those are first exported as strict-side objects with an explicit
   role-transfer theorem.

## Hard limits

`T144` must not claim:

1. that the target is already discharged,
2. that `alpha_geo` is already strict-derived on the current repo state,
3. discharge of `T124` (sigma-int derivation) or `T123` (gauge quotient
   safety),
4. discharge of `N302` or any bridge/export-map object export,
5. admissible `S_sel_int`, selector closure, or `QW-2191` discharge,
6. ToE closure.

