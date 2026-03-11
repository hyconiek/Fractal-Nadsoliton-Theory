# N420 Current First Actual Strict Alpha-Geo Equipartition Witness (16 Microstates) Strict-Derivation/Source-Upgrade Theorem

Status: `N420_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_ALPHA_GEO_EQUIPARTITION_WITNESS_16_MICROSTATE_STRICT_DERIVATION_SOURCE_UPGRADE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`T144/P384/F299/N411` name one strict-side missing ingredient:

```text
an observer-free equipartition witness of exactly 16 microstates upgrading
alpha_geo from premise/canonical status to strict-derived source status.
```

This theorem packages the strongest honest current statement that the repo now
exports such a strict witness via `F309`, with explicit scope limits.

## Theorem-level conclusion

From `F309`, the current repo exports:

1. a strict microstate object `Omega_16_v1` with `|Omega_16_v1|=16`,
2. an observer-free symmetry action `G_bit_v1 ⟲ Omega_16_v1` that is
   transitive,
3. the unique invariant equipartition measure `mu_eq_v1` on `Omega_16_v1`,
4. a strict four-bit structure witness `Omega_16_v1 ≅ {0,1}^4`,
5. a strict Shannon computation witness:

```text
alpha_geo_strict_derived_v1 := H(mu_eq_v1) = ln(16) = 4 ln 2.
```

Therefore the repo now discharges the acceptance tests of `T144` as an **actual**
strict derivation/source-upgrade witness package.

## What N420 proves

`N420` proves only:

1. existence of an exported strict equipartition witness package of size `16`,
2. strict-derived source-upgrade of `alpha_geo` via Shannon entropy on that
   exported equipartition measure,
3. observer-free and noncyclic contract satisfaction (no `K_obs`, no `theta`,
   no populated basis-pair inputs).

## What N420 does not prove

`N420` does not prove:

1. discharge of `T147/N414` (selector-track identification),
2. discharge of `QW-2191` or strict-core selector closure,
3. export of any residual-datum bridge/export-map object (`N301`) or discharge
   of `N300`,
4. ToE closure.

