# F309 First Actual Strict Alpha-Geo Equipartition Witness (16 Microstates) Strict-Derivation/Source-Upgrade Packet

Status: `F309_EXECUTED_FIRST_ACTUAL_STRICT_ALPHA_GEO_EQUIPARTITION_WITNESS_16_MICROSTATE_STRICT_DERIVATION_SOURCE_UPGRADE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`T144/P384/F299/N411` keep one strict-side missing ingredient explicit:

```text
alpha_geo := 4 ln 2 is usable as a strategic constant,
but no strict microstate/equipartition witness is exported that upgrades it
into a strict-derived source status.
```

This packet executes the narrow discharge construction demanded by `T144`:

```text
export Omega_16_v1, export mu_eq_v1 forced by an observer-free symmetry action,
export a four-bit structure witness, and export
alpha_geo_strict_derived_v1 := H(mu_eq_v1) = ln(16) = 4 ln 2,
without importing any legacy-only operator decompositions as strict sources.
```

## Inputs reused (strict-admissible)

1. `T144`
   - acceptance tests for the strict equipartition witness discharge.
2. `S2`
   - strict-only closure priority; `4 ln 2` is permitted as strict-side premise
     but forbidden to silently promote without a strict witness.

## Strict microstate object (|Omega| = 16)

Export a strict microstate object:

```text
Omega_16_v1 := {0,1}^4
```

so `|Omega_16_v1| = 16` by construction.

Persisted artifact:

`fundamental_action_reconstruction/generated/omega_16_v1_microstate_space.json`

## Observer-free symmetry action forcing equipartition

Let:

```text
G_bit_v1 := (Z_2)^4
```

act on `Omega_16_v1` by bitwise addition (XOR):

```text
g ⋅ ω := ω ⊕ g.
```

This action is transitive: for any `ω,ω' ∈ Omega_16_v1`, choose `g = ω ⊕ ω'`
so that `g⋅ω = ω'`.

Therefore any `G_bit_v1`-invariant probability measure on `Omega_16_v1` must be
uniform (equipartition), because invariance and transitivity force:

```text
mu({ω}) = mu({ω'})  for all ω, ω',
```

and the probabilities must sum to `1`.

## Equipartition measure

Export the strict equipartition measure `mu_eq_v1` as the unique
`G_bit_v1`-invariant probability measure on `Omega_16_v1`:

```text
mu_eq_v1({ω}) = 1/16   for all ω ∈ Omega_16_v1.
```

Persisted artifact:

`fundamental_action_reconstruction/generated/mu_eq_v1_equipartition_measure_on_omega_16_v1.json`

## Four-bit structure witness

Export an explicit four-bit structure witness:

```text
Omega_16_v1 ≅ {0,1}^4
```

with the canonical identification (the identity map in this model).

This witnesses four independent binary degrees of freedom and forces
`2^4 = 16` microstates.

Persisted artifact:

`fundamental_action_reconstruction/generated/omega_16_v1_four_bit_structure_witness.json`

## Shannon computation (strict alpha-geo source upgrade)

Compute the Shannon entropy (natural log):

```text
H(mu_eq_v1)
  := - Σ_{ω ∈ Omega_16_v1} mu_eq_v1({ω}) ln(mu_eq_v1({ω}))
   = - 16 * (1/16) * ln(1/16)
   = ln(16)
   = 4 ln(2).
```

Export the strict-derived status upgrade object:

```text
alpha_geo_strict_derived_v1 := H(mu_eq_v1) = ln(16) = 4 ln 2.
```

Persisted artifact:

`fundamental_action_reconstruction/generated/alpha_geo_strict_derived_v1.json`

## Status discipline

This packet does **not** claim:

1. discharge of `T147/N414` (selector-track identification),
2. discharge of `QW-2191` or strict-core selector closure,
3. export of any residual-datum bridge/export-map object (`N301`) or discharge
   of `N300`,
4. ToE closure.

