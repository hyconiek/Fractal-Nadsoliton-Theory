# F307 First Actual Strict Sigma-Int FR-Sign Map + Source-Upgrade Packet

Status: `F307_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_FR_SIGN_MAP_AND_SOURCE_UPGRADE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `F306/N417`, the strict sigma-int sign route has two prerequisite
topology objects required by `T149`:

1. a strict configuration-space object `C_v1`, and
2. a strict witness `pi_1(C_v1) ≅ Z_2` with a generator label `gamma_pi1_v1`.

`T149` still demands two additional strict objects before any honest strict
sigma-int source-upgrade may be claimed under the FR-sign definition:

1. a strict FR-sign map object

```text
chi_FR_strict_v1 : pi_1(C_v1) -> {+1,-1},
```

with explicit provenance (strict-derived **or** explicit strict-side premise,
but never silent hybrid reuse), and

2. a strict sigma-int source-upgrade object:

```text
sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1).
```

This packet exports those two objects, strictly below any gauge-quotient safety
claim and strictly below any selector-track/`QW-2191` claim.

## Inputs reused (strict-admissible)

1. `F306/N417`
   - strict declared-domain object `C_v1`,
   - strict witness `pi_1(C_v1) ≅ Z_2`,
   - generator label `gamma_pi1_v1`.
2. `T149`
   - acceptance-test contract for the FR-sign map and sigma-int source upgrade.

## Construction (mathematical, no physics claim)

### Codomain group

Let:

```text
Sign := {+1,-1}
```

as a group under multiplication. Then `Sign ≅ Z_2`.

### Unique nontrivial character on pi_1(C_v1)

From `N417`, the declared strict domain satisfies:

```text
pi_1(C_v1) ≅ Z_2.
```

Therefore there exist exactly two group homomorphisms:

```text
pi_1(C_v1) -> Sign:
  (i) trivial, (ii) the unique nontrivial character.
```

### Exported map object: chi_FR_strict_v1

Export the strict FR-sign map object as the **unique nontrivial** character:

```text
chi_FR_strict_v1 : pi_1(C_v1) -> Sign
```

characterized by:

```text
chi_FR_strict_v1(e) = +1,
chi_FR_strict_v1(gamma_pi1_v1) = -1,
```

where `e` is the identity loop class and `gamma_pi1_v1` is the nontrivial loop
class label from `F306/N417`.

### Provenance declaration (strict-side premise, no hybrid reuse)

This packet exports `chi_FR_strict_v1` with explicit provenance:

```text
provenance(chi_FR_strict_v1) = explicit_strict_side_premise_nontrivial_character
```

Meaning:

1. no hybrid-only FR quantization sources (e.g. `QW-1622`) are used to justify
   the sign,
2. no theorem-level physical derivation is claimed here,
3. the map is introduced explicitly as the minimal nontrivial strict-side
   character enabled by `pi_1(C_v1) ≅ Z_2`.

### Exported strict sigma-int source upgrade

Define and export the strict sigma-int source-upgrade object:

```text
sigma_int_strict_derived_v1 := chi_FR_strict_v1(gamma_pi1_v1) = -1.
```

This upgrades sigma-int from candidate-only status to a strict-core **source**
datum on the declared domain, with the above explicit premise-based
provenance.

## Persisted artifacts

1. `fundamental_action_reconstruction/generated/chi_fr_strict_v1_map_on_pi1_c_v1_v1.json`
2. `fundamental_action_reconstruction/generated/sigma_int_strict_derived_v1.json`

## Status discipline

This packet does **not** claim:

1. theorem-level gauge-quotient safety (`T123/N388`),
2. selector-track identification beyond overlay-only (`T147/N414`),
3. export of any residual-datum bridge/export-map object (`N301`), nor
   discharge of `N300`,
4. admissible `S_sel_int`, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

