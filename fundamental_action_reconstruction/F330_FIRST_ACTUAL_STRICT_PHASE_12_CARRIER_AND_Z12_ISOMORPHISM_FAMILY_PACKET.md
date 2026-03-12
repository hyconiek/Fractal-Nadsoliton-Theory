# F330 First Actual Strict `Phase_12` Carrier + `Z_12` Isomorphism-Family Packet

Status: `F330_EXECUTED_FIRST_ACTUAL_STRICT_PHASE_12_CARRIER_AND_Z12_ISOMORPHISM_FAMILY_PACKET_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `F329/N450`, the repo exports a typed cyclic group carrier `Z_12_v1` and its regular action on the
12-slot scaffold.

`T163/P417/N451` keep explicit that the next typed ingredient for the `AX20/T162` direction is a typed
12-phase carrier and a phase embedding, *without* silently smuggling a generator/orientation choice as an
untracked selector slot.

`F330` executes the narrowest honest move **below canonicity**:

```text
export a typed 12-phase carrier Phase_12_v1 (12th-root-like finite phase group),
and export the full finite family of Z_12_v1 -> Phase_12_v1 group isomorphisms
(showing non-uniqueness explicitly),
without claiming any canonical choice and without claiming T163 discharge.
```

## Inputs reused (strict-admissible)

1. `F329/N450`
   - typed cyclic group object `Z_12_v1` and its regular action on the 12-slot scaffold.
2. `T163`
   - acceptance tests require a typed phase carrier and an embedding discipline with no hidden generator choice.

## Packet result

### 1) Typed phase carrier `Phase_12_v1`

Export one strict finite phase group object:

```text
Phase_12_v1 := { ζ^k | k=0..11 },  with  ζ^a · ζ^b := ζ^(a+b mod 12).
```

This is a “12th-root-like” phase carrier: it is a cyclic group of order 12 written in phase notation.

Persisted artifact:

`fundamental_action_reconstruction/generated/phase_12_v1_group.json`

### 2) Explicit isomorphism family `Iso(Z_12_v1, Phase_12_v1)`

Export the full finite family of group isomorphisms

```text
emb_u : Z_12_v1 -> Phase_12_v1
emb_u(k) := ζ^(u*k mod 12)
```

for units `u ∈ (Z/12Z)^× = {1,5,7,11}`.

This family has cardinality `4` and witnesses the non-uniqueness of a “phase embedding” even after both
carriers exist.

Persisted artifact:

`fundamental_action_reconstruction/generated/iso_z12_v1_to_phase_12_v1_isomorphism_family.json`

## Meaning

This packet means only:

1. the repo now exports a typed 12-phase carrier `Phase_12_v1`,
2. the repo now exports the **explicit 4-element family** of `Z_12_v1 -> Phase_12_v1` isomorphisms,
3. therefore any later claim of a *canonical* phase embedding must explicitly address the generator/orientation
   choice (or prove quotient invariance),
4. no downstream selector closure is implied (no Berry/holonomy, no theta export, no `O(2)` cut, no `QW-2191`
   discharge).

## Hard limits (no false pass)

`F330` does not claim:

1. a canonical selection of one `emb_u` (it exports the family and leaves canonicity open),
2. discharge of `T163` (phase-embedding canonicity/quotient-safety),
3. any density-operator rigidity forcing eigenvalues `1/2`,
4. any Berry/holonomy construction with gauge discipline,
5. discharge of `T162` (slot-free sigma-int → theta construction class),
6. strict theta export, strict-core selector closure, or `QW-2191` discharge,
7. ToE closure.

