# F333 First Actual Strict `Phase_12` / `Aut(Z_12)` Quotient Orbit-Decomposition Packet

Status: `F333_EXECUTED_FIRST_ACTUAL_STRICT_PHASE_12_AUT_Z12_QUOTIENT_ORBIT_DECOMPOSITION_PACKET_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After:

1. `F330/N452` (typed phase carrier `Phase_12_v1` and the explicit 4-element embedding family
   `Iso(Z_12_v1, Phase_12_v1)`),
2. `F331/N453` (typed gauge symmetry `Aut_Z12_v1 ≅ (Z/12Z)^×` acting on the phase layer),

the remaining honest canonicity discipline (`T163/P417/P418/P419`) is:

```text
either prove downstream invariance under Aut(Z_12),
or explicitly introduce symmetry breaking (non-strict).
```

This packet executes the narrowest strictly mathematical move that supports the **quotient-safe** option:

```text
export the explicit orbit decomposition of Phase_12_v1 under the Aut_Z12_v1 action,
and export the quotient-set object Phase_12_v1 / Aut_Z12_v1 as a typed strict carrier.
```

This does **not** claim any theta export, any Berry/holonomy ingredient, any `O(2)` cut, nor any `QW-2191`
discharge.

## Inputs reused (strict-admissible)

1. `F330/N452`
   - typed `Phase_12_v1` carrier.
2. `F331/N453`
   - typed `Aut_Z12_v1` group and its action on `Phase_12_v1`:
     `alpha_u(ζ^k) := ζ^(u*k mod 12)`.
3. `F332/N454`
   - parity character is Aut-invariant and therefore descends to the quotient.

## Packet result (quotient object + explicit orbits)

### 1) Explicit orbit decomposition

Under the action `alpha_u(ζ^k) := ζ^(u*k mod 12)`, the orbit decomposition has exactly 6 orbits:

1. `{ζ^0}` (fixed),
2. `{ζ^6}` (fixed),
3. `{ζ^1, ζ^5, ζ^7, ζ^11}` (the unit-exponent orbit),
4. `{ζ^2, ζ^10}`,
5. `{ζ^3, ζ^9}`,
6. `{ζ^4, ζ^8}`.

This follows from the standard classification by `gcd(k,12)` under multiplication by units.

### 2) Typed quotient-set object

Export the quotient carrier:

```text
Phase_12_v1_aut_z12_quotient_v1 := Phase_12_v1 / Aut_Z12_v1
```

materialized by:

1. explicit orbit IDs,
2. explicit orbit members,
3. the quotient map `q : Phase_12_v1 -> Phase_12_v1/Aut_Z12_v1`.

Persisted artifact:

`fundamental_action_reconstruction/generated/phase_12_v1_aut_z12_quotient_orbits.json`

### 3) Parity character descends to the quotient (quotient-safe invariant)

Because the parity character `chi_parity_Phase_12_v1(ζ^k)=(-1)^(k mod 2)` is Aut-invariant (`N454`),
it factors through the quotient:

```text
chi_parity_Phase_12_v1 = chi_parity_quotient ∘ q
```

so any downstream construction that depends only on `q(ζ^k)` (or on invariants that factor through the
quotient) is automatically free of the “generator/orientation choice” at this layer.

## Meaning

This packet means only:

1. the phase-embedding ambiguity is now represented as an explicit finite gauge action with an explicit finite
   quotient carrier,
2. Aut-invariant phase-layer ingredients can be stated as functions on the quotient
   `Phase_12_v1/Aut_Z12_v1` (quotient-safe discipline),
3. this does not yet provide any strict theta export, density/holonomy ingredient, or `QW-2191` selector
   closure.

## Hard limits (no false pass)

`F333` does not claim:

1. discharge of `T163` or `T162`,
2. any strict density operator forcing `1/2`,
3. any Berry/holonomy construction,
4. strict theta export, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

