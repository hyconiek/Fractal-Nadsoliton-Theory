# N454 Current First Actual Strict `Z_12` Parity Character + `Aut(Z_12)` Invariance Witness Theorem

Status: `N454_DISCHARGED_CURRENT_FIRST_ACTUAL_STRICT_Z12_PARITY_CHARACTER_AND_AUT_INVARIANCE_WITNESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Package theorem-level the narrowest honest statement enabled by `F332`, without false pass:

```text
the Z_12 parity character χ_parity(k)=(-1)^(k mod 2) is a typed strict object,
and it is invariant under Aut(Z_12_v1).
```

This theorem is intentionally scoped below any theta export and below any `QW-2191` claim.

## Theorem-level conclusion

From `F329/N450` the repo exports a typed cyclic group object `Z_12_v1`.

From `F331/N453` the repo exports:

```text
Aut_Z12_v1 ≅ (Z/12Z)^× = {1,5,7,11}
```

acting on `Z_12_v1` by multiplication:

```text
u · k := (u*k mod 12).
```

From `F332`, the repo exports the parity character:

```text
chi_parity_Z12_v1(k) := (-1)^(k mod 2).
```

Then:

```text
for all u ∈ Aut_Z12_v1 and all k ∈ Z_12_v1:
  chi_parity_Z12_v1(u · k) = chi_parity_Z12_v1(k).
```

Proof (finite arithmetic):

1. every `u ∈ Aut_Z12_v1` is odd (`u ∈ {1,5,7,11}`),
2. therefore `(u*k) mod 2 = k mod 2`,
3. hence `(-1)^(u*k mod 2) = (-1)^(k mod 2)`.

So the parity character is `Aut(Z_12)`-invariant.

## Consequence (strict lane hygiene)

This theorem supports one precise “no hidden slot” statement:

```text
any dependence on the parity class k mod 2 (equivalently on χ_parity)
does not depend on the Aut(Z_12) generator/orientation gauge.
```

This does **not** supply:

- a strict `O(2)` cut source against `QW-2191`,
- a density-operator rigidity forcing `1/2`,
- a Berry/holonomy construction,
- a slot-free theta export (`T162`).

## What N454 does not prove

`N454` does not prove:

1. discharge of `T163` (canonical embedding) beyond the stated invariance,
2. strict-core theta export,
3. strict-core selector closure or `QW-2191` discharge,
4. ToE closure.

