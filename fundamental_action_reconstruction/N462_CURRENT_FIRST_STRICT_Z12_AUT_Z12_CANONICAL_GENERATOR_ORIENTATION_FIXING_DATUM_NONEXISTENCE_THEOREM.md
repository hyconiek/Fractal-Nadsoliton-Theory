# N462 Current First Strict `Z_12` / `Aut(Z_12)` Canonical Generator/Orientation Fixing‑Datum Nonexistence Theorem

Status: `N462_DISCHARGED_CURRENT_FIRST_STRICT_Z12_AUT_Z12_CANONICAL_GENERATOR_ORIENTATION_FIXING_DATUM_NONEXISTENCE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `F329..F333`, the typed `AX20/T162` lane contains:

1. a typed cyclic group carrier `Z_12_v1` (`F329/N450`),
2. a typed phase carrier `Phase_12_v1` and a non-unique isomorphism family
   `emb_u : Z_12_v1 -> Phase_12_v1` (`F330/N452`),
3. the typed gauge symmetry `Aut_Z12_v1` acting on both carriers and on the embedding family (`F331/N453`),
4. explicit orbit/quotient infrastructure (`F333/N455`).

`T163` keeps explicit that any “slotless phase embedding” must eliminate the hidden generator/orientation
choice, either:

1. by an exported canonical-fixing datum (`T163(4a)`), or
2. by quotient invariance of downstream numeric outputs (`T163(4b)`).

The boundary theorems `N456` and `N461` already show that demanding *pure* Aut-invariance kills the common
“go once around the 12-cycle” constructions and collapses phase information to the parity-only sector.

`N462` closes one additional strict canonicity point that is frequently left implicit:

```text
there is no Aut(Z_12)-invariant way to canonically pick a generator / orientation of Z_12_v1 (or Phase_12_v1)
from the typed Z_12 / Aut(Z_12) structure alone.
```

This does **not** deny that a future symmetry-breaking datum could exist; it only proves that such a datum
cannot be obtained “for free” while treating `Aut(Z_12)` as a gauge freedom.

## Setup (exported objects)

From `F329/N450`, the repo exports the typed group:

```text
Z_12_v1 := (I_12_v1, + mod 12).
```

From `F331/N453`, the repo exports the typed automorphism group:

```text
Aut_Z12_v1 = {1,5,7,11},
alpha_u(k) := u*k mod 12.
```

The set of generators of `Z_12_v1` (elements of order `12`) is:

```text
Gen(Z_12_v1) = {1,5,7,11}.
```

## Statement

Call a “canonical generator/orientation fixing datum” **Aut-invariant** if it would, using only the typed
`Z_12_v1` object (and its exported `Aut_Z12_v1` action), canonically select a generator (or equivalently an
oriented successor step) in a way that is invariant under the full gauge symmetry.

Concretely, suppose there existed an Aut-invariant selection map:

```text
sel : (Z_12_v1 as a typed object) -> Gen(Z_12_v1)
```

such that the selected element is invariant under all automorphisms:

```text
alpha_u(sel(Z_12_v1)) = sel(Z_12_v1)    for all u ∈ Aut_Z12_v1.
```

Then no such `sel` exists.

Equivalently:

> There exists **no** Aut(`Z_12`)-invariant canonical way to select a generator / orientation of `Z_12_v1`
> from the typed `Z_12_v1` + `Aut_Z12_v1` structure alone.

## Proof (finite, strict)

Assume, for contradiction, that such an Aut-invariant `sel` exists.

Let:

```text
g := sel(Z_12_v1) ∈ Gen(Z_12_v1) = {1,5,7,11}.
```

Aut-invariance gives:

```text
alpha_u(g) = g    for all u ∈ Aut_Z12_v1.
```

But the `Aut_Z12_v1` action on `Gen(Z_12_v1)` is nontrivial. In particular:

```text
alpha_5(1) = 5,    alpha_5(5) = 1,    alpha_5(7) = 11,    alpha_5(11) = 7.
```

So no element of `{1,5,7,11}` is fixed by all of `Aut_Z12_v1`.

This contradicts the required Aut-invariance of `g`.

Therefore no such Aut-invariant canonical generator/orientation selection exists.

## Consequence (for T163/T164 discipline)

This theorem implies:

1. any attempt to treat “choose a generator/orientation of `Z_12`” as *canonical* while maintaining the full
   `Aut(Z_12)` gauge discipline is invalid in strict core,
2. therefore, if a future construction needs a nontrivial oriented 12-cycle (or a nontrivial phase embedding
   with a fixed generator), it must either:
   - introduce an additional strict datum that explicitly breaks/restricts the `Aut(Z_12)` ambiguity
     (a symmetry-breaking fixing datum; `T164`), or
   - remain quotient-safe and accept that the Aut-invariant sector collapses to parity-only information
     (`N461`), which cannot supply nontrivial theta angles.

## What N462 does not prove

`N462` does not prove:

1. impossibility in principle of adding a symmetry-breaking fixing datum (`T164`) in a declared scope,
2. discharge of `T163` or `T162`,
3. any Berry/holonomy ingredient,
4. any strict theta export or `QW-2191` discharge,
5. ToE closure.

