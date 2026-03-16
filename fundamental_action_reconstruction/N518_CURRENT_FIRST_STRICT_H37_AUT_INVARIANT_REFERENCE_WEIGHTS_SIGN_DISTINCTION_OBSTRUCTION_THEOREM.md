# N518 Current First Strict H37 Aut‑Invariant Reference Weights Sign‑Distinction Obstruction Theorem

Status: `N518_DISCHARGED_CURRENT_FIRST_STRICT_H37_AUT_INVARIANT_REFERENCE_WEIGHTS_SIGN_DISTINCTION_OBSTRUCTION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

Record the strongest honest strict obstruction, above `N517`, explaining why direction‑free (`Aut(Z_12)`‑invariant) reference weights cannot supply a strict sign‑distinction observable on `pair1` of the scalar form `Σ_x w(x) u_1(x)` on the current exported `pair1` sine axis.

## Theorem (strict obstruction)

Let `n=12` and let `Z_n` be written additively.

Assume:

1. `w : Z_n -> R` is `Aut(Z_n)`‑invariant, i.e. `w(a·x)=w(x)` for every unit `a ∈ (Z_n)^×` and every `x ∈ Z_n`;
2. `u : Z_n -> R` is odd under inversion: `u(-x) = -u(x)` for every `x ∈ Z_n`.

Then:

```text
Σ_{x∈Z_n} w(x) u(x) = 0.
```

In particular, for the current exported strict `pair1` representative `u_1 = ± s1` (odd under inversion) and any direction‑free (`Aut(Z_12)`‑invariant) reference weight family `w`, the scalar observable candidate `Σ_x w(x) u_1(x)` cannot distinguish `u_1` from `-u_1` (it vanishes).

## Proof (one line)

Since `-1 ∈ (Z_n)^×`, inversion `x↦-x` is an automorphism, hence Aut‑invariance implies `w(-x)=w(x)`.
Therefore:

```text
Σ_x w(x)u(x) = 1/2 Σ_x ( w(x)u(x) + w(-x)u(-x) )
            = 1/2 Σ_x ( w(x)u(x) + w(x)(-u(x)) )
            = 0.
```

## Consequence for `H37`

Any strict attempt to extract a sign‑distinction observable on the current `pair1` sine axis via a direction‑free (`Aut(Z_12)`‑invariant) scalar weight family and a linear scalar of the form `Σ_x w(x) u_1(x)` is blocked.

Therefore, discharging `H37` in strict core requires at least one of:

1. an explicit reflection‑breaking / orientation source (a new strict internal datum not invariant under inversion),
2. a different observable class not reducible to direction‑free linear scalar weights on `Z_12`.

This theorem does **not** claim either (1) or (2) exists on current repo state.

## Realization artifact (optional)

`F472` exports an explicit obstruction object on the current exported strict instance, enumerating `Aut(Z_12)` orbits and verifying cancellation for the orbit‑basis weight family.

## What N518 does not prove

`N518` does not prove:

1. discharge of `H37` (it proves an obstruction, not a sign-sensitive physical datum),
2. that no sign‑distinction observable exists in strict core in every possible formalism,
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.

