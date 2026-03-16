# N517 Current First Strict `H37` Even Reference Weights Sign‑Distinction Obstruction Theorem (No False‑PASS)

Status: `N517_DISCHARGED_CURRENT_FIRST_STRICT_H37_EVEN_REFERENCE_WEIGHTS_SIGN_DISTINCTION_OBSTRUCTION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

`H37` asks whether strict core exports any sign‑sensitive state object or observable distinguishing `u` from `-u` on `pair1`.

One tempting strict route is:

```text
use the strict internal ord-reference shape (ord_Z12 / r_ord) to define a sign-sensitive scalar S(u)=Σ_x w(x)u(x)
and use its sign to distinguish u from -u.
```

This theorem closes that route on the **current exported strict instance** where the `pair1` axis is exported as `u_1 = ± s1`:
for any **even** weight function `w(x)=w(-x)`, the weighted sum against the odd mode `s1` cancels exactly.
In particular, both `ord_Z12(x)` and `r_ord(x) ∝ exp(-alpha_geo*ord_Z12(x))` are even under reflection, so they cannot supply
a sign-distinction observable of this form on `pair1` for this instance.

This does **not** prove `H37` is impossible in strict core. It only blocks one specific recurrent attempt.

## Strict‑admissible inputs reused

1. `F446`
   - exports `ord_Z12(x)` values in `r_ord_z12_v1_reference_distribution`.
2. `F456`
   - exports the current strict representative `u_1` and its `(c1,s1)` coordinates (projector-level sign-gauge discipline remains explicit).
3. `F471`
   - exports an explicit obstruction object recording the cancellation on the current exported instance (as an artifact).
4. `A10`
   - anti‑overclaim boundary.

## Theorem (even weights cannot distinguish sign on the exported `pair1` sine axis)

Let `n=12` and define the reflection map on the site index:

```text
ρ(x) := (-x) mod 12.
```

Let `w : Z_12 -> R` be any weight function satisfying:

```text
w(ρ(x)) = w(x)  (even under reflection).
```

Let `s1` be the canonical real Fourier sine mode on `Z_12`, so:

```text
s1(ρ(x)) = -s1(x)  (odd under reflection).
```

Then the weighted sum cancels:

```text
Σ_{x∈Z_12} w(x) s1(x) = 0.
```

**Proof.** Pair each `x` with `ρ(x)`.
Using evenness of `w` and oddness of `s1`:

```text
w(ρ(x)) s1(ρ(x)) + w(x) s1(x) = w(x)(-s1(x)) + w(x)s1(x) = 0.
```

All terms cancel in pairs; the fixed points `x=0` and `x=6` also contribute zero since `s1(0)=s1(6)=0`. ∎

On the current exported strict instance, `F456` exports `u_1 = ± s1` (via `u_1_coords_in_c1_s1 = (0, ±1)`), so the same cancellation holds:

```text
Σ_x w(x) u_1(x) = 0
```

for every even weight `w`.

In particular:

1. `ord_Z12(x)` is even under reflection (`ord_Z12(ρ(x))=ord_Z12(x)`), hence cannot distinguish `u_1` from `-u_1` via `Σ_x ord_Z12(x) u_1(x)`,
2. `r_ord(x) ∝ exp(-alpha_geo*ord_Z12(x))` is also even, hence cannot distinguish sign via `Σ_x r_ord(x) u_1(x)`.

Therefore the strict `ord`/`r_ord` reference family does not, by itself, supply a sign-distinction observable of the form `Σ_x w(x)u_1(x)` on `pair1`
for the current exported instance. ∎

## Consequence (next honest requirement for `H37`)

To obtain a strict sign-distinction observable on `pair1` for this axis, one must go beyond even reference weights:

- introduce an explicit reflection-breaking / orientation source, or
- use a different observable class not reducible to `Σ_x w(x)u(x)` with even `w`,

while keeping all such symmetry-breaking premises explicit (no false‑PASS).

## What `N517` does not claim

`N517` does not claim:

1. discharge of `H37`,
2. a sign-sensitive physical orientation datum derived in strict core,
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.
