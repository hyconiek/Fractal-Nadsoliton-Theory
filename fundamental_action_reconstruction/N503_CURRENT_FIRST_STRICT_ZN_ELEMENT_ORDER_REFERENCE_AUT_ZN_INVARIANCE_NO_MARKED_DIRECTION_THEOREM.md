# N503 Current First Strict `Z_n` Element‑Order Reference `Aut(Z_n)`‑Invariance ⇒ No Marked‑Direction Slot Theorem

Status: `N503_DISCHARGED_CURRENT_FIRST_STRICT_ZN_ELEMENT_ORDER_REFERENCE_AUT_ZN_INVARIANCE_NO_MARKED_DIRECTION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`N479` established (for `n=12`) that the element-order profile `ord_{Z_12}(x)` is `Aut(Z_12)`‑invariant and therefore does not
smuggle an untracked **marked-direction / generator** choice slot.

`N503` records the corresponding **general cyclic-group fact**:

```text
For any n≥2, ord_{Z_n}(x) is invariant under Aut(Z_n).
Therefore any reference weights of the form f(ord_{Z_n}(x)) are direction-free under Aut(Z_n).
```

This theorem is intentionally narrow. It does **not** claim any theorem-level `O(2) -> Z2` cut, does not claim any physical
admissibility of such a reference, and does not discharge `QW-2191`.

## Strict-admissible inputs reused

Only elementary finite-group arithmetic.

(No typed `Phase_n`/`Aut(Z_n)` infrastructure is required for the core proof.)

## Setup

Let `n≥2` and let `Z_n` be the additive cyclic group.

Let `Aut(Z_n)` be the unit group `(Z/nZ)^×`, acting on elements by multiplication:

$$
u\cdot x := (u x)\bmod n,\qquad u\in (Z/nZ)^{\times}.
$$

Define the element order:

$$
\operatorname{ord}_{Z_n}(x)
:=
\min\{m\ge 1:\ m x \equiv 0 \pmod{n}\}.
$$

Equivalently:

$$
\operatorname{ord}_{Z_n}(0)=1,\qquad
\operatorname{ord}_{Z_n}(x)=\frac{n}{\gcd(x,n)}\quad (x\neq 0).
$$

## Theorem-level claims

### Claim 1. `ord_{Z_n}` is `Aut(Z_n)`‑invariant.

For any unit $u\in (Z/nZ)^{\times}$ and any $x\in Z_n$:

$$
\operatorname{ord}_{Z_n}(u\cdot x)=\operatorname{ord}_{Z_n}(x).
$$

**Proof.** Since $u$ is a unit modulo $n$, $\gcd(u,n)=1$.
Hence $\gcd(ux,n)=\gcd(x,n)$:

1. $\gcd(x,n)\mid ux$ and $\gcd(x,n)\mid n$, so $\gcd(x,n)\mid \gcd(ux,n)$,
2. if $d\mid ux$ and $d\mid n$ with $\gcd(u,n)=1$, then $d\mid x$ as well (because $u$ is invertible modulo $d$),
   so $\gcd(ux,n)\mid \gcd(x,n)$.

Therefore $\gcd(ux,n)=\gcd(x,n)$, and the formula $\operatorname{ord}_{Z_n}(x)=n/\gcd(x,n)$ yields the claim (with the trivial $x=0$ case). ∎

### Claim 2. Any reference weights of the form `w_x := f(ord_{Z_n}(x))` are `Aut(Z_n)`‑invariant.

Immediate from Claim 1: for all units $u$,

$$
w_{u\cdot x}=f(\operatorname{ord}(u\cdot x))=f(\operatorname{ord}(x))=w_x.
$$

In particular, the Shannon-type reference weights

$$
w_x := \exp(-\alpha\,\operatorname{ord}_{Z_n}(x))
$$

are `Aut(Z_n)`‑invariant for any real $\alpha$ (including the strict-side source value
`alpha_geo_strict_derived_v1 := 4 ln 2` when used as the coefficient in the strict Shannon lane).

## Consequence (what this does and does not buy)

1. **No marked-direction slot:** `ord_{Z_n}` (and any `f(ord_{Z_n})` reference) does not introduce a generator/direction
   choice; it is `Aut(Z_n)`‑invariant.
2. **Marked-site remains:** the reference is still **not translation-invariant** on the regular action scaffold; it
   distinguishes the identity orbit `{0}`. Any strict admissibility story must therefore still track this as an explicit
   non-translation-invariant datum.

## What N503 does not prove

`N503` does not prove:

1. that `f(ord_{Z_n})` is strict-core admissible as physics on any `n`,
2. that any defect `F_{2m}(ord_{Z_n})` is nonzero (even for a fixed `n`),
3. that a KL/Shannon objective using this reference has a theorem-level unique minimizer on any pair plane,
4. any discharge of `QW-2191`,
5. strict-core selector closure / admissible `S_sel_int`,
6. ToE closure.

