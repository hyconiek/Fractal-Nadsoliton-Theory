# N479 Current First Strict `Z_12` Element‑Order Reference `Aut(Z_12)`‑Invariance ⇒ No Marked‑Direction Slot Theorem

Status: `N479_DISCHARGED_CURRENT_FIRST_STRICT_Z12_ELEMENT_ORDER_REFERENCE_AUT_Z12_INVARIANCE_NO_MARKED_DIRECTION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`T165/P439` highlight one concrete candidate class of non‑translation‑invariant reference data on the `n=12` carrier:

```text
r(x) ∝ exp(-alpha_geo * ord_Z12(x)).
```

One sharp hygiene question is whether this reference shape smuggles an untracked **marked-direction / generator choice**
slot.

`N479` records a strict structural fact:

```text
ord_Z12(x) is invariant under Aut(Z_12) and therefore does not fix a generator/direction.
So reference weights of the form f(ord_Z12(x)) are Aut(Z_12)-invariant and are “direction-free” on Z_12.
```

This theorem is intentionally narrow. It does **not** claim the reference is already strict‑core admissible as physics,
and it does not discharge `T165`.

## Strict-admissible inputs reused

1. `QW-2190` / `F329`
   - typed `Z_12` carrier language,
2. `N455`
   - explicit `Aut(Z_12)` orbit decomposition (used only for the orbit listing; not required for the core proof).

The proof itself is elementary finite‑group arithmetic.

## Setup

Let `n=12` and let `Z_12` be the additive cyclic group.

Let `Aut(Z_12)` be the unit group `(Z/12Z)^× = {1,5,7,11}` acting on elements by multiplication:

$$
u\cdot x := (u x)\bmod 12.
$$

Define the element order:

$$
\operatorname{ord}_{Z_{12}}(x)
:=
\min\{m\ge 1:\ m x \equiv 0 \pmod{12}\}.
$$

Equivalently:

$$
\operatorname{ord}_{Z_{12}}(0)=1,\qquad
\operatorname{ord}_{Z_{12}}(x)=\frac{12}{\gcd(x,12)}\quad (x\neq 0).
$$

## Theorem-level claims

### Claim 1. `ord_Z12` is `Aut(Z_12)`‑invariant.

For any unit $u\in (Z/12Z)^{\times}$ and any $x\in Z_{12}$:

$$
\operatorname{ord}_{Z_{12}}(u\cdot x)=\operatorname{ord}_{Z_{12}}(x).
$$

**Proof.** Since $u$ is a unit modulo $12$, $\gcd(u,12)=1$, hence:

$$
\gcd(u x,12)=\gcd(x,12).
$$

Using $\operatorname{ord}_{Z_{12}}(x)=12/\gcd(x,12)$ for $x\neq 0$ (and the trivial $x=0$ case), we obtain the claim. ∎

### Claim 2. Any reference weights of the form `w_x := f(ord_Z12(x))` are `Aut(Z_12)`‑invariant.

Immediate from Claim 1: for all units $u$,

$$
w_{u\cdot x}=f(\operatorname{ord}(u\cdot x))=f(\operatorname{ord}(x))=w_x.
$$

In particular, the weights used in `P439`:

$$
w_x := \exp(-\alpha_{\mathrm{geo}}\operatorname{ord}_{Z_{12}}(x)),
$$

are `Aut(Z_12)`‑invariant for any real $\alpha_{\mathrm{geo}}$ (including the strict source value
`alpha_geo_strict_derived_v1 := 4 ln 2` from `N420`).

## Orbit form on `n=12` (explicit)

By `N455`, the `Aut(Z_12)` orbits in `Z_12` are:

$$
\{0\},\ \{6\},\ \{1,5,7,11\},\ \{2,10\},\ \{3,9\},\ \{4,8\}.
$$

Since `ord_Z12(x)` depends only on `gcd(x,12)`, it is constant on these orbits:

```text
ord(0)=1,
ord(6)=2,
ord(4)=ord(8)=3,
ord(3)=ord(9)=4,
ord(2)=ord(10)=6,
ord(1)=ord(5)=ord(7)=ord(11)=12.
```

So the element‑order reference is quotient‑safe under the `Aut(Z_12)` orbit discipline: it does not pick a generator.

## Consequence (what this does and does not buy)

1. **No marked-direction slot:** `ord_Z12` (and any `f(ord_Z12)` reference) does not introduce a generator/direction
   choice; it is `Aut(Z_12)`‑invariant.
2. **Marked-site remains:** the reference is still **not translation-invariant** on the 12‑site regular action (it
   distinguishes the identity orbit `{0}`), so any strict-core admissibility story must still explicitly account for
   that non‑translation‑invariant datum (cf. `N464` / `T165`).

## What N479 does not prove

`N479` does not prove:

1. that any non‑translation‑invariant reference datum is strict‑core admissible as physics,
2. that a KL/Shannon objective using this reference has a theorem‑level unique minimizer (discharge of `T165`),
3. any strict-derived per-site vacuum/self-coupling provider (`T168`) or sigma six‑sum instantiation (`T167`),
4. discharge of `T166` / `QW-2191`,
5. ToE closure.

