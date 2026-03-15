# N480 Current First Strict Shannon Element‑Order Reference Cross‑Entropy Cuts `pair1` `O(2)` to Residual `Z2` (Uniqueness) Theorem

Status: `N480_DISCHARGED_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR1_O2_TO_Z2_UNIQUENESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`T165/P439` isolate a concrete escape hatch from the `N463/N464` translation‑invariance obstruction:
use a **non‑translation‑invariant** reference weight on the typed `Z_12` carrier, but **without** introducing a
marked‑direction/generator slot.

This theorem proves a theorem‑level `O(2)`‑cut uniqueness result for one explicit Shannon‑typed objective built from:

1. the strict‑derived Shannon amplitude `alpha_geo_strict_derived_v1 := 4 ln 2` (`F309/N420`),
2. the typed group structure of `Z_12` via the element order `ord_Z12(x)` (`F329`),
3. the `pair1` `O(2)` family from the `n=12` ring Fourier scaffold (`QW-2190/QW-2191`).

It is intentionally narrow: it does **not** decide `T168/T167/T166`, does **not** claim global `QW-2191` discharge, and
does **not** claim ToE closure.

## Strict‑admissible inputs reused

1. `F329` / `QW-2190`
   - typed `Z_12` carrier language + the `n=12` ring Fourier scaffold (`pair1` exists as a 2D degenerate subspace),
2. `F309/N420`
   - strict‑derived Shannon amplitude `alpha_geo_strict_derived_v1 := 4 ln 2`,
3. `N479`
   - `ord_Z12(x)` is `Aut(Z_12)`‑invariant ⇒ **no marked‑direction slot** for references of the form `f(ord_Z12(x))`.

## Setup

### 1) Typed carrier and element order

Work on the typed group `Z_12` from `F329`.

Define:

$$
\operatorname{ord}_{Z_{12}}(x)
:=
\min\{m\ge 1:\ m x \equiv 0 \pmod{12}\}.
$$

Equivalently (`N479`):

$$
\operatorname{ord}_{Z_{12}}(0)=1,\qquad
\operatorname{ord}_{Z_{12}}(x)=\frac{12}{\gcd(x,12)}\quad (x\neq 0).
$$

### 2) Element‑order reference distribution (no marked direction)

Fix `alpha_geo := alpha_geo_strict_derived_v1 = 4 ln 2 > 0`.

Define unnormalized weights:

$$
w(x):=\exp(-\alpha_{\mathrm{geo}}\,\operatorname{ord}_{Z_{12}}(x)).
$$

Define the normalized reference distribution:

$$
r_{\mathrm{ord}}(x):=\frac{w(x)}{\sum_{y\in Z_{12}} w(y)}.
$$

By `N479`, `ord_Z12(x)` is `Aut(Z_12)`‑invariant, hence this reference does **not** fix a generator/direction on `Z_12`.
It remains **non‑translation‑invariant** on the regular action (it distinguishes the identity orbit `{0}`); this is an
explicit part of the symmetry‑breaking premise class.

### 3) `pair1` `O(2)` family and induced site‑probability distribution

Let `(c_1,s_1)` be any fixed orthonormal basis of the `pair1` 2D mode subspace on `n=12`.

Parameterize the `O(2)` orbit by a rotation angle $\theta$:

$$
u_\theta := \cos\theta\,c_1+\sin\theta\,s_1.
$$

Define the squared‑amplitude site distribution on the 12‑slot scaffold:

$$
p_\theta(x):=\frac{u_\theta(x)^2}{\sum_{y\in Z_{12}} u_\theta(y)^2}.
$$

This is the same “site probability distribution” functional class used in `N463/N464`; the only new ingredient below
will be the non‑translation‑invariant reference `r_{\mathrm{ord}}`.

For the `n=12` ring Fourier scaffold (`QW-2190`), one may choose the standard representative basis (real mode pair)
so that for $x\in\{0,\dots,11\}$:

$$
u_\theta(x)=\cos\!\left(\frac{2\pi}{12}x-\theta\right).
$$

Then $\sum_x u_\theta(x)^2=6$ and:

$$
p_\theta(x)=\frac{1}{6}\cos^2\!\left(\frac{2\pi}{12}x-\theta\right).
$$

## Shannon‑typed objective (cross‑entropy to the element‑order reference)

Define the cross‑entropy objective:

$$
J_{\mathrm{ord}}(\theta):= -\sum_{x\in Z_{12}} p_\theta(x)\,\log r_{\mathrm{ord}}(x).
$$

Since $\log r_{\mathrm{ord}}(x)= -\alpha_{\mathrm{geo}}\operatorname{ord}(x)-\log Z$ with $Z=\sum_y w(y)$, we have:

$$
J_{\mathrm{ord}}(\theta)
=
\alpha_{\mathrm{geo}}\sum_{x\in Z_{12}} p_\theta(x)\,\operatorname{ord}(x)
\log Z.
$$

So the $\theta$‑dependence is entirely in the **linear** expectation
$\mathbb{E}_\theta[\operatorname{ord}] := \sum_x p_\theta(x)\operatorname{ord}(x)$.

## Theorem (nontriviality and `Z2`‑unique minimizer)

### Claim 1. `J_ord(θ)` is not constant on the `O(2)` family.

Using $\cos^2 t = \frac{1+\cos(2t)}{2}$:

$$
\mathbb{E}_\theta[\operatorname{ord}]
=
\frac{1}{6}\sum_{x=0}^{11} \operatorname{ord}(x)\,\cos^2\!\left(\frac{2\pi}{12}x-\theta\right)
$$

$$
=
\frac{1}{12}\sum_{x=0}^{11}\operatorname{ord}(x)
\frac{1}{12}\sum_{x=0}^{11}\operatorname{ord}(x)\cos\!\left(\frac{\pi}{3}x-2\theta\right).
$$

Define the mode‑2 Fourier coefficient (same phase convention as `N467`):

$$
F_2(\operatorname{ord})
:=
\sum_{x=0}^{11}\operatorname{ord}(x)\,e^{i\pi x/3}.
$$

Then:

$$
\sum_{x=0}^{11}\operatorname{ord}(x)\cos\!\left(\frac{\pi}{3}x-2\theta\right)
=
\Re\!\left(e^{-i2\theta}F_2(\operatorname{ord})\right).
$$

So:

$$
\mathbb{E}_\theta[\operatorname{ord}]
=
C_0+\frac{1}{12}\Re\!\left(e^{-i2\theta}F_2(\operatorname{ord})\right),
$$

with constant $C_0=\frac{1}{12}\sum_x\operatorname{ord}(x)$.

It remains to show $F_2(\operatorname{ord})\neq 0$.

Compute `ord` explicitly on `Z_12` (`N479`):

```text
ord(0)=1,
ord(6)=2,
ord(4)=ord(8)=3,
ord(3)=ord(9)=4,
ord(2)=ord(10)=6,
ord(1)=ord(5)=ord(7)=ord(11)=12.
```

Use opposite‑pair reduction (`x` and `x+6` have the same phase factor $e^{i\pi x/3}$):

Let $S_k:=\operatorname{ord}(k)+\operatorname{ord}(k+6)$ for $k=0,\dots,5$. Then:

```text
S0 = 1+2  = 3
S1 = 12+12 = 24
S2 = 6+3   = 9
S3 = 4+4   = 8
S4 = 3+6   = 9
S5 = 12+12 = 24
```

and:

$$
F_2(\operatorname{ord})=\sum_{k=0}^{5}S_k\,e^{i\pi k/3}.
$$

Therefore:

$$
\Re(F_2)= (S_0-S_3)+\frac12(S_1-S_2-S_4+S_5)
= (3-8)+\frac12(24-9-9+24)= -5+\frac12\cdot 30 = 10,
$$

$$
\Im(F_2)=\frac{\sqrt3}{2}(S_1+S_2-S_4-S_5)=\frac{\sqrt3}{2}(24+9-9-24)=0.
$$

So:

$$
F_2(\operatorname{ord})=10\neq 0.
$$

Hence $\mathbb{E}_\theta[\operatorname{ord}]$ (and therefore $J_{\mathrm{ord}}(\theta)$) is not constant on the `O(2)`
family. ∎

### Claim 2. `J_ord(θ)` has a `Z2`‑unique global minimizer on `θ ∈ [0,2π)`.

Since $F_2(\operatorname{ord})=10$ is real:

$$
\mathbb{E}_\theta[\operatorname{ord}] = C_0 + \frac{10}{12}\cos(2\theta).
$$

Thus:

$$
J_{\mathrm{ord}}(\theta)= C_1 + \alpha_{\mathrm{geo}}\frac{10}{12}\cos(2\theta),
$$

for a constant $C_1$.

Because $\alpha_{\mathrm{geo}}>0$ and $\frac{10}{12}>0$, the global minima occur exactly at:

$$
\cos(2\theta)=-1
\quad\Longleftrightarrow\quad
\theta=\frac{\pi}{2}\ (\mathrm{mod}\ \pi).
$$

So the minimizer set on $[0,2\pi)$ is:

$$
\{\theta^\*,\theta^\*+\pi\}
\quad\text{with}\quad
\theta^\*=\frac{\pi}{2}.
$$

This is exactly a residual `Z2` ambiguity (sign flip $u_{\theta+\pi}=-u_\theta$ yields the same $p_\theta$). ∎

## Consequence for `T165` (what is discharged and what is not)

1. This provides an explicit Shannon‑typed objective `J_ord` satisfying:
   - **no marked direction slot** (by `N479`),
   - **nontriviality** (Claim 1),
   - **theorem‑level uniqueness up to residual `Z2`** (Claim 2).
2. The ingredient is **not kernel‑alone**: the non‑translation‑invariant reference `r_ord` is an additional structural
   datum derived from the typed `Z_12` group (it marks the identity orbit). This is exactly the intended “symmetry
   breaking ingredient” class, and it must remain explicit.

## What N480 does not prove

`N480` does not prove:

1. strict-derived per‑site vacuum/self‑coupling arrays (`T168`),
2. strict-derived sigma opposite‑pair sums (`T167`) or the diagonal mode‑2 defect decision (`T166`),
3. strict-core selector closure or global `QW-2191` discharge,
4. ToE closure.

