# N496 Current First Strict Shannon Element‑Order Reference Cross‑Entropy Cuts `pair3..pair5` `O(2)` to Residual `Z2` (Uniqueness) Theorem

Status: `N496_DISCHARGED_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR3_TO_PAIR5_O2_TO_Z2_UNIQUENESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`F446/N480` export and discharge a strict Shannon-typed `O(2)->Z2` cut on `pair1` for the element‑order reference
cross‑entropy objective (`T165`).

`N488` extends the same direction‑free reference datum to `pair2`.

This theorem extends the same strict reference/objective family to the remaining degenerate Fourier pairs on the strict
`n=12` scaffold:

```text
pair3 = span{c3,s3},  pair4 = span{c4,s4},  pair5 = span{c5,s5}.
```

It does **not** claim sigma-int corridor upgrade (`T159`) by itself, does **not** claim strict-core selector closure,
does **not** claim global `QW-2191` discharge, and does **not** claim ToE closure.

## Strict‑admissible inputs reused

1. `F329` / `QW-2190`
   - typed `Z_12` scaffold + real Fourier pair planes,
2. `F309/N420`
   - strict-derived Shannon amplitude `alpha_geo_strict_derived_v1 := 4 ln 2`,
3. `N479`
   - `ord_Z12(x)` is `Aut(Z_12)`‑invariant ⇒ no marked direction slot for references of the form `f(ord_Z12(x))`,
4. `F446`
   - strict element‑order reference datum `r_ord(x) ∝ exp(-alpha_geo * ord_Z12(x))`,
5. `N480`
   - proof template reducing the cross-entropy objective on `pair_m` to a cosine term with Fourier defect `F_{2m}(ord)`,
6. `N488`
   - executed example of the same reduction and explicit `F_4(ord)` computation.

## Setup (pair m = 3,4,5)

Work on `n=12` and let `(c_m,s_m)` be an orthonormal basis of `pair_m` on the strict `QW-2190` scaffold.

Parameterize the `O(2)` family by `θ ∈ [0,2π)`:

$$
u_{m,\theta} := \cos\theta\,c_m+\sin\theta\,s_m.
$$

Define the squared site‑amplitude distribution:

$$
p_{m,\theta}(x):=\frac{u_{m,\theta}(x)^2}{\sum_{y\in Z_{12}} u_{m,\theta}(y)^2}.
$$

Define the strict element‑order reference distribution (as in `F446`):

$$
r_{\mathrm{ord}}(x)\propto \exp(-\alpha_{\mathrm{geo}}\,\operatorname{ord}_{Z_{12}}(x)),
\qquad \alpha_{\mathrm{geo}}:=4\ln 2>0.
$$

Define the Shannon‑typed objective (cross‑entropy to `r_ord`):

$$
J_{\mathrm{ord},m}(\theta):=-\sum_{x\in Z_{12}} p_{m,\theta}(x)\,\log r_{\mathrm{ord}}(x).
$$

As in `N480/N488`, since $\log r_{\mathrm{ord}}(x)= -\alpha_{\mathrm{geo}}\operatorname{ord}(x) + \mathrm{const}$,
minimizing `J_ord,m` is equivalent to minimizing:

$$
E_m(\theta) := \sum_{x\in Z_{12}} p_{m,\theta}(x)\,\operatorname{ord}(x).
$$

## Claim 1. Reduction to `F_{2m}(ord)` and nontriviality on `pair3..pair5`.

As in `N480/N488`, one obtains:

$$
E_m(\theta) = C_0 + \frac{1}{12}\Re\!\left(e^{-i2\theta}F_{2m}(\operatorname{ord})\right),
$$

where:

$$
F_{2m}(\operatorname{ord})
:=
\sum_{x=0}^{11}\operatorname{ord}(x)\,e^{i\frac{4\pi m}{12}x}.
$$

It remains to show `F_{2m}(ord) ≠ 0` for `m=3,4,5`.

Using the explicit element orders on `Z_12` (as in `N488`), compute:

### `m=3` (`F_6`)

Since $e^{i(4\pi\cdot 3/12)x} = e^{i\pi x} = (-1)^x$:

$$
F_6(\operatorname{ord}) = \sum_{x=0}^{11}\operatorname{ord}(x)\,(-1)^x
= -35 \neq 0.
$$

### `m=4` (`F_8`)

Since $e^{i(4\pi\cdot 4/12)x} = e^{i\frac{4\pi}{3}x} = \overline{e^{i\frac{2\pi}{3}x}}$ and `ord(x)` is real:

$$
F_8(\operatorname{ord}) = \overline{F_4(\operatorname{ord})} = F_4(\operatorname{ord}) = -22 \neq 0,
$$

using `N488`’s explicit computation `F_4(ord)=-22`.

### `m=5` (`F_10`)

Similarly $e^{i(4\pi\cdot 5/12)x} = e^{i\frac{5\pi}{3}x} = \overline{e^{i\frac{\pi}{3}x}}$, hence:

$$
F_{10}(\operatorname{ord}) = \overline{F_2(\operatorname{ord})} = F_2(\operatorname{ord}) = 10 \neq 0,
$$

using `N480`’s nontriviality (`F_2(ord)≠0`) and the explicit value recorded in that proof family.

Therefore, for each `m=3,4,5`, the objective is not constant and cuts the continuous `O(2)` family down to residual `Z2`.
∎

## Claim 2. `Z2`‑unique global minimizer sets on `pair3..pair5`.

In each case above, `F_{2m}(ord)` is real, so:

$$
E_m(\theta) = C_0 + \frac{F_{2m}(\operatorname{ord})}{12}\cos(2\theta).
$$

Since `alpha_geo>0`, minimizing `J_ord,m` is equivalent to minimizing `E_m`.

Therefore:

1. for `m=3`: `F_6(ord)=-35<0` ⇒ minimizers are `θ = 0 (mod π)`,
2. for `m=4`: `F_8(ord)=-22<0` ⇒ minimizers are `θ = 0 (mod π)`,
3. for `m=5`: `F_10(ord)=10>0` ⇒ minimizers are `θ = π/2 (mod π)`.

On `[0,2π)`, each minimizer set has exactly two points differing by `π`, i.e. residual `Z2` (the unavoidable sign flip
`u_{θ+π}=-u_{θ}` leaves the squared‑amplitude distribution invariant). ∎

## Consequence

Combining:

1. `N480` (`pair1` cut),
2. `N488` (`pair2` cut),
3. this theorem (`pair3..pair5` cuts),

the strict Shannon element‑order reference cross‑entropy objective cuts `O(2)` down to residual `Z2` on **all**
Fourier-degenerate pair planes `pair_m (m=1..5)` on the strict `n=12` scaffold.

This enables exporting a full strict-core mode-index assignment basis object on the Shannon lane (`F454`), without
importing per-site vacuum/self-coupling providers (`T168/T169`) and without implying global selector closure.

## What `N496` does not claim

`N496` does not claim:

1. sigma-int corridor strict-core upgrade (`T159`) by itself,
2. strict-core selector closure or admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.

