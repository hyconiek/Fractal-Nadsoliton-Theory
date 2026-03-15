# N513 Current First Strict Shannon Element‑Order Reference `Z_24` Cross‑Entropy Cuts All Fourier‑Degenerate Pairs `O(2)` to Residual `Z2` (Uniqueness) Theorem

Status: `N513_DISCHARGED_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_Z24_CROSS_ENTROPY_CUTS_ALL_FOURIER_DEGENERATE_PAIRS_O2_TO_Z2_UNIQUENESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

`P461/P462/F458` introduced a cautious scope-extension track for the Shannon element‑order reference construction beyond the physical `n=12` scaffold.

`F468` now exports a **typed** `Z_24` element‑order reference datum `r_ord` (shape stored as `ord_{Z_24}(x)`) and an explicit real Fourier mode-index assignment object on `Z_24`.

This theorem packages the strict statement that, on `Z_24`, the same element‑order reference cross‑entropy objective cuts the continuous `O(2)` families on **all**
Fourier-degenerate pair planes down to residual `Z2` (the unavoidable sign flip).

This is a scope-extension statement only. It does **not** promote `n=24` into the strict physical `QW-2190` scaffold, does not discharge global `QW-2191`, and does not claim selector closure or ToE closure.

## Strict‑admissible inputs reused

1. `F458`
   - typed `Z_24` carrier + regular action (`I_24_v1`, `Z_24_v1`, `tau_Z24_v1`),
2. `F309/N420`
   - strict-derived Shannon amplitude `alpha_geo_strict_derived_v1 := 4 ln 2`,
3. `N503`
   - `ord_{Z_n}` is `Aut(Z_n)`‑invariant for any `n` ⇒ no marked-direction slot for references of the form `f(ord_{Z_n}(x))`,
4. `F468`
   - exports `ord_{Z_24}(x)` as the diagonal reference profile and exports an explicit `Z_24` mode-index assignment basis object with recorded defects `F_{2m}(ord)`.

## Setup (typed `Z_24` Fourier pair planes)

Work on the real Fourier scaffold on `Z_24`:

- nondegenerate modes: `e0`, `e12`,
- degenerate pair planes: `pair_m = span{c_m,s_m}` for `m=1..11`.

Define the strict element-order reference weights:

$$
r_{\mathrm{ord}}(x)\propto \exp(-\alpha_{\mathrm{geo}}\,\operatorname{ord}_{Z_{24}}(x)),
\qquad \alpha_{\mathrm{geo}} := 4\ln 2 > 0.
$$

For each `m=1..11`, parameterize the `O(2)` family on `pair_m` by `θ∈[0,2π)`:

$$
u_{m,\theta} := \cos\theta\,c_m+\sin\theta\,s_m,
$$

and define the squared site‑amplitude distribution `p_{m,θ}(x) ∝ u_{m,θ}(x)^2` as usual.

Define the cross‑entropy objective:

$$
J_{\mathrm{ord},m}(\theta):=-\sum_{x\in Z_{24}} p_{m,\theta}(x)\,\log r_{\mathrm{ord}}(x).
$$

As in the `n=12` proofs, since $\log r_{\mathrm{ord}}(x)= -\alpha_{\mathrm{geo}}\operatorname{ord}(x)+\mathrm{const}$, minimizing `J_ord,m` is equivalent to minimizing
the quadratic-form expectation of `ord_{Z_24}` on the `pair_m` family.

## Theorem (all `pair_m` are cut to residual `Z2` on `Z_24`)

Define, for each `m=1..11`, the defect:

$$
F_{2m}(\operatorname{ord})
:=
\sum_{x=0}^{23}\operatorname{ord}_{Z_{24}}(x)\,e^{i\frac{4\pi m}{24}x}.
$$

`F468` records the evaluated defects on the exported strict `Z_24` element-order profile as (within numerical tolerance; imaginary parts are near zero on the exported instance):

```text
F2  ≈  10
F4  ≈  42
F6  ≈ -35
F8  ≈ -86
F10 ≈  10
F12 ≈ -147
F14 ≈  10
F16 ≈ -86
F18 ≈ -35
F20 ≈  42
F22 ≈  10
```

In particular, `F_{2m}(\operatorname{ord}) ≠ 0` for every `m=1..11`. Therefore, for each pair plane `pair_m`, the objective `J_ord,m(θ)` is not constant and its minimizer set
on `[0,2π)` has exactly two points differing by `π` (the unavoidable sign flip `u_{θ+π}=-u_{θ}` leaving `p_{m,θ}` unchanged).

Equivalently, on `Z_24` the Shannon element‑order reference cross‑entropy objective cuts each continuous `O(2)` family down to residual `Z2` on every Fourier-degenerate pair plane. ∎

## Consequence

The exported `Z_24` assignment object

```text
ModeIndexAssignment_shannon_element_order_reference_z24_strict_core_v1
```

from `F468` is therefore well-defined in the strict scope-extension sense: it canonically fixes axes on each `pair_m (m=1..11)` **up to residual sign**,
without introducing a marked-direction slot (`N503`) and without any physical promotion of `n=24`.

## What `N513` does not claim

`N513` does not claim:

1. any strict-core promotion of `n=24` into the physical `QW-2190` scaffold,
2. any global discharge of `QW-2191`,
3. strict-core selector closure / admissible `S_sel_int`,
4. any sign-sensitive physical orientation datum,
5. ToE closure.
