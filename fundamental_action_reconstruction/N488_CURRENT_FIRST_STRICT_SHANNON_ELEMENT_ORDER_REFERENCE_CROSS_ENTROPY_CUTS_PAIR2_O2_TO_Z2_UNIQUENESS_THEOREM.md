# N488 Current First Strict Shannon Element‑Order Reference Cross‑Entropy Cuts `pair2` `O(2)` to Residual `Z2` (Uniqueness) Theorem

Status: `N488_DISCHARGED_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR2_O2_TO_Z2_UNIQUENESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`N480` proves a theorem‑level `O(2) -> Z2` cut on `pair1` for the strict Shannon element‑order reference
cross‑entropy objective (`F446` / `T165`).

This theorem extends the same **direction‑free** reference datum to `pair2 = span{c2,s2}` on the strict `n=12`
Fourier scaffold (`QW-2190`), so the same strict Shannon ingredient can supply a second phase coordinate
`theta_2` (up to residual `Z2`) without introducing any `eps` / `delta_d` selector slots.

This theorem is intentionally narrow:

1. it does **not** claim sigma‑int corridor strict‑core upgrade (`T159`) by itself,
2. it does **not** discharge `T168/T169` nor decide any diagonal/local defect,
3. it does **not** claim global `QW-2191` discharge or ToE closure.

## Strict‑admissible inputs reused

1. `F329` / `QW-2190`
   - typed `Z_12` scaffold + real Fourier pair planes (`pair_m` exist as degenerate 2D subspaces),
2. `F309/N420`
   - strict‑derived Shannon amplitude `alpha_geo_strict_derived_v1 := 4 ln 2`,
3. `N479`
   - `ord_Z12(x)` is `Aut(Z_12)`‑invariant ⇒ no marked direction slot for references of the form `f(ord_Z12(x))`,
4. `F446`
   - element‑order reference datum `r_ord(x) ∝ exp(-alpha_geo * ord_Z12(x))` (direction‑free; not translation‑invariant).

## Setup (pair2)

Work on `n=12` and the `pair2` `O(2)` orbit.

Let `(c_2,s_2)` be an orthonormal basis of `pair2` on the strict `QW-2190` scaffold.

Parameterize the `O(2)` family by `θ ∈ [0,2π)`:

$$
u_{2,\theta} := \cos\theta\,c_2+\sin\theta\,s_2.
$$

Define the squared site‑amplitude distribution:

$$
p_{2,\theta}(x):=\frac{u_{2,\theta}(x)^2}{\sum_{y\in Z_{12}} u_{2,\theta}(y)^2}.
$$

Define the strict element‑order reference distribution (as in `F446`):

$$
r_{\mathrm{ord}}(x)\propto \exp(-\alpha_{\mathrm{geo}}\,\operatorname{ord}_{Z_{12}}(x)),
\qquad \alpha_{\mathrm{geo}}:=4\ln 2>0.
$$

Define the Shannon‑typed objective (cross‑entropy to `r_ord`):

$$
J_{\mathrm{ord},2}(\theta):=-\sum_{x\in Z_{12}} p_{2,\theta}(x)\,\log r_{\mathrm{ord}}(x).
$$

As in `N480`, since $\log r_{\mathrm{ord}}(x)= -\alpha_{\mathrm{geo}}\operatorname{ord}(x)-\log Z$,
the $\theta$‑dependence reduces to the linear expectation:

$$
J_{\mathrm{ord},2}(\theta) = C + \alpha_{\mathrm{geo}}\sum_{x\in Z_{12}} p_{2,\theta}(x)\,\operatorname{ord}(x),
$$
for a constant `C`.

So nontriviality/uniqueness is controlled by the Fourier coefficient
`F_4(ord) = F_{2m}(ord)` with `m=2`.

## Theorem (nontriviality and `Z2`‑unique minimizer on pair2)

### Claim 1. `J_ord,2(θ)` is not constant on the `pair2` `O(2)` family.

On the `QW-2190` ring scaffold one may take:

$$
u_{2,\theta}(x)=\cos\!\left(\frac{4\pi}{12}x-\theta\right),
$$

so (up to the normalization constant) the same identity as in `N480` yields:

$$
\sum_{x=0}^{11} p_{2,\theta}(x)\,\operatorname{ord}(x)
=
C_0 + \frac{1}{12}\Re\!\left(e^{-i2\theta}F_4(\operatorname{ord})\right),
$$

where:

$$
F_4(\operatorname{ord})
:=
\sum_{x=0}^{11}\operatorname{ord}(x)\,e^{i\frac{8\pi}{12}x}
=
\sum_{x=0}^{11}\operatorname{ord}(x)\,e^{i\frac{2\pi}{3}x}.
$$

It remains to show `F_4(ord) ≠ 0`.

Compute `ord` on `Z_12` (`N479`):

```text
ord(0)=1,
ord(6)=2,
ord(4)=ord(8)=3,
ord(3)=ord(9)=4,
ord(2)=ord(10)=6,
ord(1)=ord(5)=ord(7)=ord(11)=12.
```

Use opposite‑pair sums `S_k := ord(k)+ord(k+6)` for `k=0..5`:

```text
S0 = 3,  S1 = 24,  S2 = 9,  S3 = 8,  S4 = 9,  S5 = 24.
```

Since `e^{i(2π/3)(k+6)} = e^{i(2π/3)k}`, we have:

$$
F_4(\operatorname{ord})=\sum_{k=0}^{5} S_k\,e^{i\frac{2\pi}{3}k}.
$$

Now `cos(2πk/3) ∈ {1,-1/2}` and `sin(2πk/3) ∈ {0,±\sqrt3/2}`.
Compute:

$$
\Re(F_4) = (S_0+S_3) - \frac12(S_1+S_2+S_4+S_5)
= (3+8) - \frac12(24+9+9+24) = 11 - 33 = -22,
$$

$$
\Im(F_4) = \frac{\sqrt3}{2}(S_1 - S_2 + S_4 - S_5)
=\frac{\sqrt3}{2}(24-9+9-24) = 0.
$$

Therefore:

$$
F_4(\operatorname{ord})=-22\neq 0.
$$

Hence `J_ord,2(θ)` is not constant on the `pair2` `O(2)` family. ∎

### Claim 2. `J_ord,2(θ)` has a `Z2`‑unique global minimizer set on `θ ∈ [0,2π)`.

Since `F_4(ord)` is real, we have:

$$
\sum_x p_{2,\theta}(x)\,\operatorname{ord}(x)
=
C_0 + \frac{-22}{12}\cos(2\theta).
$$

Because `alpha_geo > 0`, minimizing `J_ord,2` is equivalent to minimizing that expectation, hence to maximizing
`cos(2θ)`.

Therefore the global minimizers are exactly:

$$
\cos(2\theta)=1
\quad\Longleftrightarrow\quad
\theta=0\ (\mathrm{mod}\ \pi),
$$

so the minimizer set on `[0,2π)` is:

$$
\{0,\pi\}
=
\{0\ (\mathrm{mod}\ \pi)\}.
$$

This is exactly residual `Z2` (sign flip `u_{θ+π}=-u_θ` leaves the squared‑amplitude distribution invariant). ∎

## Consequence (what this enables and what it does not)

This theorem provides a theorem‑level, direction‑free `pair2` `O(2) -> Z2` cut for the strict Shannon element‑order
reference objective.

It does **not** by itself:

1. discharge sigma‑int corridor strict‑core upgrade (`T159`),
2. discharge `T162` (slot‑free sigma‑int → theta construction class),
3. discharge `T168/T169` nor any ToE closure target.

