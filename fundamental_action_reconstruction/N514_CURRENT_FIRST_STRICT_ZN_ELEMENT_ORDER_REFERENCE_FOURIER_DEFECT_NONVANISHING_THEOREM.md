# N514 Current First Strict `Z_n` Element‑Order Reference Fourier Defect Nonvanishing Theorem

Status: `N514_DISCHARGED_CURRENT_FIRST_STRICT_ZN_ELEMENT_ORDER_REFERENCE_FOURIER_DEFECT_NONVANISHING_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

The strict Shannon element‑order reference lane uses the defect

$$
F_{2m}(\operatorname{ord})
:=
\sum_{x=0}^{n-1}\operatorname{ord}_{Z_n}(x)\,e^{i\frac{4\pi m}{n}x}
$$

to cut each Fourier‑degenerate `pair_m` continuous `O(2)` family down to residual `Z2` (axis‑only uniqueness).

For `n=12`, this is discharged in the repo on all `pair_m (m=1..5)` (`N480`, `N488`, `N496`). For typed `Z_24`, the same statement is packaged
in a strict scope‑extension sense (`F468`, `N513`). A computational scan (`P461`) shows nonzero defects on multiple `n`, but remained probe‑only.

This theorem removes the remaining scope‑extension **false‑PASS** risk by proving:

```text
for every n≥2 and every k∈{1,…,n−1}, the Fourier defect F_k(ord_{Z_n}) is nonzero.
```

Equivalently, for every `n≥2` and every Fourier‑degenerate `pair_m` on `Z_n`, the defect `F_{2m}(ord_{Z_n})` is nonzero.

This is a purely arithmetic result. It does **not** promote any `n≠12` carrier into the strict physical `QW-2190` scaffold, does **not** discharge
global `QW-2191`, and does **not** claim strict‑core selector closure or ToE closure.

## Strict‑admissible inputs reused

1. Elementary arithmetic of cyclic groups,
2. (Hygiene only) `N503`: `ord_{Z_n}` is `Aut(Z_n)`‑invariant ⇒ no marked‑direction/generator slot.

## Setup

Let `n≥2`. Work on the additive cyclic group `Z_n`.

Define the element order:

$$
\operatorname{ord}_{Z_n}(0)=1,\qquad
\operatorname{ord}_{Z_n}(x)=\frac{n}{\gcd(x,n)}\quad (x\neq 0).
$$

For any integer `k`, define the discrete Fourier defect:

$$
F_k(n)
:=
\sum_{x=0}^{n-1}\operatorname{ord}_{Z_n}(x)\,e^{2\pi i k x/n}.
$$

Note: the repo’s `pair_m` defect convention uses `e^{i(4\pi m/n)x}=e^{2\pi i (2m)x/n}`, so:

$$
F_{2m}(\operatorname{ord}_{Z_n}) = F_{k}(n)\quad\text{with}\quad k:=2m.
$$

Define the Ramanujan sum:

$$
c_q(k)
:=
\sum_{\substack{1\le a\le q\\ \gcd(a,q)=1}} e^{2\pi i a k/q}.
$$

## Lemma 1 (gcd‑class decomposition ⇒ Ramanujan representation)

For all `n≥2` and integers `k`:

$$
F_k(n)=\sum_{q\mid n} q\,c_q(k).
$$

**Proof.** Group the sum over $x\in\{0,\dots,n-1\}$ by $t:=\gcd(x,n)$.
Write $q:=n/t$, so $t=n/q$, and every such $x$ can be written uniquely as $x=(n/q)a$ with $0\le a<q$ and $\gcd(a,q)=1$.
Then $\operatorname{ord}_{Z_n}(x)=n/t=q$ and:

$$
\sum_{\substack{x\bmod n\\ \gcd(x,n)=n/q}} e^{2\pi i kx/n}
\;=\;
\sum_{\substack{0\le a<q\\ \gcd(a,q)=1}} e^{2\pi i k (n/q)a/n}
\;=\;
\sum_{\substack{0\le a<q\\ \gcd(a,q)=1}} e^{2\pi i k a/q}
\;=\; c_q(k).
$$

Multiplying by the common value `q` of `ord(x)` on this gcd‑class and summing over all divisors `q|n` gives the claim. ∎

## Lemma 2 (prime‑power Ramanujan sums)

Let `q=p^t` with prime `p` and `t≥1`. Let `v:=v_p(k)` (the `p`‑adic valuation, with `v_p(0)=+∞`).
Then:

$$
c_{p^t}(k)=
\begin{cases}
\varphi(p^t)=p^t-p^{t-1}, & \text{if } v\ge t,\\
-p^{t-1}, & \text{if } v=t-1,\\
0, & \text{if } v\le t-2.
\end{cases}
$$

**Proof.** Consider the full geometric sum

$$
S:=\sum_{a=0}^{p^t-1} e^{2\pi i a k/p^t}.
$$

If $p^t\mid k$ then every term is $1$ and $S=p^t$. Otherwise $S=0$.

Now split the full sum into “units” and “multiples of $p$”:

$$
S = c_{p^t}(k) + \sum_{b=0}^{p^{t-1}-1} e^{2\pi i (pb)k/p^t}
  = c_{p^t}(k) + \sum_{b=0}^{p^{t-1}-1} e^{2\pi i b k/p^{t-1}}.
$$

Let

$$
S_p:=\sum_{b=0}^{p^{t-1}-1} e^{2\pi i b k/p^{t-1}}.
$$

If $p^{t-1}\mid k$ then $S_p=p^{t-1}$, else $S_p=0$.
Hence $c_{p^t}(k)=S-S_p$, which gives:

- if $p^t\mid k$: $c_{p^t}(k)=p^t-p^{t-1}=\varphi(p^t)$,
- if $p^{t-1}\mid k$ but $p^t\nmid k$: $c_{p^t}(k)=0-p^{t-1}=-p^{t-1}$,
- if $p^{t-1}\nmid k$: $c_{p^t}(k)=0-0=0$.

∎

## Lemma 3 (explicit prime‑power defects are nonzero)

Let `n=p^a` with prime `p` and `a≥1`. Let `v:=v_p(k)`.
Then:

1. If `v < a`, then

$$
F_k(p^a) = -\frac{p^{2v+2}-1}{p+1}\neq 0.
$$

2. If `v \ge a`, then

$$
F_k(p^a) = 1 + (p-1)p\sum_{t=0}^{a-1} p^{2t}\neq 0.
$$

**Proof.** By Lemma 1:

$$
F_k(p^a) = \sum_{t=0}^{a} p^t\,c_{p^t}(k),
$$

with $c_{p^0}(k)=c_1(k)=1$.

- If `v < a`, then by Lemma 2:
  - for `t=1..v`: $c_{p^t}(k)=\varphi(p^t)=p^t-p^{t-1}$,
  - for `t=v+1`: $c_{p^{v+1}}(k)=-p^{v}$,
  - for `t>v+1`: $c_{p^t}(k)=0$.

  Therefore:

  $$
  \begin{aligned}
  F_k(p^a)
  &= 1 + \sum_{t=1}^{v} p^t(p^t-p^{t-1}) + p^{v+1}(-p^{v})\\
  &= 1 + \sum_{t=1}^{v} (p^{2t}-p^{2t-1}) - p^{2v+1}\\
  &= 1 + (p-1)\sum_{t=1}^{v} p^{2t-1} - p^{2v+1}.
  \end{aligned}
  $$

  The geometric sum gives $\sum_{t=1}^{v} p^{2t-1}=p\frac{p^{2v}-1}{p^2-1}$, so:

  $$
  F_k(p^a)=1 + (p-1)p\frac{p^{2v}-1}{p^2-1} - p^{2v+1}
  = 1 + \frac{p(p^{2v}-1)}{p+1} - p^{2v+1}
  = -\frac{p^{2v+2}-1}{p+1}.
  $$

  This is nonzero.

- If `v ≥ a`, then by Lemma 2, for all `t=1..a` one has $c_{p^t}(k)=\varphi(p^t)=p^t-p^{t-1}$, so:

  $$
  \begin{aligned}
  F_k(p^a)
  &= 1 + \sum_{t=1}^{a} p^t(p^t-p^{t-1})
   = 1 + \sum_{t=1}^{a} (p^{2t}-p^{2t-1})\\
  &= 1 + (p-1)\sum_{t=1}^{a} p^{2t-1}
   = 1 + (p-1)p\sum_{t=0}^{a-1}p^{2t},
  \end{aligned}
  $$

  which is strictly positive, hence nonzero.

∎

## Lemma 4 (multiplicativity in `n`)

For fixed integer `k`, the function `n ↦ F_k(n)` is multiplicative:
if `gcd(n_1,n_2)=1` then:

$$
F_k(n_1 n_2)=F_k(n_1)\,F_k(n_2).
$$

**Proof.** First record a standard divisor formula for the Ramanujan sum.
Let `μ` be the Möbius function. Use the identity

$$
\mathbf{1}_{\gcd(a,q)=1}=\sum_{d\mid \gcd(a,q)}\mu(d).
$$

Then:

$$
\begin{aligned}
c_q(k)
&=\sum_{a=0}^{q-1} e^{2\pi i a k/q}\,\mathbf{1}_{\gcd(a,q)=1}
 =\sum_{a=0}^{q-1} e^{2\pi i a k/q}\sum_{d\mid \gcd(a,q)}\mu(d)\\
&=\sum_{d\mid q}\mu(d)\sum_{\substack{0\le a<q\\ d\mid a}} e^{2\pi i a k/q}
 =\sum_{d\mid q}\mu(d)\sum_{b=0}^{q/d-1} e^{2\pi i (db)k/q}\\
&=\sum_{d\mid q}\mu(d)\sum_{b=0}^{q/d-1} e^{2\pi i b k/(q/d)}.
\end{aligned}
$$

The inner geometric sum equals `q/d` if `(q/d) | k` and equals `0` otherwise.
Therefore, writing `e:=q/d`, we obtain:

$$
c_q(k)=\sum_{e\mid \gcd(q,k)} e\,\mu(q/e).
$$

From this divisor formula, if `gcd(q_1,q_2)=1` then `c_{q_1 q_2}(k)=c_{q_1}(k)c_{q_2}(k)` by factoring each divisor
`e|gcd(q_1 q_2,k)` uniquely as `e=e_1 e_2` with `e_i|gcd(q_i,k)` and using multiplicativity of `μ`.

Therefore the function $f(q):=q\,c_q(k)$ is multiplicative in `q`, and the divisor sum from Lemma 1,

$$
F_k(n)=\sum_{q\mid n} f(q),
$$

is multiplicative in `n`. ∎

## Theorem (nonvanishing of `F_k(ord_{Z_n})`)

Let `n≥2` and let `k∈{1,…,n-1}`.
Then:

$$
F_k(n)\neq 0.
$$

**Proof.** Write $n=\prod_i p_i^{a_i}$.
By Lemma 4,

$$
F_k(n)=\prod_i F_k(p_i^{a_i}).
$$

Each prime‑power factor is nonzero by Lemma 3, hence the product is nonzero. ∎

## Corollary (Shannon element‑order reference cuts every Fourier‑degenerate pair plane on every `Z_n`)

For any `n≥2` and any Fourier‑degenerate pair plane `pair_m` on `Z_n` (repo convention `k:=2m`), the defect
`F_{2m}(ord_{Z_n})` is nonzero. Therefore (by the same reduction template used in `N480/N488/N496`),
the element‑order reference cross‑entropy objective on `pair_m`:

$$
J_{\mathrm{ord},m}(\theta):=-\sum_{x\in Z_n} p_{m,\theta}(x)\,\log r_{\mathrm{ord}}(x),
$$

is not constant in `θ` and has a residual‑`Z2` unique minimizer set `θ* (mod π)`. Equivalently, the Shannon element‑order
reference mechanism cuts the continuous `O(2)` family on `pair_m` down to residual `Z2` (axis‑only uniqueness) on every `Z_n`.

This remains a **lane‑scoped arithmetic statement**: any strict export/promotion of a particular `n` (typed carrier, basis object,
and any physical scaffold identification) must still be handled separately.

## What `N514` does not claim

`N514` does not claim:

1. any strict physical promotion of `n≠12` into `QW-2190`,
2. any global discharge of `QW-2191`,
3. strict‑core selector closure / admissible `S_sel_int`,
4. any sign‑sensitive physical orientation datum,
5. ToE closure.
