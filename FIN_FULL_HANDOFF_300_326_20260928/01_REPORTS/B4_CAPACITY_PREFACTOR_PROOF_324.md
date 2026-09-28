# FIN — 324 B4-CAPACITY/PREFACTOR

Date: 2026-09-27

Status: **CONDITIONAL THEOREM for the capacity exponent; prefactor still open.**

This task continues the parallel metastability proof lane after 318 and 322. No new N was opened.

## 1. Main result

For the exact finite-N leave-one-out Gibbs heat-bath chain used in the current FIN lane, the following implication can now be proved from the microscopic form of the model:

> If the **global continuous communication height** between the declared Z3 metastable sectors is
>
> \[
> \Gamma=B_4=0.6622191371274597,
> \]
>
> then
>
> \[
> \boxed{
> \lim_{N\to\infty}-\frac1N\log\operatorname{cap}_N(A,B)=B_4.
> }
> \]

If the valley mass stays order one, as required by the symmetry/metastable-sector construction, then the corresponding exit-rate exponent is also `B4`.

This removes an important earlier ambiguity: **potential theory does not introduce a different exponential barrier once the global communication height is known.**

The still-missing theorem is the global statement `Gamma=B4` itself.

## 2. Uniform stationary large-deviation bound from type classes

For a count state `n`, `p=n/N`, the exact invariant weight is

\[
\pi_N(n)\propto
\frac{N!}{\prod_i n_i!}
\exp\left[\frac{g}{2N}n^T A n\right].
\]

The elementary method-of-types inequality gives uniformly, including simplex boundaries,

\[
(N+1)^{-12}e^{NH(p)}
\le
\frac{N!}{\prod_i n_i!}
\le
e^{NH(p)}.
\]

Hence the finite-state weight has the uniform exponential form

\[
\pi_N(n)
=
\exp\{-N[\Phi(p)-\min\Phi]+O(\log N)\},
\]

with

\[
\Phi(p)=\sum_i p_i\log p_i-\frac g2p^TAp.
\]

Because the number of count states is

\[
\binom{N+11}{11}=N^{11+o(1)},
\]

normalization itself contributes only `O(log N)` to `-log pi`.

Thus the large-deviation exponent needed by capacity theory is already an exact consequence of the finite-N model; no Gaussian saddle approximation is needed for the exponent.

## 3. Transition rates are subexponential

For an allowed one-particle move `i -> j`,

\[
q_N(n,n-e_i+e_j)
=n_i\,\mathrm{softmax}_j\left[\frac gN A(n-e_i)\right].
\]

Since `p` lies in a compact simplex and `A` is fixed, the field-coordinate range is uniformly bounded.

For the present kernel the explicit crude bound is

\[
\operatorname{range}(h)\le10.19801472,
\]

so every softmax component is at least

\[
\boxed{3.10368\times10^{-6}}.
\]

Therefore on every allowed edge with `n_i>=1`,

\[
c_*\le q_N\le N,
\]

with `c_*>0` independent of `N`.

Rates can modify only polynomial/subexponential prefactors, not the exponential barrier.

## 4. Lower bound on capacity from a path

For a reversible chain, edge conductance is

\[
c_N(x,y)=\pi_N(x)q_N(x,y).
\]

Assume the global communication height is `Gamma`. For any `epsilon>0`, choose a continuous path between the two wells with

\[
\max\Phi\le\Phi_{min}+\Gamma+\epsilon.
\]

The count lattice has mesh `1/N`, so a lattice path can approximate it with polynomially many nearest-neighbor steps and only `o(1)` change in the maximal potential.

Sending unit flow along this path and using Thomson's principle gives

\[
\operatorname{cap}_N(A,B)
\ge
N^{-C_1}
\exp[-N(\Gamma+\epsilon+o(1))].
\]

Thus

\[
\limsup_N -\frac1N\log\operatorname{cap}_N(A,B)
\le\Gamma.
\]

## 5. Upper bound from a sublevel cut

Fix `epsilon>0`. Let `C_epsilon` be the connected component of the source well in

\[
\{p:\Phi(p)<\Phi_{min}+\Gamma-\epsilon\}.
\]

By definition of global communication height, the target well is not in that component.

Use as Dirichlet test function the indicator of the corresponding lattice component, smoothed only across its edge boundary if desired. Every crossing edge must reach potential at least

\[
\Phi_{min}+\Gamma-\epsilon-o(1).
\]

The number of states/edges is polynomial and each rate is at most `N`, hence

\[
\operatorname{cap}_N(A,B)
\le
N^{C_2}
\exp[-N(\Gamma-\epsilon-o(1))].
\]

Therefore

\[
\liminf_N -\frac1N\log\operatorname{cap}_N(A,B)
\ge\Gamma.
\]

Letting `epsilon -> 0` proves the theorem.

## 6. Consequence for B4

Conditional on the still-missing global statement

\[
\Gamma=B_4,
\]

we have

\[
\boxed{
-\frac1N\log\operatorname{cap}_N\to B_4.
}
\]

Thus the proof boundary is sharper than in reports 228/233/260:

### No longer an independent blocker for the exponent

- whether the exact finite-N invariant measure has the right LDP exponent;
- whether heat-bath rates change the exponent;
- whether capacity could acquire a different exponential rate than communication height.

These are controlled by the argument above.

### Still the central blocker

\[
\boxed{
\text{Global communication-height proof: no path below }B_4.
}
\]

The mapped `d3/d4` saddle graph does not by itself exclude an unseen route elsewhere in the full simplex.

## 7. Prefactor remains open

324 does **not** prove an Eyring-Kramers prefactor.

A prefactor theorem would require at least:

1. nondegenerate local saddle geometry at the true global gate;
2. mobility/conductance quadratic form of the leave-one-out chain near that saddle;
3. a tube theorem showing all reactive flow is captured by the named gate(s);
4. control of boundary/simplex effects;
5. a uniform remainder.

The exact finite-N rates remain compatible with `B4` plus a polynomial prefactor, but the fitted power is not stable enough to promote:

```
fit window N>=3: alpha ~0.619
N>=6:            alpha ~0.617
N>=7:            alpha ~0.563
N>=8:            alpha ~0.504
```

This drift is precisely why a fitted prefactor is not a theorem.

## 8. Relation to current finite-N evidence

The exact deep-core local exponents continue:

```
6->7   0.54030
7->8   0.56339
8->9   0.58583
9->10  0.60517
10->11 0.62090
11->12 0.63334
```

against

\[
B_4=0.6622191371.
\]

This trend is consistent with the theorem's conditional conclusion but is not used to prove it.

## Verdict

**324 PASS** for a conditional proof of the exponential capacity law:

\[
\Gamma=B_4\quad\Longrightarrow\quad
-\frac1N\log\operatorname{cap}_N\to B_4.
\]

**324 OPEN** for:

- proving globally that `Gamma=B4`;
- the Eyring-Kramers/polynomial prefactor;
- explicit finite-N remainder bounds.

The research target is therefore narrower and clearer: the next proof campaign should attack global communication-path exclusion, not fit another finite-N exponent.
