# 326 — CONSTRAINED-SEPARATOR-INTERVAL-EXHAUSTION
## The mod-3 communication height is the d=4 barrier B4 in the declared FIN potential

Date: 2026-09-28

Scientific status: **PASS — computer-assisted global separator certificate in the declared finite FIN model.**

Arithmetic note: the two global covers use conservative long-double enclosure formulas with explicit downward inflation `1e-8` on every Taylor lower bound and accept an energy exclusion only with margin `>1e-6`. The local root and path certificates use high-precision `mpmath.iv` interval arithmetic. A future directed-rounding MPFI/Arb replay would be a valuable independent formal-audit upgrade, but no state-space box remains unresolved in the present certificate.

This report closes the only remaining geometric blocker left by reports 324–325.

At the working gain

\[
g=5.145228719489142,
\]

the previously mapped inter-sector barrier was

\[
B_4=0.6622191371274597.
\]

Report 324 proved that if the **global** continuous communication height between the relevant metastable sectors equals `B4`, then the declared leave-one-out Gibbs capacity has exponential rate `B4`. Report 325 reduced the missing global statement to a separator problem but stopped before interval exhaustion. Report 326 completes that exhaustion.

---

## 1. Exact separator theorem inherited from 325

Define the three mod-3 masses

\[
P_a(p)=\sum_{j\equiv a\pmod3}p_j,
\qquad a=0,1,2.
\]

The three closed metastable Voronoi sectors are

\[
C_a=\{p:P_a\ge P_b\;\forall b\}.
\]

Every continuous path from the interior of `C0` to a distinct mod-3 sector must cross one of

\[
P_0=P_1\ge P_2,
\qquad
P_0=P_2\ge P_1.
\]

By `D12` symmetry it is sufficient to certify the representative half-wall

\[
\boxed{P_0=P_1\ge P_2.}
\]

The known `d=4` saddle lies in its interior and has

\[
V_{d4}-V_{\rm loc}=B_4.
\]

Thus proving that this saddle is the global minimum of the representative separator gives the lower mountain-pass inequality `Gamma >= B4`. A certified path through the same saddle gives `Gamma <= B4`.

---

## 2. Exact six-dimensional equality-wall dual

Let

\[
c_j=\begin{cases}
+1,&j\equiv0\pmod3,\\
-1,&j\equiv1\pmod3,\\
0,&j\equiv2\pmod3.
\end{cases}
\]

The normal direction lies entirely in the retained `k=4` Fourier plane. Choose six orthonormal mediator directions tangent to `P0=P1`, with label-space field

\[
h_j(y)=(By)_j,
\qquad y\in\mathbb R^6.
\]

Write

\[
S_a(y)=\sum_{j\equiv a\pmod3}e^{h_j(y)}.
\]

The normal multiplier enforcing `P0=P1` is eliminated **exactly**:

\[
\boxed{
z(y)=\frac12\log\frac{S_1(y)}{S_0(y)}.
}
\]

The equality-wall Gibbs state is therefore

\[
p_j(y)=\frac{\exp(h_j(y)+z(y)c_j)}
{2\sqrt{S_0S_1}+S_2}.
\]

The corresponding smooth dual potential is

\[
\boxed{
\Phi_{\rm sep}(y)
=\frac{\|y\|^2}{2g}
-\log\!\left(\frac{2\sqrt{S_0(y)S_1(y)}+S_2(y)}{12}\right).
}
\]

Its gradient is

\[
\nabla\Phi_{\rm sep}(y)
=\frac yg-B^Tp(y),
\]

so every interior constrained extremum satisfies

\[
\boxed{y=gB^Tp(y).}
\]

Any such stationary point lies in the explicit compact box

```text
[-2.9417983842,  2.9417983842]
[-2.9417983842,  2.9417983842]
[-3.1152854030,  1.5576427015]
[-3.1846473912,  3.1846473912]
[-3.1846473912,  3.1846473912]
[-2.2731305848,  2.2731305848].
```

---

## 3. Exact internal C4 quotient

Translation by three labels preserves `P0=P1>=P2` and leaves the potential invariant. On the `k=3` Fourier pair it acts as a quarter-turn.

Hence every internal `C4` orbit intersects the closed wedge

\[
\boxed{y_0\ge |y_1|.}
\]

The global cover can therefore be run on this fundamental wedge without losing any candidate minimum. The four symmetry-related `d=4` minima reduce to a single representative local box.

---

## 4. Root-preserving interval-style contractor

For each six-dimensional box `[l,u]` the certificate computes rigorous enclosure formulas for the equality-wall probabilities.

Linear field bounds are exact box extrema. Class sums bound

\[
S_a=\sum_{j\equiv a}e^{h_j}.
\]

With

\[
q=\sqrt{S_0S_1},
\]

the common class mass satisfies

\[
P_0=P_1=\frac{q}{2q+S_2}.
\]

For each label the within-class softmax is bounded separately, producing `p_j` intervals. Linear programming over these component bounds and normalization gives enclosures of every component of `g B^T p`.

A box is root-free whenever one coordinate interval of `y` is disjoint from the corresponding interval of `gB^Tp`.

### Important boundary correction

The representative half-wall has the extra inequality

\[
P_0=P_1\ge P_2.
\]

Since `P0=P1=m` and `P2=1-2m`, this is exactly

\[
\boxed{m\ge1/3.}
\]

Boxes with `m_U<1/3` are outside the half-wall. Boxes that can touch `m=1/3` are **not** excluded by the interior fixed-point equation, because a boundary minimizer can carry an additional KKT multiplier. They are instead handled by the energy bound and by the separate five-dimensional triple-junction certificate in Section 7.

---

## 5. Certified energy lower bound on every box

The separator Hessian has the exact structure

\[
\nabla^2\Phi_{\rm sep}
=\frac1gI-K,
\]

where `K` is the covariance of the tangent features after conditioning out the normal equality direction.

Because the conditioning correction is positive semidefinite,

\[
K\preceq B^T(\operatorname{diag}p-pp^T)B.
\]

For a box, let `p_i in [p_i^L,p_i^U]`. Exact small linear programs give

- an upper bound on `E||B_J||^2`;
- coordinate intervals for `mu=E[B_J]`.

Therefore

\[
\operatorname{tr}K
\le
E\|B_J\|^2-\|E B_J\|^2
\le T_B,
\]

and throughout the box

\[
\nabla^2\Phi_{\rm sep}\succeq m_B I,
\qquad
m_B=\frac1g-T_B.
\]

At box center `c`, with half-width vector `w`, Taylor's theorem gives

\[
\Phi(c+d)
\ge
\Phi(c)+\nabla\Phi(c)\cdot d
+\frac{m_B}{2}\|d\|^2.
\]

The right-hand side is minimized analytically coordinate-by-coordinate over `|d_i|<=w_i`.

For arithmetic robustness the implementation then subtracts an additional

\[
10^{-8}
\]

from every computed lower bound, and a box is energy-excluded only if the remaining margin above the `d=4` comparison value is larger than

\[
\boxed{10^{-6}}.
\]

---

## 6. Complete interior-half-wall cover

The boundary-safe global replay processed

\[
\boxed{1,548,513\text{ boxes}}.
\]

Results:

```text
boxes discarded:       772,626
boxes inside d4 local cube: 1,631
unresolved boxes:             0
maximum depth:               59
minimum accepted discard margin: 1.0029449815e-6
```

Every box outside the certified local `d=4` cube is therefore either

1. outside the representative half-wall / `C4` wedge;
2. incompatible with the stationary fixed-point equation when interior;
3. or has a certified potential lower bound strictly above the `d=4` level.

No second interior region survives.

---

## 7. Triple-junction boundary certificate

The remaining logical boundary is

\[
P_0=P_1=P_2=1/3.
\]

On this boundary the whole `k=4` Fourier plane is normal and disappears from the tangent dynamics. The problem reduces exactly to **five dimensions**, with tangent field built only from

```text
k=3 cosine/sine,
k=5 cosine/sine,
k=6.
```

For `h=By`, each mod-3 class is normalized independently:

\[
p_j(y)=\frac{e^{h_j(y)}}{3S_a(y)},
\qquad j\equiv a\pmod3.
\]

The exact five-dimensional dual is

\[
\boxed{
\Phi_{\rm tri}(y)
=
\frac{\|y\|^2}{2g}
-\frac13\sum_{a=0}^2\log S_a(y)
+\log4.
}
\]

A separate `C4`-quotiented root/energy cover processed

\[
\boxed{54,203\text{ boxes}}
\]

with

```text
feasibility exclusions:  9,306
energy exclusions:      17,796
unresolved boxes:            0
maximum depth:              20
minimum energy-discard margin: 3.5946418672e-6
```

No triple-junction stationary point can have energy at or below the `d=4` saddle level. Numerically its actual lowest branch is the already-known level

\[
-0.0623598338\ldots,
\]

more than `0.0807` above `V_d4`, but the proof only needs the certified exclusion relative to `V_d4`.

---

## 8. Local d=4 isolation and convexity

High-precision `mpmath.iv` Krawczyk arithmetic isolates the representative `d=4` root at

```text
y = (
  2.57030485631592353839,
 -2.6101854e-16,
  1.29925166569059554588,
  0.68179265329213179194,
 -1.18089951572916437736,
  1.97598637879081091047
)
```

inside a radius

\[
10^{-18}
\]

box, with strict Krawczyk inclusion in all six coordinates.

A separate high-precision probability-box / LP Hessian enclosure on the full coordinate cube

\[
|y-y^*|_\infty\le0.02
\]

gives

\[
\boxed{
\lambda_{\min}(\nabla^2\Phi_{\rm sep})
\ge0.0551088593996>0.
}
\]

Thus every box accepted into that local cube contains no lower competing point: the isolated root is the unique local minimum there.

By exact `C4` symmetry this produces exactly the four symmetry-related separator minima found numerically in 325.

---

## 9. Explicit path gives the upper mountain-pass bound

Let `p_loc^(0)` be a localized minimum, `p_d4` the adjacent `d=4` saddle, and `p_loc^(4)` the translated neighboring minimum.

Consider the two straight segments

\[
p(t)=(1-t)p_{\rm loc}^{(0)}+tp_{d4},
\]

and

\[
q(t)=(1-t)p_{d4}+tp_{\rm loc}^{(4)}.
\]

High-precision interval subdivision certifies that

\[
\frac{d}{dt}V(p(t))>0
\quad(0<t<1),
\]

and

\[
\frac{d}{dt}V(q(t))<0
\quad(0<t<1).
\]

On the interval-arithmetic middle region the worst signed derivative bounds are

\[
\boxed{+0.00239277474965}
\]

and

\[
\boxed{-0.00239277474965}.
\]

Endpoint neighborhoods are handled by second-derivative signs. For the first segment:

```text
near t=0:  V'' in [43.2971, 54.7616] > 0
near t=1:  V'' in [-1.61593, -1.61313] < 0
```

and the signs reverse appropriately on the second segment.

Therefore the explicit continuous route has maximum **exactly at the `d=4` saddle**.

Hence

\[
\Gamma\le V_{d4}-V_{\rm loc}=B_4.
\]

---

## 10. Global communication-height theorem

Sections 1–8 prove that any path changing mod-3 sector must cross a separator whose global minimum consists only of the four `d=4` images. Therefore

\[
\Gamma\ge B_4.
\]

Section 9 supplies a path attaining that level, so

\[
\Gamma\le B_4.
\]

Thus, in the declared continuous FIN potential at the working gain,

\[
\boxed{
\Gamma=B_4.
}
\]

The high-precision numerical value from the same stationary points is

\[
\boxed{
B_4
=0.6622191371274619267259530929\ldots
}
\]

(the tiny difference from the earlier stored `0.6622191371274597` is ordinary floating representation of the same barrier).

This result no longer depends on completeness of the previously mapped `d3/d4/d5` saddle graph.

---

## 11. Combination with report 324

Report 324 proved for the declared leave-one-out finite-N Gibbs chain that if the continuous communication height is `Gamma`, then

\[
\lim_{N\to\infty}
-\frac1N\log \operatorname{cap}_N(A,B)
=\Gamma,
\]

because

- the number of count states is only polynomial in `N`;
- stationary type weights have exponent `Phi` with only `O(log N)` corrections;
- allowed single-step rates are subexponential;
- Thomson path bounds and Dirichlet cut bounds squeeze the exponent to the communication height.

Combining 324 and 326 gives

\[
\boxed{
\lim_{N\to\infty}
-\frac1N\log \operatorname{cap}_N(A,B)
=B_4.
}
\]

Under the same declared source/target and order-one metastable-valley-mass conditions used in 324, the exit-rate exponent is therefore also `B4`.

This is the first proof-grade closure of the **exponential metastable scale** in this FIN lane.

---

## 12. What is still not proved

326 does **not** prove an Eyring–Kramers prefactor.

Still open:

1. the polynomial / determinant prefactor of the capacity and exit time;
2. a controlled finite-N remainder sharp enough to predict moderate N from the asymptotic formula alone;
3. the full six-observable analytic preparation coefficients left open by 323;
4. any dimensional physical clock, temperature, laboratory realization, Standard Model / GR identification, or ToE claim.

The result is a theorem about the declared dimensionless finite FIN Gibbs potential and its specified microscopic heat-bath chain.

---

## 13. Reproducibility

Principal artifacts:

- `GLOBAL_B4_CERTIFICATE_326.json`
- `interval_exhaust_326_boundarysafe.cpp`
- `interval_exhaust_326_safe3.cpp`
- `triple_boundary_exhaust_326.cpp`
- `verify_local_and_path_326.py`
- `local_hessian_custom_326.py`
- global and triple-junction replay stdout/stderr files.

The final global cover is the boundary-safe version; earlier intermediate executables are retained only for audit provenance.

## Final verdict

\[
\boxed{\textbf{326 PASS}}
\]

with the arithmetic qualification stated at the beginning of the report.
