# 325 — GLOBAL-B4-COMMUNICATION-HEIGHT-CERTIFICATE
## Exact mod-3 separator reduction plus a 6D constrained-stationary formulation

Date: 2026-09-28

Status: **PARTIAL — exact topological reduction + strong numerical global evidence; final interval exhaustion still open.**

Microscopic/statistical contract: the same rank-7 FIN potential and working gain

\[
g=5.145228719489142.
\]

The mapped inter-sector saddle barrier is

\[
B_4=0.6622191371274597.
\]

This report asks the remaining question left by 324: is the **global** continuous communication height between distinct mod-3 metastable sectors really `B4`, rather than merely the height in the already mapped saddle graph?

---

## 1. Exact sector separator

Let

\[
P_a(p)=\sum_{j\equiv a\pmod 3}p_j,
\qquad a=0,1,2.
\]

Equivalently these are the three projections encoded by the `k=4` Fourier order parameter. Define the three closed Voronoi sectors

\[
C_a=\{p: P_a\ge P_b\ \forall b\}.
\]

A continuous path from the interior of `C0` to a different mod-3 sector must cross one of the two walls

\[
S_{01}=\{P_0=P_1\ge P_2\},
\qquad
S_{02}=\{P_0=P_2\ge P_1\}.
\]

By `D12` symmetry it is sufficient to study

\[
\boxed{S_{01}:\ P_0=P_1\ge P_2.}
\]

This reduction is exact and topological; it does not assume the existing saddle census is complete.

---

## 2. The known d=4 saddle lies on the wall

The localized minimum has

\[
V_{\rm loc}=-0.8053194621423083.
\]

A representative `d=4` saddle has

\[
V_{d4}=-0.1431003250148470,
\]

hence

\[
V_{d4}-V_{\rm loc}
=\boxed{0.6622191371274613},
\]

numerically identical to the previously recorded `B4` to floating precision.

Its mod-3 masses satisfy

\[
P_0=P_1>P_2,
\]

so it lies in the interior of the relevant separator wall rather than at the triple-sector junction.

---

## 3. Full-simplex constrained multistart test

The potential was minimized directly in all twelve probability coordinates subject only to

\[
\sum_jp_j=1,
\qquad P_0=P_1,
\qquad P_0\ge P_2,
\qquad p_j\ge0.
\]

No reflection symmetry, saddle ansatz, or low-dimensional branch parametrization was imposed.

From **257 feasible independent starts**, all successful constrained minimizations gave

\[
\min_{S_{01}}V
=-0.1431003250148477,
\]

with difference from the `d=4` value below `7e-16`.

This is strong global numerical evidence, not a proof of globality.

---

## 4. d=4 is a strong local minimum on the separator

Although the same point is an index-one saddle in the full simplex, its Hessian restricted to the affine tangent of

\[
\sum p_j=1,
\qquad P_0=P_1
\]

is strictly positive.

The smallest restricted eigenvalue is

\[
\boxed{7.2593212142}.
\]

Thus the `d=4` state is a strongly isolated local minimum of the separator problem.

This is the correct mountain-pass geometry: unstable in the full state space, stable inside the separating wall.

---

## 5. Explicit route gives the B4 upper direction numerically

The straight segments

\[
p_{\rm loc}^{(0)}\to p_{d4}\to p_{\rm loc}^{(4)}
\]

were sampled at 20,001 points per segment.

For both segments the largest sampled potential is exactly the endpoint saddle value to floating precision; no overshoot was found.

Therefore a completely explicit path realizes a maximum numerically equal to `B4`.

This is currently a numerical path certificate, not yet outward-interval monotonicity certification.

---

## 6. Exact relative-entropy identity around the saddle

Because `p* = p_d4` is a full stationary point, for every simplex state `p`

\[
\boxed{
V_g(p)-V_g(p^*)
= D(p\Vert p^*)
-\frac g2(p-p^*)^TA_7(p-p^*).
}
\]

The replay verifies the identity to machine precision on the sampled constrained solutions.

Hence the global separator theorem is equivalent to proving

\[
D(p\Vert p^*)
\ge
\frac g2(p-p^*)^TA_7(p-p^*)
\]

for every `p` on the separator, with equality at the four symmetry-related `d=4` representatives.

This is a much sharper formulation than a generic barrier search.

---

## 7. Constrained KKT census

A second calculation did not minimize the energy directly. It solved the constrained stationarity equations from **3013 independent starts**.

On the representative half-wall it found **28 distinct stationary solutions** numerically.

The four lowest are the `C4` orbit of `d=4` and all have

- separator Morse index `0`;
- separator nullity `0`;
- barrier `B4`.

The next stationary-energy levels are approximately

\[
-0.06235983379,
\qquad
-0.02266475237,
\qquad
-1.3511\times10^{-5},
\qquad
0,
\]

all strictly above the `d=4` wall minimum candidate.

The numerical separator-Morse histogram of all found roots is

- index 0: `7` roots;
- index 1: `13` roots;
- index 2: `8` roots.

This finite census is **not** yet an interval proof that there is no 29th root.

---

## 8. Exact reduction of the constrained stationary problem to six dimensions

This is the most useful constructive result of 325.

Let

\[
c_j=
\begin{cases}
+1,&j\equiv0\pmod3,\\
-1,&j\equiv1\pmod3,\\
0,&j\equiv2\pmod3.
\end{cases}
\]

The normal direction `c` lies entirely in the retained `k=4` Fourier plane.
Choose six orthonormal mediator directions tangent to the wall and write their label-space field as

\[
h_j(y)=(By)_j,
\qquad y\in\mathbb R^6.
\]

Adding the normal multiplier gives field

\[
h_j(y)+z c_j.
\]

Define

\[
S_a(y)=\sum_{j\equiv a\pmod3}e^{h_j(y)}.
\]

The equality `P0=P1` determines the normal multiplier **exactly**:

\[
\boxed{
z(y)=\frac12\log\frac{S_1(y)}{S_0(y)}.
}
\]

Therefore the constrained KKT problem contains no free Lagrange multiplier. It is exactly the six-dimensional fixed point

\[
\boxed{
y=g B^Tp(y),}
\]

where

\[
p_j(y)\propto e^{h_j(y)+z(y)c_j}.
\]

Moreover every root lies in the explicit compact coordinate box recorded in `GLOBAL_B4_SEPARATOR_325.json`, because each component of `y` is `g` times an expectation of a finite label feature.

This is the correct target for the next interval/Krawczyk campaign.

---

## 9. Negative results — shortcuts that do not prove globality

Several tempting simpler routes were explicitly tested and rejected.

### 9.1 Independent-coordinate Bregman quadratic bound

Using the best coordinatewise global quadratic lower bounds for

\[
D(p\Vert p^*)
\]

still leaves the restricted comparison matrix with minimum eigenvalue about

\[
-5.225.
\]

So a diagonal KL bound is too lossy.

### 9.2 Reflection symmetrization by global convexity

The reflection exchanging sectors 0 and 1 does not make the potential globally convex in all antisymmetric directions. A numerical minimization of the corresponding Hessian lower model reaches approximately

\[
-7.11.
\]

Therefore symmetry of the minimizer cannot be asserted from a simple Jensen argument.

### 9.3 CRT + Pinsker

The exact CRT decomposition

\[
\mathbb Z_{12}\simeq\mathbb Z_3\times\mathbb Z_4
\]

splits entropy into mod-3 entropy, mod-4 entropy, and mutual information. But controlling the mixed `k=5` character only by Pinsker is far too weak: the resulting reduced lower bound descends to roughly

\[
-0.776,
\]

well below the required `-0.1431` wall energy.

### 9.4 Generic 12D Taylor/LP branch-and-bound

A first generic branch-and-bound over probability boxes also remains too loose. After about ten thousand processed boxes its best unresolved lower bound was still around `-0.762`.

The failure is useful: the next proof should work on the exact six-dimensional constrained fixed-point map, not return to a generic 12D energy cover.

---

## 10. Verdict

### Proved exactly in 325

1. the correct mod-3 sector separator;
2. every inter-sector path must cross one of the symmetry-equivalent separator walls;
3. the exact relative-entropy identity around `d=4`;
4. the exact six-dimensional constrained KKT reduction and analytic elimination of the normal multiplier.

### Strong numerical evidence

1. the global wall minimum is the `d=4` orbit;
2. its value is exactly the previously mapped `B4` within numerical precision;
3. no lower constrained stationary root was found in 3013-start root census;
4. the explicit localized-to-saddle-to-localized route has no sampled overshoot.

### Still open

An **interval-assisted exhaustion of the complete six-dimensional fixed-point box**. Until that is done, the statement

\[
\Gamma=B_4
\]

must remain unproved.

Therefore 325 is recorded as

\[
\boxed{\textbf{PARTIAL — major reduction, final global exclusion open.}}
\]

Combined with 324, a successful interval exhaustion in the next task would immediately imply the proof-grade capacity exponent

\[
-\frac1N\log\operatorname{cap}_N\to B_4.
\]
