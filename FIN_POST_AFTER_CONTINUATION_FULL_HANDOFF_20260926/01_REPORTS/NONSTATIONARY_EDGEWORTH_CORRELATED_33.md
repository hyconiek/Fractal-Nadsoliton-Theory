# NONSTATIONARY-EDGEWORTH-CORRELATED-33
## Correlated hidden preparation is quadratically identifiable

Date: 2026-09-26

Repository baseline inspected:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
**PROOF-GRADE CONDITIONAL inside the declared Gaussian-projected finite-N
heat-bath/ME7 lane.**

This extends `NONSTATIONARY-EDGEWORTH-INITIAL-SLIP-32` and uses the exact
four-channel hidden multiplication structure already established by
`HIDDEN-MULTIPLICATION-RANK-14` and `EDGEWORTH-STATIONARY-15`.

It is NOT a derivation of physical time, a laboratory detector model, or a
microscopic ontology.

---

## 1. Correlated preparation

Let the leading hidden residual have conditional mean

\[
E[z_0\mid x_0=x]=m+Mx,
\]

with

- \(m\in\mathbb R^4\),
- \(M\in\mathbb R^{4\times 7}\).

At the first nonstationary order \(\varepsilon=N^{-1/2}\), the hidden part of
the visible backward generator is

\[
L^{(1)}_{\rm hid}f(x)
=
\frac12\sum_{a=1}^4 \bigl(m_a+(Mx)_a\bigr)\,O_a f(x),
\]

where

\[
O_a=T_a:\nabla^2,\qquad
T_a=X^T\operatorname{diag}(Y_a)X.
\]

Therefore correlated hidden preparation produces **state-dependent diffusion**.

It does not create an ordinary first-derivative drift term at this order.

---

## 2. Linear observables are exactly blind

For every linear observable

\[
\ell_v(x)=v^T x,
\]

\[
\nabla^2\ell_v=0,
\]

hence

\[
L^{(1)}_{\rm hid}\ell_v=0.
\]

So all seven coordinate means are blind to this \(N^{-1/2}\) correlated-hidden
correction at generator level.

This gives a clean diagnostic distinction:

- a genuine visible drift \(b(x)\cdot\nabla\) acts on linear observables;
- the correlated-hidden correction above does not.

Thus the two mechanisms are not operator-equivalent once both linear and
quadratic probes are available.

---

## 3. Exact linear independence of the four hidden diffusion tensors

Use the Euclidean-orthonormal real Fourier basis
\(e_{3c},e_{3s},e_{4c},e_{4s},e_{5c},e_{5s},e_6\) for the retained space and
\(e_{1c},e_{1s},e_{2c},e_{2s}\) for the hidden space.

The strict feature columns satisfy

\[
X_{kc}=\sqrt{\lambda_k}\,e_{kc},\qquad
X_{ks}=\sqrt{\lambda_k}\,e_{ks},\qquad
X_6=\sqrt{\lambda_6}\,e_6.
\]

The exact product identities inherited from
`HIDDEN-MULTIPLICATION-RANK-14` give four isolating entries:

\[
(T_{1c})_{3c,4c}
=
\sqrt{\lambda_3\lambda_4}\frac{\sqrt6}{12},
\]

\[
(T_{1s})_{3c,4s}
=
\sqrt{\lambda_3\lambda_4}\frac{\sqrt6}{12},
\]

\[
(T_{2c})_{4c,6}
=
\sqrt{\lambda_4\lambda_6}\frac{\sqrt3}{6},
\]

\[
(T_{2s})_{4s,6}
=
-\sqrt{\lambda_4\lambda_6}\frac{\sqrt3}{6}.
\]

Since the strict retained eigenvalues are positive, all four coefficients are
nonzero. These four matrix entries isolate the four hidden directions.

Therefore

\[
T_1,T_2,T_3,T_4
\]

are linearly independent.

This is an exact algebraic statement; it does not rely on numerical rank
thresholds.

---

## 4. Gram matrix and conditioning for the frozen strict operator

Define

\[
G_{ab}=\operatorname{tr}(T_aT_b).
\]

For the frozen strict kernel,

\[
G \approx
\operatorname{diag}(
2.458991086258,\,
2.458991086258,\,
2.050348028205,\,
2.050348028205).
\]

Its eigenvalues are therefore

\[
2.050348028205,\,
2.050348028205,\,
2.458991086258,\,
2.458991086258,
\]

and

\[
\kappa_2(G)\approx1.19930424.
\]

So the four channels are not only identifiable; the quadratic detector system
is numerically well conditioned in the frozen strict model.

---

## 5. Dual quadratic detector theorem

Define

\[
Q_b=\sum_{c=1}^4 (G^{-1})_{bc}T_c
\]

and the four quadratic observables

\[
f_b(x)=x^TQ_bx.
\]

Because

\[
\nabla^2 f_b=2Q_b,
\]

we obtain

\[
O_a f_b
=
2\,\operatorname{tr}(T_aQ_b)
=
2\delta_{ab}.
\]

Substituting into the hidden generator gives the exact identity

\[
\boxed{
L^{(1)}_{\rm hid} f_b(x)
=
m_b+(Mx)_b
}
\]

for \(b=1,\dots,4\).

Thus four fixed quadratic probes directly read out the complete conditional
hidden mean.

No cubic or quartic probe is needed for this first-order nonstationary
identifiability problem.

---

## 6. Minimal polynomial degree

Degree 1 is insufficient because every linear observable has zero Hessian.

Degree 2 is sufficient by the dual detector construction.

Therefore the minimal polynomial degree required to detect the
\(N^{-1/2}\) correlated-hidden preparation is exactly

\[
\boxed{2}.
\]

This is stronger than merely saying that hidden modes are quadratically
reachable.

---

## 7. Reconstruction of the full \(m,M\)

For each detector \(b\),

\[
R_b(x):=L^{(1)}_{\rm hid} f_b(x)
=
m_b+M_bx.
\]

Hence the \(b\)-th row of \(M\) and the intercept \(m_b\) form one affine
function on the seven-dimensional visible state space.

Using preparations

\[
x^{(0)}=0,\quad
x^{(j)}=\delta e_j,\ j=1,\dots,7,
\]

one gets

\[
m_b=R_b(0),
\]

\[
M_{bj}
=
\frac{R_b(\delta e_j)-R_b(0)}{\delta}.
\]

Therefore, in the ideal declared model:

- four quadratic observables,
- evaluated at eight affinely independent visible preparations,

are sufficient to recover all

\[
4+4\times7=32
\]

parameters of the conditional hidden mean.

The same eight preparations can be used for all four detectors, so this is
eight preparation settings with four measured generator responses each.

This is a mathematical identifiability statement, not a sample-size or
laboratory-feasibility result.

---

## 8. Drift-versus-hidden-memory separation

Suppose an alternative visible correction contains ordinary drift

\[
D f=b(x)\cdot\nabla f.
\]

For coordinate probes \(f=x_i\),

\[
D x_i=b_i(x),
\qquad
L^{(1)}_{\rm hid}x_i=0.
\]

Hence if a correction changes first moments, it cannot be explained solely by
the correlated-hidden diffusion term above.

Conversely, a correction that leaves all linear probes unchanged but modifies
the four dual quadratic detectors exactly as
\(m+Mx\) has the signature of the declared correlated-hidden mechanism.

Thus there is no exact drift/diffusion alias at generator level once linear
and dual quadratic probes are both admitted.

---

## 9. Effect on higher-degree observables

For a polynomial of degree \(d\),

\[
O_a=T_a:\nabla^2
\]

lowers degree by two.

Multiplication by \(m_a\) leaves the result at degree \(d-2\);
multiplication by \((Mx)_a\) raises it once, yielding degree \(d-1\).

Therefore

\[
L^{(1)}_{\rm hid}:\mathcal P_d\to\mathcal P_{d-1}.
\]

In particular:

- degree 1 -> 0,
- degree 2 -> affine,
- degree 3 -> quadratic,
- degree 4 -> cubic.

So cubic and quartic observables carry consistency checks and nonlinear
structure, but add no new identifiability requirement for \(m,M\).

---

## 10. Physical interpretation boundary

What this supports:

A hidden relational state can leave a measurable signature in the evolution
of second-order visible structure even while all first-order visible means are
unchanged.

That is a concrete mechanism by which "memory of relations" can survive
coarse-graining.

What it does NOT support:

- that these four hidden directions are physical extra dimensions;
- that the heat-bath clock is physical time;
- that the hidden state is a physical field;
- that FIN has derived inertia, quantum mechanics, gravity or matter.

---

## 11. Consequence for the relational-time programme

The earlier clock-free escape-profile result used the exponent

\[
\gamma=
\frac{d\log\|r-r_0\|}{d\log S}.
\]

The present theorem says that nonstationary hidden preparation is not an
unobservable nuisance: in principle it can be independently estimated from
quadratic responses.

This creates a new route:

1. estimate/reject correlated hidden preparation using the four dual
   quadratic detectors;
2. condition the clock-free dynamic classifier on the resulting hidden-state
   estimate;
3. only then interpret the escape-profile exponent.

This is substantially stronger than assuming stationarity by fiat.

---

## 12. Next research atom

### RELATIONAL-CLOCK-HIDDEN-CORRECTION-34

Given an operational escape/readout vector \(p(t,x)\), augment its short-time
expansion with the reconstructed hidden preparation \(m+Mx\).

Derive whether a calibrated subtraction based on the four quadratic detectors
restores the intrinsic projective exponent \(\gamma\).

Acceptance:

- a theorem proving first-order de-aliasing of \(\gamma\) under the declared
  finite-N heat-bath preparation class; or
- an explicit obstruction showing that the hidden correction enters the
  escape record through directions not determined by \(m,M\).

This is the next direct bridge between FIN hidden-memory mathematics and an
observable notion of emergent relational time.
