# RELATIONAL-CLOCK-HIDDEN-BRIDGE-34
## Typed gap, minimal endpoint bridge, and memory-induced exponent aliasing

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`,
visible HEAD `fe14a6f4e436815635df54429102f23c22862296`.

Status:
- **PROOF-GRADE typed-gap / nonidentifiability result** for the current repository objects;
- **EXACT conditional classification** of endpoint-local D12-equivariant linear edge bridges;
- **CONDITIONAL candidate theorem** for one minimal fractional bridge;
- **NUMERICAL replay** of the predicted clock-free exponents under that candidate.

This does not promote the candidate bridge to a physical FIN law.

---

## 1. The present repository contains two distinct dynamic objects

The finite-N Edgeworth lane uses visible/hidden fluctuation variables

\[
(x,z)\in\mathbb R^7\times\mathbb R^4
\]

and produces the nonstationary hidden-preparation correction

\[
L^{(1)}_{\rm hid}
=
\frac12\sum_a (m_a+(Mx)_a)\,T_a:\nabla_x^2.
\]

The operational channel lane instead uses fixed objects \(W,A\) and records

- heat: \(e^{-tA}\),
- unitary-like: \(e^{-itA}\),
- zero-velocity wave-like: \(\cos(t\sqrt A)\).

The current channel formulas contain no \(x,z,m,M\).

Therefore the repository currently supplies no typed map

\[
(x,z)\quad\longrightarrow\quad
\delta A,\ \delta W,\ \delta{\rm preparation},\ \delta{\rm detector}
\]

connecting the finite-N hidden state to the channel record.

---

## 2. Exact typed-gap consequence

Read literally, with the current declared channel definitions held fixed,

\[
\frac{\partial p_{\rm channel}}{\partial m}
=
\frac{\partial p_{\rm channel}}{\partial M}
=
0.
\]

So the newly reconstructed hidden preparation cannot currently be "subtracted"
from the escape-profile record, because no nonzero coupling has been defined.

Conversely, once an additional coupling is admitted, the same reconstructed
\((m,M)\) can be sent into many inequivalent channel perturbations.

Hence \((m,M)\) alone does not determine a correction to the relational-clock
curve.

This is a missing-object theorem, not a numerical limitation.

---

## 3. Smallest endpoint-local bridge class

Let

\[
u\in H_4
\]

be the hidden scalar field on the twelve labels.

Assume an infinitesimal edge-weight perturbation has all of the following
properties:

1. linear in the endpoint field values;
2. local to the two endpoints of the edge;
3. symmetric under exchanging the endpoints;
4. D12-equivariant;
5. depends on the base edge only through its unoriented cyclic shell
   \(d=1,\dots,6\).

For an edge \(i\leftrightarrow j\) in shell \(d\), the most general rule is

\[
\delta w_{ij}=c_d\,(u_i+u_j).
\]

Thus the complete endpoint-local class is six-dimensional:

\[
(c_1,\ldots,c_6).
\]

This classification follows directly from linearity, endpoint exchange
symmetry, and D12 orbit classification of undirected edges.

No desired physical answer enters the derivation.

---

## 4. One-parameter fractional bridge candidate

If one further requires the **same fractional response on every shell**,

\[
\frac{\delta w_{ij}}{w_{ij}}
=
\frac{\alpha}{2}(u_i+u_j),
\]

then

\[
c_d=\frac{\alpha}{2}w_d.
\]

The bridge is unique up to the overall scalar \(\alpha\).

Set \(\alpha=1\) as a dimensionless normalization convention for the following
candidate calculations.

For

\[
D_u=\operatorname{diag}(u),
\]

the induced symmetric Laplacian perturbation is

\[
\boxed{
J(u)=
\frac12(D_uA+AD_u)
-\frac12\operatorname{diag}(Au)
}
\]

and satisfies exactly

\[
J(u)^T=J(u),\qquad
J(u)\mathbf 1=0.
\]

It is D12-equivariant:

\[
J(Pu)=P\,J(u)\,P^T.
\]

For sufficiently small coupling it also preserves positivity of the underlying
edge weights.

This makes it a mathematically clean candidate bridge, but not a sourced FIN
law.

---

## 5. The candidate bridge sees all four hidden modes

For the frozen strict operator and the orthonormal hidden Fourier basis
\(Y=(k=1c,1s,2c,2s)\), define

\[
J_a=J(Y_a).
\]

Their Frobenius Gram matrix has eigenvalues

\[
1.111876834758,\,
1.111876834758,\,
2.123766145338,\,
2.123766145338.
\]

So the map

\[
u\in H_4\mapsto J(u)
\]

is injective.

Its Gram condition number is about

\[
1.91007.
\]

Hence this candidate does not erase any of the four hidden directions.

---

## 6. Hidden field is visible already in the leading destination profile

Fix source node \(i\).

### Heat

The leading conditional destination profile is

\[
r_H(j)=\frac{w_{ij}}{\sum_{k\ne i}w_{ik}}.
\]

Under the fractional bridge, the Jacobian from the four hidden coordinates to
the eleven-dimensional normalized profile has full rank four.

For the strict kernel its singular values are approximately

\[
0.0913398,\ 0.0521640,\ 0.0194286,\ 0.0121173.
\]

### Unitary / wave leading profile

Both use

\[
r_{UW}(j)=
\frac{A_{ji}^2}{\sum_{k\ne i}A_{ki}^2}
\]

at leading order.

Its hidden-profile Jacobian also has rank four, with singular values

\[
0.241696,\ 0.0622075,\ 0.0148105,\ 0.00441345.
\]

Thus, inside this candidate bridge, one full source-resolved leading profile
already contains local information about all four hidden directions.

This is bridge-dependent identifiability; it is not available before the
bridge is chosen.

---

## 7. Static hidden preparation does not change the intrinsic exponent

Let

\[
A_h=A+\varepsilon J(h)
\]

with \(h\) fixed during one channel run.

If the leading profile \(r_0\) is defined self-consistently from \(A_h\), then
the previously derived projective exponents remain

\[
\gamma_H=1,\qquad
\gamma_U=1,\qquad
\gamma_W=\frac12
\]

for the declared heat, unitary and wave categories.

A static hidden field therefore acts like operator calibration.

However, comparing to the *unperturbed* \(r_0(A)\) would create a constant
\(O(\varepsilon)\) profile offset and formally drive the measured asymptotic
slope to zero.

So re-centering by the actual leading profile is essential.

---

## 8. Relaxing hidden memory changes the exponent

Now use the OU-like preparation decay already present in the nonstationary
Edgeworth lane:

\[
h(t)=e^{-t}h_0.
\]

Under the candidate bridge,

\[
A(t)=A+\varepsilon e^{-t}J(h_0).
\]

This creates a qualitatively new effect.

### 8.1 Heat

The integrated jump rate contains

\[
\int_0^t e^{-s}ds
=
t-\frac{t^2}{2}+O(t^3).
\]

After re-centering by the \(t=0\) perturbed leading profile, the first profile
change is \(O(t)\).

Since heat escape is \(S_H=O(t)\),

\[
\boxed{\gamma_H=1}.
\]

The exponent is unchanged.

### 8.2 Unitary-like channel

The first-order transition amplitude contains

\[
\int_0^t A(s)ds
=
[A+\varepsilon J]t
-\frac{\varepsilon J}{2}t^2
+O(t^3).
\]

After re-centering by the static perturbed leading profile,
the normalized destination shape acquires an \(O(\varepsilon t)\) transverse
term.

But

\[
S_U=O(t^2).
\]

Therefore generically

\[
\boxed{\gamma_U^{\rm relaxing\ hidden}=\frac12}.
\]

This is crucial: a relaxing hidden state can make a unitary-like record carry
the same projective exponent that a static wave-like record had previously.

So \(\gamma=1/2\) by itself is not a stable channel label in the presence of
unspecified dynamic hidden coupling.

### 8.3 Wave-like channel

For

\[
q''(t)+A(t)q(t)=0,\qquad q'(0)=0,
\]

the off-diagonal amplitude expands as

\[
q_j(t)
=
-\frac12(A+\varepsilon J)_{ji}t^2
+\frac{\varepsilon J_{ji}}{6}t^3
+O(t^4).
\]

After re-centering, the profile therefore changes as \(O(\varepsilon t)\).

Since

\[
S_W=O(t^4),
\]

generically

\[
\boxed{\gamma_W^{\rm relaxing\ hidden}=\frac14}.
\]

Thus dynamic hidden memory lowers the projective exponent by exposing a
lower-order transverse shape correction.

---

## 9. Numerical replay on the strict operator

For one nontrivial hidden direction

\[
h=(0.7,-0.4,0.5,0.2)
\]

and test coupling \(\varepsilon=0.05\), direct ODE integration gives local
log-log slopes approaching:

- heat: \(\gamma\approx1.00002\),
- unitary: \(\gamma\approx0.504\),
- wave: \(\gamma\approx0.265\),

with the latter two moving toward the analytic asymptotes

\[
1/2,\qquad1/4
\]

as the window is reduced.

These numbers are numerical checks of the analytic short-time orders, not
physical predictions.

---

## 10. The important physics-side conclusion

The current evidence now separates three statements:

### A. Relational time without scale

The pair \((S,r)\) can encode clock-free ordering information.

### B. Absolute time remains absent

Global rate rescaling remains a gauge; no SI clock emerges.

### C. Dynamic hidden memory can change the relational exponent

Therefore the projective exponent is not a property of the static operator
alone.

It belongs to a typed structure such as

\[
(\text{generator},\text{preparation},\text{hidden coupling},
 \text{observation law}).
\]

This is a stronger and more precise version of the earlier statement that
"time may emerge from transformations of relations": the *form* of temporal
change can be relational, but only after the memory/coupling law is specified.

---

## 11. What the candidate bridge does and does not accomplish

It does provide:

- a D12-equivariant hidden-to-generator map;
- conservation of the zero mode;
- symmetric Laplacian structure;
- all-four-mode injectivity;
- a concrete prediction for how relaxing hidden memory changes clock-free
  exponents.

It does not provide:

- a derivation of the coupling constant;
- a proof that endpoint fractional modulation is selected by FIN;
- a physical clock;
- a Born rule;
- a unique microscopic interpretation;
- QW-2191 closure;
- role-bearing \(L_{\rm total}\);
- SM/GR or ToE closure.

---

## 12. Next decisive atom

### DYNAMIC-HIDDEN-BRIDGE-SOURCE-35

Do not optimize the candidate bridge further.

Instead ask whether the endpoint fractional rule can be derived from an
already admitted FIN object.

Two concrete routes should be tested:

1. **state-to-multiplication route**
   \[
   u\mapsto D_u,\qquad
   (A,D_u)\mapsto J(u),
   \]
   checking whether conservation, symmetry, and a variational/Dirichlet
   principle force the Laplacianized anticommutator above;

2. **tree-response route**
   derive the first boundary-operator variation caused by hidden storage or
   conductance modulation and compare it with \(J(u)\).

Acceptance:
- a theorem deriving the bridge up to one scalar from an admitted principle;
  or
- a no-go showing at least two inequivalent admissible bridges survive the
  same principles.

That result is more important than further fitting of the coupling strength.
