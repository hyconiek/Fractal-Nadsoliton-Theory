# N484 Current First Strict Diagonal Sector Pair‑m O(2)‑Cut Criterion + Canonical Eigenbasis Reconstruction From `F_{2m}` (n=12) Theorem

Status: `N484_DISCHARGED_CURRENT_FIRST_STRICT_DIAGONAL_SECTOR_PAIR_M_O2_CUT_CRITERION_AND_EIGENBASIS_FROM_F2M_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Generalize the already exported `pair1` diagonal‑sector `O(2)` cut facts (`N466/N468`) to an arbitrary Fourier‑degenerate
pair `pair_m` on the strict `n=12` ring scaffold (`QW-2190/QW-2191`):

```text
for a diagonal profile d_k (k=0..11),
the diagonal sector cuts O(2) on pair_m  ⇔  the Fourier defect F_{2m}(d) is nonzero,
and the canonical diagonalization angle/eigenvalues are explicit functions of F_{2m}(d).
```

This theorem is purely about the **diagonal** operator (a typed diagonal/local sector).
It does not claim global `QW-2191` discharge by itself.

## Setup (n=12 Fourier pairs + diagonal operator)

Let `n=12` and let `d = (d_0,...,d_{11}) ∈ R^{12}`.
Define the diagonal operator on `R^{12}`:

```text
D := diag(d_0,...,d_{11}).
```

For `m ∈ {1,2,3,4,5}`, define the normalized real Fourier pair:

```text
c_m(k) := sqrt(2/n) * cos(2π m k / n),
s_m(k) := sqrt(2/n) * sin(2π m k / n),
```

and the `pair_m` plane:

```text
V_m := span{c_m, s_m} ⊂ R^n.
```

Define the complex Fourier defect at mode `2m`:

$$
F_{2m}(d)
:=
\sum_{k=0}^{n-1} d_k\,e^{i\,2\pi(2m)k/n}
\in \mathbb{C}.
$$

Write `F_{2m}(d) = Re(F_{2m}) + i Im(F_{2m})`, and define the diagonal mean:

$$
\mu(d) := \frac{1}{n}\sum_{k=0}^{n-1} d_k.
$$

## Theorem (restriction to `pair_m`)

### Claim 1. The restriction `D|_{V_m}` has an explicit `2×2` matrix in the basis `(c_m,s_m)`.

In the ordered basis `(c_m, s_m)`, the matrix of `D|_{V_m}` is:

$$
\begin{pmatrix}
\mu(d) + \frac{1}{n}\,Re(F_{2m}(d)) & \frac{1}{n}\,Im(F_{2m}(d)) \\
\frac{1}{n}\,Im(F_{2m}(d)) & \mu(d) - \frac{1}{n}\,Re(F_{2m}(d))
\end{pmatrix}.
$$

*Proof sketch.* Use:

$$
\cos^2\theta = \frac12(1+\cos 2\theta),\quad
\sin^2\theta = \frac12(1-\cos 2\theta),\quad
\cos\theta\sin\theta = \frac12\sin 2\theta,
$$

with $\theta = 2\pi m k/n$, and the normalization $c_m(k)^2+s_m(k)^2 = 2/n$. ∎

### Claim 2. `D` cuts `O(2)` on `pair_m` iff `|F_{2m}(d)| ≠ 0`.

If `F_{2m}(d) = 0`, the above matrix reduces to `μ(d) I_2`, hence the diagonal sector is isotropic on `V_m`
and cannot cut the internal `O(2)` family on `pair_m`.

If `|F_{2m}(d)| ≠ 0`, the eigenvalues are distinct, so the eigenbasis is unique up to the residual `Z2` sign
convention in each eigenvector; therefore the diagonal sector cuts `O(2)` down to `Z2` on `pair_m`.
∎

### Claim 3. Canonical eigenvalues and diagonalization angle are explicit.

If `|F_{2m}(d)| ≠ 0`, then:

1. eigenvalues:

$$
\lambda_\pm
=
\mu(d) \pm \frac{1}{n}|F_{2m}(d)|,
$$

2. diagonalization angle (one canonical representative):

$$
\theta_*(m)
=
\frac12\,\operatorname{atan2}\!\left(Im(F_{2m}(d)),\,Re(F_{2m}(d))\right),
$$

so that the eigenvectors are obtained from `(c_m,s_m)` by a rotation by `θ_*(m)`.
∎

## What N484 does not claim

`N484` does not claim:

1. that any particular diagonal/local profile `d_k` exists in strict core (that is a separate provider task),
2. global `QW-2191` discharge,
3. strict-core theta export (`T159`) or sigma-int slot elimination (`T160/T161/T162`),
4. ToE closure.

