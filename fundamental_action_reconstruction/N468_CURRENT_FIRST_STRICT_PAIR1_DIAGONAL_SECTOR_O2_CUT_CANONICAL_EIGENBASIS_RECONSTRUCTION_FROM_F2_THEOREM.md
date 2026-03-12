# N468 Current First Strict `pair1` Diagonal-Sector `O(2)`-Cut Canonical Eigenbasis Reconstruction From `F2` Theorem

Status: `N468_DISCHARGED_CURRENT_FIRST_STRICT_PAIR1_DIAGONAL_SECTOR_O2_CUT_CANONICAL_EIGENBASIS_RECONSTRUCTION_FROM_F2_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After the strict diagonal criterion (`N466`) and the `n=12` reduction (`N467`), the remaining strictly relevant
connection to the `QW-2191` uniqueness obstruction is:

```text
if a strict diagonal/local sector has F2(d) ≠ 0,
it must not only “break O(2)”, but also supply a canonical axis/eigenbasis on pair1.
```

This theorem exports that **explicit canonicalization formula**, without claiming that any currently exported diagonal
profile actually has `F2(d) ≠ 0`.

## Strict-admissible evidence reused

1. `QW-2190`
   - the 12-octave ring and the real Fourier `pair1 = span{c1,s1}` scaffold.
2. `N466`
   - diagonal `pair1` `O(2)`-cut criterion and the closed-form relation between the `pair1` anisotropy signature and
     the mode-2 Fourier coefficient `F2(d)`.
3. `N467`
   - for `n=12`, the canonical diagonal mode-2 defect `F2(d)` reduces to the six opposite-pair sums.

## Setup

Let `n=12` (the strict carrier in `QW-2190`). Let:

$$
D=\mathrm{diag}(d_0,\ldots,d_{n-1})
$$

be a diagonal operator in the site basis.

Let `pair1 = span{c1,s1}` be the exported real Fourier `m=1` mode subspace, with the standard orthonormal scaffold:

$$
c_{1,k}:=\sqrt{\frac{2}{n}}\cos\!\left(\frac{2\pi k}{n}\right),\qquad
s_{1,k}:=\sqrt{\frac{2}{n}}\sin\!\left(\frac{2\pi k}{n}\right),
\qquad k=0,\ldots,n-1.
$$

Define the complex mode-2 Fourier coefficient of the diagonal profile:

$$
F_2(d):=\sum_{k=0}^{n-1} d_k\,e^{i\frac{4\pi k}{n}}\in\mathbb{C}.
$$

Write its real/imaginary parts as:

$$
F_2(d)=\mathrm{Re}\,F_2(d)+i\,\mathrm{Im}\,F_2(d).
$$

## Theorem-level claims

### Claim 1. The `pair1` restriction diagonalizes canonically when `F2(d)≠0`.

Let the `pair1` block in basis `(c1,s1)` be:

$$
\left[D\right]_{\{c_1,s_1\}}
=
\begin{pmatrix}
a_1 & b_1\\
b_1 & d_1
\end{pmatrix}.
$$

By `N466`:

$$
a_1-d_1=\frac{2}{n}\,\mathrm{Re}\,F_2(d),\qquad
b_1=\frac{1}{n}\,\mathrm{Im}\,F_2(d),
$$

and (since $c_1,s_1$ are orthonormal) the trace is:

$$
a_1+d_1
=
c_1^\top D c_1 + s_1^\top D s_1
=
\sum_{k=0}^{n-1} d_k\,(c_{1,k}^2+s_{1,k}^2)
=
\frac{2}{n}\sum_{k=0}^{n-1} d_k.
$$

Therefore the eigenvalues of the `pair1` block are:

$$
\lambda_\pm
=
\frac{a_1+d_1}{2}\pm\sqrt{\left(\frac{a_1-d_1}{2}\right)^2+b_1^2}
=
\frac{1}{n}\sum_{k=0}^{n-1} d_k \pm \frac{1}{n}\left|F_2(d)\right|.
$$

In particular, if $F_2(d)\neq 0$, then $\lambda_+\neq\lambda_-$ and the `pair1` eigenbasis is unique up to the
residual sign flips (and the trivial swap induced by ordering the two eigenvalues).

### Claim 2. The canonical `pair1` axis angle is half the argument of `F2(d)`.

Assume $F_2(d)\neq 0$ and define:

$$
\theta_* := \frac12\,\operatorname{atan2}\!\left(\mathrm{Im}\,F_2(d),\ \mathrm{Re}\,F_2(d)\right).
$$

Equivalently, $\theta_*$ is any angle satisfying:

$$
\tan(2\theta_*)=\frac{\mathrm{Im}\,F_2(d)}{\mathrm{Re}\,F_2(d)}.
$$

Define the rotated basis:

$$
\begin{pmatrix}c_1' & s_1'\end{pmatrix}

=
\begin{pmatrix}c_1 & s_1\end{pmatrix}
\begin{pmatrix}
\cos\theta_* & -\sin\theta_*\\
\sin\theta_* & \cos\theta_*
\end{pmatrix},
$$

i.e.

$$
c_1'=\cos\theta_*\,c_1+\sin\theta_*\,s_1,\qquad
s_1'=-\sin\theta_*\,c_1+\cos\theta_*\,s_1.
$$

Then the `pair1` block diagonalizes in this basis:

$$
\left[D\right]_{\{c_1',s_1'\}}
=
\begin{pmatrix}
\lambda_+ & 0\\
0 & \lambda_-
\end{pmatrix}.
$$

So, given a strict diagonal/local profile with $F_2(d)\neq 0$, the canonical `pair1` axis (and hence the `O(2)` cut)
is recovered explicitly by the eigenbasis rotation angle $\theta_*=\frac12\arg(F_2(d))$.

### Claim 3. This supplies a canonical `O(2)` cut on `pair1` (conditional on `F2(d)≠0`).

If $F_2(d)\neq 0$, then:

1. `D|pair1` is not scalar and therefore breaks `O(2)` on `pair1` (`N466`),
2. the induced eigenbasis $(c_1',s_1')$ is canonical up to:
   - the unavoidable sign flips $c_1'\mapsto -c_1'$, $s_1'\mapsto -s_1'$,
   - and the swap $(c_1',s_1')\leftrightarrow(s_1',c_1')$ (fixed canonically by ordering eigenvalues).

Therefore a strict-derived diagonal/local sector with $F_2(d)\neq 0$ would constitute a concrete strict mechanism to
cut the `QW-2191` `O(2)` family on `pair1` (in the declared scope), without any additional continuous selector slot.

## What N468 does not prove

`N468` does not prove:

1. that the canonical FIN local diagonal residual sector has $F_2(d)\neq 0$,
2. any strict-derived coefficient instantiation of $d_k$,
3. a strict-core theta export, strict-core selector closure, or global `QW-2191` discharge,
4. ToE closure.

It only gives the explicit reconstruction formula that becomes applicable **if** a strict diagonal/local mode‑2 defect
is ever derived.

