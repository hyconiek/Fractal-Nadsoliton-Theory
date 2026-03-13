# N473 Current First Strict Canonical Local-Diagonal Mode‑2 Defect Recovered From the `R18` `pair1` Entry System Theorem

Status: `N473_DISCHARGED_CURRENT_FIRST_STRICT_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_RECOVERED_FROM_R18_PAIR1_ENTRY_SYSTEM_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

Two strict lanes currently coexist around the canonical FIN local diagonal residual sector:

1. the **host-matching lane** exports the exact declared `pair1` residual pullback entry system in terms of the six
   opposite-pair sums (`R18`), and
2. the **diagonal accelerator lane** exports the strict `pair1` `O(2)`‑cut criterion via the mode‑2 defect
   `F2(d)` (`N466`) and the explicit six‑class reduction for the canonical sector (`N467/P426`).

`N473` provides the missing strict glue:

```text
the explicit mode‑2 defect (Re/Im) is exactly recoverable from the R18 pair1-entry system:
  Re(F2) = 6*(c1c1 - s1s1),
  Im(F2) = 12*(c1s1),
and substituting the R18 entry formulas reproduces the N467/P426 six-class reduction.
```

This is a pure structural identity. It does **not** decide the canonical diagonal profile (see `T166` and `N472/P431`).

## Strict-admissible evidence reused

1. `N466`
   - diagonal `pair1` `O(2)`‑cut criterion and closed-form relation between the anisotropy signature and `F2(d)`.
2. `R18`
   - exact coefficient-class reduction of the declared `pair1` block of the canonical local diagonal residual sector,
     expressed via six opposite-pair sums `Sigma_psi0_psi6, ..., Sigma_psi5_psi11`.
3. `N467/P426`
   - strict six-class reduction of `F2(d)` for the canonical FIN local diagonal residual sector (symbolic on values).

## Setup

Let `n=12`. Let the canonical FIN local diagonal residual sector be:

$$
D_{\mathrm{local\_residual}}=\mathrm{diag}(d_0,\ldots,d_{11}).
$$

Define the mode‑2 defect (as in `N466`):

$$
F_2(d):=\sum_{k=0}^{11} d_k\,e^{i\frac{4\pi k}{12}}\in\mathbb{C}.
$$

Let the `pair1` restriction matrix in basis `(c1,s1)` be:

$$
M_1(D_{\mathrm{local\_residual}})=
\begin{pmatrix}
a_1 & b_1\\
b_1 & d_1
\end{pmatrix},
\qquad
a_1:=\langle c_1, D c_1\rangle,\;
b_1:=\langle c_1, D s_1\rangle,\;
d_1:=\langle s_1, D s_1\rangle.
$$

`R18` exports (by transport reduction) the same three entries as explicit linear combinations of the six opposite-pair
sums:

$$
\Sigma_{0}:=\Sigma_{\psi0,\psi6},\;
\Sigma_{1}:=\Sigma_{\psi1,\psi7},\;
\Sigma_{2}:=\Sigma_{\psi2,\psi8},\;
\Sigma_{3}:=\Sigma_{\psi3,\psi9},\;
\Sigma_{4}:=\Sigma_{\psi4,\psi10},\;
\Sigma_{5}:=\Sigma_{\psi5,\psi11},
$$

where each $\Sigma_k$ is the residual opposite‑pair sum $d_k+d_{k+6}$ (expanded in `R18` as a symbolic coefficient-class
expression on the diagonal coefficients).

Concretely, `R18` gives:

$$
a_1
=
\frac{1}{6}\Sigma_0 + \frac{1}{8}\Sigma_1 + \frac{1}{24}\Sigma_2 + \frac{1}{24}\Sigma_4 + \frac{1}{8}\Sigma_5,
$$

$$
b_1
=
\frac{\sqrt{3}}{24}(\Sigma_1+\Sigma_2-\Sigma_4-\Sigma_5),
$$

$$
d_1
=
\frac{1}{24}\Sigma_1 + \frac{1}{8}\Sigma_2 + \frac{1}{6}\Sigma_3 + \frac{1}{8}\Sigma_4 + \frac{1}{24}\Sigma_5.
$$

## Theorem-level claims

### Claim 1. Recovering `Re(F2)` and `Im(F2)` from the `pair1` block entries

By `N466` (here `n=12`):

$$
a_1-d_1=\frac{1}{6}\mathrm{Re}\,F_2(d),\qquad b_1=\frac{1}{12}\mathrm{Im}\,F_2(d).
$$

Therefore:

$$
\mathrm{Re}\,F_2(d)=6(a_1-d_1),\qquad \mathrm{Im}\,F_2(d)=12\,b_1.
$$

### Claim 2. Substituting the `R18` entry system reproduces the `N467/P426` six-class reduction

Using the `R18` formulas for $(a_1,b_1,d_1)$:

$$
6(a_1-d_1)
=
\Sigma_0 + \frac{1}{2}\Sigma_1 - \frac{1}{2}\Sigma_2 - \Sigma_3 - \frac{1}{2}\Sigma_4 + \frac{1}{2}\Sigma_5,
$$

$$
12\,b_1
=
\frac{\sqrt{3}}{2}(\Sigma_1+\Sigma_2-\Sigma_4-\Sigma_5).
$$

So:

$$
F_2(d)
=
\left(\Sigma_0 + \frac{1}{2}\Sigma_1 - \frac{1}{2}\Sigma_2 - \Sigma_3 - \frac{1}{2}\Sigma_4 + \frac{1}{2}\Sigma_5\right)
+ i\,\frac{\sqrt{3}}{2}(\Sigma_1+\Sigma_2-\Sigma_4-\Sigma_5),
$$

which is exactly the six-class reduction exported as the canonical diagonal mode‑2 defect in `N467/P426`.

So the diagonal accelerator lane and the `R18` host-matching entry system are strictly compatible: they are two views
of the same underlying `pair1` restriction data of the canonical diagonal residual sector.

### Corollary 1. If the `pair1` residual block vanishes (host-matching target), then `F2(d)=0`

If the `R18` host-matching target conditions hold:

$$
a_1=0,\qquad b_1=0,\qquad d_1=0,
$$

then by Claim 1:

$$
\mathrm{Re}\,F_2(d)=6(a_1-d_1)=0,\qquad \mathrm{Im}\,F_2(d)=12\,b_1=0,
$$

so $F_2(d)=0$, hence by `N466` the diagonal/local sector cannot cut `O(2)` on `pair1`.

Equivalently:

```text
any strict diagonal pair1 O(2)-cut witness from the canonical FIN D_local_residual implies
that at least one of the R18 pair1 entry equations fails (the residual block does not vanish).
```

This is a conditional structural statement only. It does not claim any of the `R18` zero equations hold today.

## What N473 does not prove

`N473` does not prove:

1. any strict-derived values for the canonical diagonal residual profile (see `T166`, `N472/P431`),
2. that the canonical FIN local diagonal residual sector cuts `O(2)` on `pair1` (`F2(d)≠0`),
3. that the canonical FIN local diagonal residual sector fails to cut `O(2)` on `pair1` (`F2(d)=0`),
4. strict-core theta export, strict-core selector closure, or global `QW-2191` discharge,
5. ToE closure.

