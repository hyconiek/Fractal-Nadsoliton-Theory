# N471 Current First Strict `Aut(Z_12)`‑Invariant Diagonal Profile Mode‑2 Defect Orbit‑Reduction ⇒ `pair1` `O(2)`‑Cut Criterion Theorem

Status: `N471_DISCHARGED_CURRENT_FIRST_STRICT_AUT_Z12_INVARIANT_DIAGONAL_PROFILE_MODE2_DEFECT_ORBIT_REDUCTION_PAIR1_O2_CUT_CRITERION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

`N466` reduces any strict diagonal/local `pair1` accelerator attempt to one explicit condition:

```text
diagonal/local sector cuts O(2) on pair1  <=>  F2(d) ≠ 0.
```

Separately, the repo exports the typed `Phase_12_v1/Aut(Z_12_v1)` quotient-orbit discipline (`N455`), which is a
slot‑hygiene move: it allows stating Aut‑invariant data on a canonical 6‑point quotient carrier without smuggling a
hidden generator/orientation choice.

This theorem connects the two lanes in a strictly checkable way:

```text
If a diagonal profile d is Aut(Z_12)-invariant (i.e. constant on the 6 quotient orbits from N455),
then F2(d) is forced to be real and reduces to one explicit orbit‑linear combination.
Therefore the diagonal pair1 O(2)-cut criterion becomes a single explicit scalar check on 6 orbit values.
```

This is a pure structural statement. It does **not** claim that the canonical FIN local diagonal residual sector is
Aut‑invariant, and it does not provide strict-derived values for any orbit class.

## Strict-admissible evidence reused

1. `N466`
   - diagonal `pair1` `O(2)`‑cut criterion via $F_2(d)$.
2. `N455`
   - explicit 6‑orbit decomposition of `Phase_12_v1` under `Aut(Z_12_v1)`.
3. `QW-2190`
   - strict `n=12` ring carrier.

## Setup

Let `n=12`. Let $D=\mathrm{diag}(d_0,\ldots,d_{11})$ be a diagonal operator in the site basis.

Define the mode‑2 Fourier coefficient (as in `N466`):

$$
F_2(d):=\sum_{k=0}^{n-1} d_k\,e^{i\frac{4\pi k}{n}}\in\mathbb{C}.
$$

We say the diagonal profile is **Aut-invariant** iff it is constant on the 6 orbits exported by `N455`, i.e.:

$$
d_{u k} = d_k \quad\text{for all }u\in (\mathbb{Z}/12\mathbb{Z})^{\times}=\{1,5,7,11\}\text{ and all }k\in\{0,\ldots,11\},
$$

equivalently (re-writing the `N455` orbits on exponents):

$$
\{0\},\ \{6\},\ \{1,5,7,11\},\ \{2,10\},\ \{3,9\},\ \{4,8\}.
$$

Define the corresponding orbit values:

$$
\begin{aligned}
&t_0:=d_0,\qquad t_6:=d_6,\qquad t_1:=d_1=d_5=d_7=d_{11},\\
&t_2:=d_2=d_{10},\qquad t_3:=d_3=d_9,\qquad t_4:=d_4=d_8.
\end{aligned}
$$

## Theorem-level claims

### Claim 1. For an Aut-invariant diagonal profile on `n=12`, $F_2(d)$ is real and equals:

$$
F_2(d)=t_0+t_6+2t_1-t_2-2t_3-t_4\in\mathbb{R}.
$$

**Proof.** Write:

$$
F_2(d)=\sum_{k=0}^{11} d_k e^{i\frac{4\pi k}{12}}=\sum_{k=0}^{11} d_k e^{i\frac{\pi k}{3}}.
$$

Group by the 6 `Aut(Z_12)` orbits and substitute the orbit constants.
The exponential sums over each orbit are:

$$
\sum_{k\in\{0\}} e^{i\frac{\pi k}{3}} = 1,\qquad
\sum_{k\in\{6\}} e^{i\frac{\pi k}{3}} = e^{i2\pi}=1,
$$

$$
\sum_{k\in\{1,5,7,11\}} e^{i\frac{\pi k}{3}}
=
\left(e^{i\frac{\pi}{3}}+e^{i\frac{7\pi}{3}}\right)+\left(e^{i\frac{5\pi}{3}}+e^{i\frac{11\pi}{3}}\right)
=2\left(e^{i\frac{\pi}{3}}+e^{i\frac{5\pi}{3}}\right)=2,
$$

$$
\sum_{k\in\{2,10\}} e^{i\frac{\pi k}{3}}
=e^{i\frac{2\pi}{3}}+e^{i\frac{10\pi}{3}}
=e^{i\frac{2\pi}{3}}+e^{i\frac{4\pi}{3}}=-1,
$$

$$
\sum_{k\in\{3,9\}} e^{i\frac{\pi k}{3}}
=e^{i\pi}+e^{i3\pi}=-2,
$$

$$
\sum_{k\in\{4,8\}} e^{i\frac{\pi k}{3}}
=e^{i\frac{4\pi}{3}}+e^{i\frac{8\pi}{3}}
=e^{i\frac{4\pi}{3}}+e^{i\frac{2\pi}{3}}=-1.
$$

Therefore:

$$
F_2(d)=t_0\cdot 1+t_6\cdot 1+t_1\cdot 2+t_2\cdot (-1)+t_3\cdot (-2)+t_4\cdot (-1)
=t_0+t_6+2t_1-t_2-2t_3-t_4,
$$

which is real. ∎

### Claim 2. Orbit-reduced strict diagonal `pair1` `O(2)`‑cut criterion

By `N466`:

$$
D\ \text{cuts `O(2)` on `pair1`}\quad\Longleftrightarrow\quad F_2(d)\neq 0.
$$

Under Aut-invariance, Claim 1 gives the explicit scalar test:

$$
D\ \text{cuts `O(2)` on `pair1`}
\quad\Longleftrightarrow\quad
t_0+t_6+2t_1-t_2-2t_3-t_4\neq 0.
$$

So, for quotient-safe Aut-invariant diagonal profiles, the strict diagonal accelerator question reduces to one explicit
check on the **six** orbit values.

## Corollary (parity-only profiles cannot cut `O(2)` on `pair1`)

If the diagonal profile depends **only** on the parity class $k\bmod 2$ (equivalently, only on
`chi_parity` from `N454`), then:

$$
t_0=t_6=t_2=t_4\ (=t_{\mathrm{even}}),\qquad t_1=t_3\ (=t_{\mathrm{odd}}),
$$

and therefore:

$$
F_2(d)=t_{\mathrm{even}}+t_{\mathrm{even}}+2t_{\mathrm{odd}}-t_{\mathrm{even}}-2t_{\mathrm{odd}}-t_{\mathrm{even}}=0.
$$

So any parity-only diagonal/local dependence is **too weak** to cut `O(2)` on `pair1` by the strict diagonal criterion
`N466` (consistent with the special sign-mask closure `N469`).

## What N471 does not prove

`N471` does not prove:

1. that the canonical FIN local diagonal residual sector is Aut-invariant,
2. any strict-derived orbit values $t_0,t_6,t_1,t_2,t_3,t_4$ for the canonical diagonal residual sector,
3. that the canonical diagonal residual sector cuts `O(2)` on `pair1`,
4. any strict-core theta export, strict-core selector closure, or global `QW-2191` discharge,
5. ToE closure.
