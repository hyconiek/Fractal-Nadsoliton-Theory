# N474 Current First Strict Vacuum EoM Eliminates Yukawa From the Canonical Local-Diagonal Residual (Conditional) Theorem

Status: `N474_DISCHARGED_CURRENT_FIRST_STRICT_VACUUM_EOM_ELIMINATES_YUKAWA_FROM_CANONICAL_LOCAL_DIAGONAL_RESIDUAL_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

`T166` (and more broadly the strict `QW-2191` frontier on `pair1`) is blocked because the **canonical** FIN local diagonal
residual profile

$$
D_{\mathrm{local\_residual}}=\mathrm{diag}(d_0,\ldots,d_{11})
$$

is currently exported only as a coefficient class (`R15/R18`), with `F2(d)` underdetermined (`N472/P431`).

This theorem executes a narrow honest *reduction* step:

```text
under the exported canonical vacuum EoM (QW-2165), the Yukawa contribution cancels out of the diagonal Hessian entry
used to define d_k (QW-2166/R15), provided vpsi_k ≠ 0.
```

So the canonical diagonal residual entries can be rewritten in a **Yukawa-free** form once vacuum stationarity is imposed.

This is structural only: it does **not** decide the canonical diagonal profile, and it does **not** discharge `T166`
or `QW-2191`.

## Strict-admissible evidence reused

1. `QW-2165`
   - exhaustive canonical EoM for `12xPsi + Phi` (exports each `eom_psi_k`).
2. `QW-2166`
   - exhaustive canonical Hessian / linearized EoM (exports the diagonal coefficient
     `3*g4_psik*vpsik**2 + 5*g6_psik*vpsik**4 + 2*gYk*vphi**2 + m2_psik` in the `eta_k` row stencil).
3. `R15`
   - defines the canonical diagonal decomposition
     $D_{\mathrm{canonical}} = m_0^2 I + D_{\mathrm{local\_residual}}$
     with $d_k = (D_{\mathrm{canonical}})_{kk} - m_0^2$.

## Setup (canonical constant vacuum specialization)

For each site index $k\in\{0,\ldots,11\}$, `QW-2165` exports a canonical EoM of the form:

$$
0 = \sum_{j\neq k}\frac{K_{k,j}+K_{j,k}}{2}\,\psi_j(x)
  + g4_{\psi k}\,\psi_k(x)^3
  + g6_{\psi k}\,\psi_k(x)^5
  + 2\,gY_k\,\phi(x)^2\,\psi_k(x)
  + m2_{\psi k}\,\psi_k(x)
  + \frac{d^2}{dx^2}\psi_k(x).
$$

Consider a *constant vacuum* specialization:

$$
\psi_k(x)\equiv v_{\psi k},\qquad \phi(x)\equiv v_\phi,
$$

so $\frac{d^2}{dx^2}\psi_k(x)\equiv 0$ and the stationarity condition becomes:

$$
0 = \sum_{j\neq k}\frac{K_{k,j}+K_{j,k}}{2}\,v_{\psi j}
  + g4_{\psi k}\,v_{\psi k}^3
  + g6_{\psi k}\,v_{\psi k}^5
  + 2\,gY_k\,v_\phi^2\,v_{\psi k}
  + m2_{\psi k}\,v_{\psi k}.
$$

Also, `QW-2166` exports the diagonal Hessian (linearized) coefficient on the `eta_k` row:

$$
(D_{\mathrm{canonical}})_{kk}
=
3\,g4_{\psi k}\,v_{\psi k}^2
+
5\,g6_{\psi k}\,v_{\psi k}^4
+
2\,gY_k\,v_\phi^2
+
m2_{\psi k}.
$$

By `R15`, the canonical local diagonal residual entry is:

$$
d_k := (D_{\mathrm{local\_residual}})_{kk} = (D_{\mathrm{canonical}})_{kk} - m_0^2.
$$

## Theorem (Yukawa elimination from the diagonal entry)

Fix $k\in\{0,\ldots,11\}$ and assume $v_{\psi k}\neq 0$.

### Claim 1. The vacuum EoM solves for the combination $(m2_{\psi k} + 2 gY_k v_\phi^2)$.

Divide the constant-vacuum stationarity equation by $v_{\psi k}$:

$$
0
=
\sum_{j\neq k}\frac{K_{k,j}+K_{j,k}}{2}\,\frac{v_{\psi j}}{v_{\psi k}}
 g4_{\psi k}\,v_{\psi k}^2
 g6_{\psi k}\,v_{\psi k}^4
 (m2_{\psi k}+2\,gY_k\,v_\phi^2).
$$

So:

$$
m2_{\psi k}+2\,gY_k\,v_\phi^2
=
-\sum_{j\neq k}\frac{K_{k,j}+K_{j,k}}{2}\,\frac{v_{\psi j}}{v_{\psi k}}
 -g4_{\psi k}\,v_{\psi k}^2
 -g6_{\psi k}\,v_{\psi k}^4.
$$

### Claim 2. Substituting into the diagonal Hessian entry cancels Yukawa.

Start from the exported diagonal entry:

$$
(D_{\mathrm{canonical}})_{kk}
=
3\,g4_{\psi k}\,v_{\psi k}^2
+
5\,g6_{\psi k}\,v_{\psi k}^4
+
(m2_{\psi k}+2\,gY_k\,v_\phi^2).
$$

Substitute Claim 1:

$$
(D_{\mathrm{canonical}})_{kk}
=
3\,g4_{\psi k}\,v_{\psi k}^2
5\,g6_{\psi k}\,v_{\psi k}^4
\;-\sum_{j\neq k}\frac{K_{k,j}+K_{j,k}}{2}\,\frac{v_{\psi j}}{v_{\psi k}}
\;-\;g4_{\psi k}\,v_{\psi k}^2
\;-\;g6_{\psi k}\,v_{\psi k}^4.
$$

Therefore:

$$
(D_{\mathrm{canonical}})_{kk}
=
-\sum_{j\neq k}\frac{K_{k,j}+K_{j,k}}{2}\,\frac{v_{\psi j}}{v_{\psi k}}
\;+\;2\,g4_{\psi k}\,v_{\psi k}^2
\;+\;4\,g6_{\psi k}\,v_{\psi k}^4.
$$

In particular, the Yukawa parameters $(gY_k, v_\phi)$ no longer appear.

### Corollary (Yukawa-free residual diagonal entry)

By `R15`:

$$
d_k
=
(D_{\mathrm{canonical}})_{kk}-m_0^2
=
-\sum_{j\neq k}\frac{K_{k,j}+K_{j,k}}{2}\,\frac{v_{\psi j}}{v_{\psi k}}
\;+\;2\,g4_{\psi k}\,v_{\psi k}^2
\;+\;4\,g6_{\psi k}\,v_{\psi k}^4
\;-\;m_0^2.
$$

So, *under constant-vacuum stationarity and $v_{\psi k}\neq 0$*, the canonical local diagonal residual is Yukawa-free.

## What N474 does not prove

`N474` does not prove:

1. that the canonical vacuum satisfies $v_{\psi k}\neq 0$ for all $k$,
2. any strict-derived vacuum vector $(v_{\psi 0},\ldots,v_{\psi 11})$,
3. any strict-derived numeric instantiation of $K_{k,j}$ or the self-coupling families $g4_{\psi k}, g6_{\psi k}$,
4. any decision of the diagonal mode‑2 defect $F_2(d)$ (`T166` remains open),
5. any strict-core `pair1` `O(2)` cut, strict theta export, strict-core selector closure, or `QW-2191` discharge,
6. ToE closure.

## Consequence (most honest use)

This theorem is a reduction tool: any future strict-derived vacuum/EoM solving (or a value-instantiation object)
deciding the canonical diagonal residual profile can be checked against a Yukawa-free diagonal formula, reducing the
unknown dependency set in the `T166` diagonal-defect decision target.

