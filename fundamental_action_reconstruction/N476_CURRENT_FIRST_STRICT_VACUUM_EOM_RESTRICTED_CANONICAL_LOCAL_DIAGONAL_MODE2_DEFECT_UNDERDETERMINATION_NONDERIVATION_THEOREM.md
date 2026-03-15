# N476 Current First Strict Vacuum‑EoM‑Restricted Canonical Local‑Diagonal Mode‑2 Defect Underdetermination (Nonderivation) Theorem

Status: `N476_DISCHARGED_CURRENT_FIRST_STRICT_VACUUM_EOM_RESTRICTED_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_UNDERDETERMINATION_NONDERIVATION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

`T166` requires a strict-derived decision of the diagonal/local mode‑2 defect:

$$
F_2(d)=\sum_{k=0}^{11} d_k\,e^{i\pi k/3}
$$

for the **canonical** FIN local diagonal residual profile $d_k=(D_{\mathrm{local\_residual}})_{kk}$.

`N472/P431` prove that $F_2(d)$ is underdetermined at the exported coefficient‑class level (`R15`).

After `N474/N475`, a common temptation remains:

```text
maybe imposing constant-vacuum stationarity (canonical EoM) plus vpsi_k≠0 already decides F2(d)
```

`N476` closes that temptation honestly:

```text
even after restricting to constant-vacuum stationarity (so N474/N475 apply),
the canonical diagonal/local mode‑2 defect remains underdetermined on the current export class.
```

So `T166` is still not discharged without exporting additional strict-derived constraints/value objects beyond
stationarity.

## Strict-admissible evidence reused

1. `QW-2165`
   - canonical EoM (for constant-vacuum stationarity premise).
2. `R15`
   - canonical diagonal residual sector definition $d_k$.
3. `N474`
   - under stationarity and $v_{\psi k}\neq 0$, Yukawa cancels from $d_k$.
4. `N475`
   - induced Yukawa-free rewrite at the opposite‑pair sum level.
5. `N467`
   - reduction: the `pair1` diagonal/local `O(2)` cut depends only on the six opposite‑pair sums.

## Setup

Assume constant-vacuum stationarity in the sense of `N474`:

$$
\psi_k(x)\equiv v_{\psi k},\qquad \phi(x)\equiv v_\phi,
$$

and assume $v_{\psi k}\neq 0$ for the indices under discussion so the Yukawa-free rewrite of `N474/N475` applies.

## Theorem (stationarity does not decide $F_2(d)$ on current exports)

Even under the stationarity/nonzero premises above, the diagonal residual profile still depends on free local
coefficient families (e.g. $(g4_{\psi k})$, $(g6_{\psi k})$) and on vacuum data $(v_{\psi k})$ for which the repo does
not export strict-derived constraints sufficient to decide the six opposite‑pair sums
`Sigma_psi0_psi6..Sigma_psi5_psi11` and therefore does not decide $F_2(d)$.

More concretely, there exist two **explicit** toy instantiations satisfying constant‑vacuum stationarity with nonzero
$v_{\psi k}$, one with $F_2(d)=0$ and one with $F_2(d)\neq 0$.

### Witness A (stationary, $F_2(d)=0$)

Take any symmetric translation-invariant kernel mixing channel (so the symmetrized mixing row‑dot is the same at each
site on a constant vacuum), and choose:

1. a constant nonzero vacuum: $v_{\psi k}\equiv v\neq 0$,
2. site‑uniform local coefficient families: $g4_{\psi k}\equiv g4$, $g6_{\psi k}\equiv g6$ (and any $(gY_k)$, $v_\phi$),
3. define $m2_{\psi k}$ per-site from the stationarity equation (possible because the EoM is affine in $m2_{\psi k}$).

Then the resulting diagonal profile is constant: $d_k\equiv d$, hence:

$$
F_2(d)=d\sum_{k=0}^{11}e^{i\pi k/3}=0.
$$

So stationarity is compatible with $F_2(d)=0$.

### Witness B (stationary, $F_2(d)\neq 0$)

Keep the same symmetric translation-invariant mixing channel and the same constant nonzero vacuum $v_{\psi k}\equiv v$,
but introduce one local defect in the quartic family, e.g.:

$$
g4_{\psi 0}=g4+\delta,\qquad g4_{\psi k}=g4\ \text{for }k\neq 0,
\qquad \delta\neq 0,
$$

with the rest of the families unchanged (take $g6_{\psi k}$ uniform for simplicity).

Again define $m2_{\psi k}$ per-site from stationarity.

Then (by the Yukawa-free diagonal rewrite of `N474`) the diagonal profile inherits a nonzero single-site defect:

$$
d_0-d_k = 2\,\delta\,v^2\neq 0\qquad (k\neq 0),
$$

so its mode‑2 coefficient is nonzero:

$$
F_2(d)=\sum_{k=0}^{11} d_k\,e^{i\pi k/3}
     = (d_0-d)\,e^{i\cdot 0}
     = 2\,\delta\,v^2 \neq 0.
$$

So stationarity is compatible with $F_2(d)\neq 0$ as well.

### Conclusion

Since both $F_2(d)=0$ and $F_2(d)\neq 0$ remain compatible even under constant‑vacuum stationarity and $v_{\psi k}\neq 0$
premises, no strict-derived decision of the canonical diagonal mode‑2 defect follows from stationarity alone on the
current export class.

Therefore `T166` remains open unless additional strict-derived constraints/value objects are exported which decide the
canonical diagonal residual sector (equivalently: decide the six opposite‑pair sums used by `N467/P434`).

## Audit probe

`P436` provides a toy-level executable check of the two-witness construction above (no strict-derived value claims).

## What N476 does not prove

`N476` does not prove:

1. any strict-derived canonical vacuum vector $(v_{\psi 0},\ldots,v_{\psi 11})$,
2. any strict-derived relations fixing the local coefficient families $(g4_{\psi k})$, $(g6_{\psi k})$,
3. any strict decision of $F_2(d)$ for the canonical FIN diagonal residual sector (`T166` remains open),
4. any strict-core `pair1` `O(2)` cut, strict theta export, strict-core selector closure, or `QW-2191` discharge,
5. ToE closure.

