# N483 Current First Strict `T169` Power‑Law Element‑Order Constrained Lift Existence + Uniqueness Theorem (No False‑PASS)

Status: `N483_DISCHARGED_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`T169` demands a **strict** constrained lift/value‑provider chain from strict scalar vacuum closure (`QW-2122`) into the
canonical per‑site arrays required by `T168`:

```text
vpsi[0..11], g4[0..11], g6[0..11]
```

so that `P437` can compute the six opposite‑pair sums (`T167`) and `P434` can evaluate the diagonal/local `pair1`
mode‑2 defect (`T166`).

`F447` exports one explicit constrained lift of exactly this type, but the downstream harness (`P437`) flags:

```text
PASS_COMPUTED_FROM_INPUTS_REQUIRES_PROVENANCE_REVIEW
```

This theorem provides the missing theorem‑level support for the `F447` lift:

1. it states the constrained lift rule as an explicit strict construction,
2. it proves existence and uniqueness of the resulting per‑site arrays (up to explicitly fixed residuals),
3. it records the remaining “non‑forced” choices as explicit **selector ingredients** (not hidden slots),
4. therefore it supplies the provenance anchor required for promoting the `F447` provider from “computable” to “strict
   value provider exported” (still without any global closure claim).

## Strict-admissible inputs reused

1. `QW-2122`
   - scalar vacuum closure outputs:
     - `rho_star_sq`,
     - `lambda_psi_strict`.
2. `F309/N420`
   - strict-derived Shannon amplitude:
     - `alpha_geo_strict_derived_v1 := 4 ln 2`.
3. `F329`
   - typed `Z_12` carrier language.
4. `N479`
   - `ord_Z12(x)` is `Aut(Z_12)`‑invariant ⇒ no marked generator/direction slot for `f(ord_Z12(x))` reference data.
5. `R14/R15`
   - frozen strict kernel specialization `K_total` and diagonal floor `m0^2 I` for the host‑matching route.
6. `T169`
   - acceptance criteria (existence + uniqueness, no hidden slots).

## Setup (reference datum + constrained lift definition)

Let `Z_12 = {0,1,...,11}` with the typed structure from `F329`.

Define the element order:

$$
\operatorname{ord}_{Z_{12}}(0)=1,\qquad
\operatorname{ord}_{Z_{12}}(x)=\frac{12}{\gcd(x,12)}\quad(x\neq 0).
$$

Fix the strict scalar inputs:

$$
\rho_*^2 := \texttt{rho\_star\_sq} \gt 0,
\qquad
\lambda_\psi := \texttt{lambda\_psi\_strict} \gt 0,
\qquad
\alpha_{\mathrm{geo}} := 4\ln 2 \gt 0.
$$

### 1) Power‑law element‑order reference distribution (no marked direction)

Define weights:

$$
w(x) := \operatorname{ord}_{Z_{12}}(x)^{-\alpha_{\mathrm{geo}}} \gt 0,
$$

and normalize:

$$
r_{\mathrm{ordpow}}(x) := \frac{w(x)}{\sum_{y\in Z_{12}} w(y)}.
$$

By `N479`, this reference datum is `Aut(Z_12)`‑invariant (depends only on `ord_Z12(x)`), hence it introduces **no marked
generator/direction slot**. It is explicitly **not** translation‑invariant (it distinguishes the identity orbit
`{0}`); this is an explicit symmetry‑breaking datum class and must remain visible.

### 2) Magnitude lift (scalar norm)

Define per‑site magnitudes:

$$
|v_{\psi i}| := \sqrt{\rho_*^2\,r_{\mathrm{ordpow}}(i)}.
$$

Then the scalar norm constraint holds identically:

$$
\sum_{i=0}^{11} |v_{\psi i}|^2 = \rho_*^2.
$$

### 3) Quartic/sextic lift (uniform minimal matching along the lifted vacuum ray)

Define sextics:

$$
g6_{\psi i} := 0.
$$

Define a uniform quartic:

$$
g4_{\psi i} := g4_{\mathrm{uniform}}
            := \frac{\lambda_\psi}{\sum_{i=0}^{11} r_{\mathrm{ordpow}}(i)^2}.
$$

This matches the scalar quartic coefficient along the lifted vacuum ray:
if $q_i := r_{\mathrm{ordpow}}(i)$ and $v_i := |v_{\psi i}|$, then
$v_i^2 = \rho_*^2 q_i$ and the induced quartic along the normalized direction has coefficient $\lambda_\psi$ by
construction.

### 4) Deterministic sign selection (finite search; explicit residual fixing)

Let `K_total[i,j]` be the strict frozen host kernel matrix from `R14`.

Consider sign vectors $s\in\{\pm 1\}^{12}$ and define the mixing energy:

$$
E_{\mathrm{mix}}(s)
:=
\sum_{i\lt j} K_{\mathrm{total}}[i,j]\,(s_i|v_{\psi i}|)\,(s_j|v_{\psi j}|).
$$

Define the minimizer set:

$$
S_{\min} := \operatorname*{argmin}_{s\in\{\pm 1\}^{12}} E_{\mathrm{mix}}(s).
$$

Fix the unavoidable global `Z2` ambiguity deterministically by requiring:

$$
s_0 = +1.
$$

If multiple minimizers remain after this constraint, choose the lexicographically smallest sign vector (with respect to
the fixed index order `0..11` of the typed carrier).

Finally define:

$$
v_{\psi i} := s_i\,|v_{\psi i}|.
$$

## Theorem (existence + uniqueness of the lift output)

The constrained lift above defines a unique per‑site triple:

$$
(v_\psi, g4_\psi, g6_\psi)\in\mathbb{R}^{12}\times\mathbb{R}^{12}\times\mathbb{R}^{12}
$$

as a function of the declared strict inputs.

### Claim 1. `r_ordpow` exists and is unique.

All weights $w(x)=\operatorname{ord}(x)^{-\alpha_{\mathrm{geo}}}$ are strictly positive and finite, so the
normalization constant $\sum_y w(y)$ is strictly positive and finite. Therefore `r_ordpow` is well-defined and unique.
∎

### Claim 2. The magnitudes `|vpsi_i|` exist and are unique.

Since $\rho_*^2>0$ and $r_{\mathrm{ordpow}}(i)>0$, each magnitude is a uniquely defined positive real number.
∎

### Claim 3. The uniform `g4` and zero `g6` arrays exist and are unique.

Since $\sum_i r_{\mathrm{ordpow}}(i)^2>0$, the uniform value $g4_{\mathrm{uniform}}$ is a uniquely defined finite real.
The definition $g6_{\psi i}=0$ is unique by construction. ∎

### Claim 4. The sign selection outputs a unique sign vector.

The finite set $\{\pm 1\}^{12}$ has $2^{12}=4096$ elements, so the minimum of $E_{\mathrm{mix}}$ exists, hence
$S_{\min}\neq\varnothing$.

Applying the explicit residual-fixing rule $s_0=+1$ and then the explicit lexicographic tie-break produces exactly one
sign vector. Therefore the final signed vacuum vector $v_\psi$ is uniquely defined. ∎

### Conclusion

Combining Claims 1–4, the `F447` constrained lift defines a unique per‑site value object `(vpsi,g4,g6)` from strict
inputs, with all residual choices explicitly fixed by the exported rules. Therefore it satisfies the *existence +
uniqueness* requirement of `T169` (with no hidden selector slot).

## What N483 does not claim

`N483` does not claim:

1. global `QW-2191` discharge (this only supports a per‑site value provider needed by the `T166` diagonal lane),
2. that the constrained lift is the only conceivable mapping on the current export class (it is one explicit exported
   strict selection ingredient; alternatives must remain explicit, not smuggled),
3. ToE closure.

