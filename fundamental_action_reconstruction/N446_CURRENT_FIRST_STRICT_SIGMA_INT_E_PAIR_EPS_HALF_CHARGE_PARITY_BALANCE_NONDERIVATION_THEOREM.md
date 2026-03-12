# N446 Current First Strict Sigma-Int `E_pair` eps=1/2 “Charge-Parity-Balance” Nonderivation Theorem

Status: `N446_DISCHARGED_CURRENT_FIRST_STRICT_SIGMA_INT_E_PAIR_EPS_HALF_CHARGE_PARITY_BALANCE_NONDERIVATION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

Audit, at theorem level (not probe level), one recurring proposed strict move:

```text
derive eps = 1/2 from a strict “charge parity balance / zero-charge point” principle
applied to the Z2-graded 12-slot generator weights w_{i,k}.
```

This theorem does **not** deny that a future strict-derived eps law could exist.

It proves only the strongest honest current statement:

```text
on the current exported strict sigma-int lane, no such strict derivation exists,
and the natural Z2 “balance” constraints do not uniquely select eps=1/2
(and in fact force eps=0 if imposed literally).
```

## Strict-admissible evidence reused

1. `T117`
   - weight law
     $$
       w_{i,k} := \frac{1 + \sigma_{int}^{in}\,\varepsilon\, b_{i,k}}{12},
       \qquad \varepsilon = eps \in [0,1],
     $$
     with one admissible `Z2` mask choice $b_{1,k}=(-1)^k,\ b_{2,k}=(-1)^{k+1}$.
2. `F317/N428`
   - exported eps value object
     `eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2`,
     explicitly classified as `strict_source_upgraded` (premise-based).
3. `P407/N441`
   - theta-pair depends on admissible eps choices (eps is a real selector slot).
4. `P408`
   - strict-admissibility audit: no exported strict theorem derives eps=1/2 from a “charge parity split/balance” law.
5. `T160/P410/N444`
   - strict-derived eps selection target is named and is not discharged.

## Theorem-level claims

### Claim 1. Normalization does not select eps.

For each fixed pair slot $i$ and sigma-int input $\sigma_{int}^{in}$,

$$
\sum_{k=0}^{11} w_{i,k}
  = \sum_{k=0}^{11}\frac{1 + \sigma_{int}^{in}\,eps\, b_{i,k}}{12}
  = 1 + \frac{\sigma_{int}^{in}\,eps}{12}\sum_{k=0}^{11} b_{i,k}.
$$

For the admissible strict mask $b_{1,k}=(-1)^k$ (and hence also $b_{2,k}=(-1)^{k+1}$),
there are exactly six $+1$ and six $-1$ values over $k\in\{0,\dots,11\}$, hence:

$$
\sum_{k=0}^{11} b_{i,k}=0.
$$

Therefore:

$$
\sum_{k=0}^{11} w_{i,k} = 1
$$

for **all** admissible $eps \in [0,1]$.

So any “charge balance” phrased only as “weights sum to 1” cannot uniquely select eps.

### Claim 2. Literal Z2 parity-balance constraints force eps=0 (not eps=1/2).

If one formalizes “charge parity balance / zero-charge point” as either:

1. equal total weight on the two Z2 parity classes:
   $$
   \sum_{k:\,b_{i,k}=+1} w_{i,k} \;=\; \sum_{k:\,b_{i,k}=-1} w_{i,k},
   $$
   or equivalently (since both sides sum to 1) as:
   $$
   \sum_{k:\,b_{i,k}=+1} w_{i,k} \;=\; \frac12,
   $$
2. vanishing signed first moment (a literal “zero-charge point” condition):
   $$
   \sum_{k=0}^{11} b_{i,k}\, w_{i,k} \;=\; 0,
   $$

then under the `T117` weight law and the strict admissible mask one computes:

$$
\sum_{k=0}^{11} b_{i,k}\, w_{i,k}
  = \sum_{k=0}^{11}\frac{b_{i,k} + \sigma_{int}^{in}\,eps\, b_{i,k}^2}{12}
  = \frac{1}{12}\sum_{k=0}^{11}b_{i,k}
    + \frac{\sigma_{int}^{in}\,eps}{12}\sum_{k=0}^{11}1
  = \sigma_{int}^{in}\,eps.
$$

Therefore:

$$
\sum_{k=0}^{11} b_{i,k}\, w_{i,k} = 0
\quad\Longleftrightarrow\quad
eps = 0
$$

(since $\sigma_{int}^{in}\in\{+1,-1\}$).

Equivalently, the even/odd weight sums are:

$$
\sum_{k:\,b_{i,k}=+1} w_{i,k}=\frac{1+\sigma_{int}^{in}\,eps}{2},
\qquad
\sum_{k:\,b_{i,k}=-1} w_{i,k}=\frac{1-\sigma_{int}^{in}\,eps}{2},
$$

and equality again forces $eps=0$.

So the most literal Z2 parity-balance/zero-charge constraints do **not** derive eps=1/2.
They either:

1. hold automatically for all eps (normalization), or
2. force eps=0 (and thus erase sigma-int dependence in the generator).

### Claim 3. Therefore eps=1/2 is not strict-derived from “charge parity balance” on the current lane.

From `F317/N428`, the repo exports only:

```text
eps_sigma_int_E_pair_amplitude_strict_provenance_v1 := 1/2
classification = strict_source_upgraded (premise-based)
```

and from `P407/N441` the strict sigma-int → theta candidate pipeline is eps-sensitive.

The repo exports **no** typed strict “charge parity balance” objective/law whose theorem-level unique
consequence is $eps=1/2$. (`P408`).

Hence, on the current repo state:

```text
eps = 1/2 is NOT strict-derived.
```

## What N446 does not prove

`N446` does not prove:

1. that no future strict-derived eps law can exist,
2. strict-core theta export,
3. strict-core selector closure or `QW-2191` discharge,
4. discharge of post-`T148` object-support targets (`N302/N395/T130`),
5. ToE closure.

## Consequence (next honest step)

If one wants to eliminate the eps selector slot in strict core, the next honest move is **not**
to re-label “charge parity balance” rhetoric as a strict derivation of eps=1/2.

It must be either:

1. export a genuinely `strict_derived` eps selection law/value object satisfying `T160`, or
2. change the sigma-int → theta construction class so the eps slot is absent by design,

while keeping `QW-2191` discipline explicit.

