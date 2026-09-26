# GEOMETRY-ROOT-REORGANIZATION-18 — where the 8-phase heat-capacity peak lives

Status: **EXACT_VARIANCE_IDENTITY + FINITE-SIZE NUMERICAL PRODUCER**.

Scope: the declared 8-phase balanced hierarchy partition model.  This is not a
thermodynamic-limit theorem or a physical phase-transition claim.

## Exact identity

Write the root partition as a sum of branch contributions
`Z(alpha)=sum_r exp(l_r(alpha))`, and let `p_r` be their normalized weights.
Then

    d2 log Z / d alpha2
      = Var_p(l_r') + E_p(l_r'').

Since `alpha=n beta`,

    C_H = alpha^2 d2 log Z / d alpha2.

Thus `C_H` splits exactly into

- a **between-root** term: reweighting among root splits;
- a **within-root** term: fluctuations internal to already selected branches.

The equal-child correction is retained as an additional root branch.  Its mass
is negligible at the reported M=32/64 points but was not dropped.

## Results near each finite-size peak

| total M | alpha used | C_H total | root/between | internal/within | root fraction |
|---:|---:|---:|---:|---:|---:|
| 16 | 0.7340 | 11.09167220 | 7.29997930 | 3.79169289 | 0.65814957 |
| 32 | 0.7450 | 63.16814006 | 49.23004976 | 13.93809030 | 0.77934936 |
| 64 | 0.8245 | 300.60981483 | 260.38871756 | 40.22109727 | 0.86620165 |

Five-point derivative replays close against direct `log Z` curvature with
residuals about `1.1e-8`, `2.3e-5`, and `3.0e-7` in `C_H`, respectively.

The monotone 3-point pattern is therefore:

    65.8% -> 77.9% -> 86.6%

of the peak coming from the top-level split choice.

## Concentration of the root distribution

At the same representative peak points the effective number of ordered root
branches is approximately

- M=16: `157.74 / 555 = 28.4%` of available branches;
- M=32: `5509.55 / 38165 = 14.4%`;
- M=64: `105272.9 / 2306025 = 4.57%`.

So the absolute number of contributing branches grows, while their fraction of
the available root state space decreases sharply.

## Parent cost versus child free-energy gain

Let `E_r=-l_r'(alpha)` be the scaled conditional energy of a root branch and
`Delta_root` its bare Ward split cost.  Weighted correlations at the three
sizes are

- M=16: `corr(E_r,Delta_root)=-0.6510`, `R^2=0.4238`;
- M=32: `corr=-0.7216`, `R^2=0.5207`;
- M=64: `corr=-0.7701`, `R^2=0.5930`.

The sign is important.  More strongly segregated root splits pay a larger local
Ward cost, but their child subtrees gain even more in conditional free energy.
The peak is therefore a competition between hierarchy levels, not a simple
preference for low root cost.

At M=64 the bare root cost alone explains only about 59% of the between-branch
energy variance; about 41% remains in child-structure information after the
best linear root-cost fit.

## Boundary

This establishes a finite-size hierarchical mechanism inside the supplied
model.  It does not prove a first-order thermodynamic transition, a continuum
limit, physical geometry, `D_H=2`, or a laboratory phenomenon.
