# FISHER-KINETIC-METRIC-44
## Fisher geometry is canonical but state-dependent and does not supply a universal kinetic metric

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Inputs used:
- `fin_rank7_followup/src/model.py`;
- `R7P-026_equal_energy_event.json`;
- `EDGEWORTH-STATIONARY-15.md`.

Status:
- exact Fisher-pullback formulas inside the supplied exponential family;
- numerical evaluation at the interval-certified localized coexistence center;
- exact symmetry statements;
- NO claim that Fisher geometry is a physical kinetic metric.

---

## 1. Fisher object

For

    p(theta)=softmax(X7 theta),

the Fisher information in a logit direction f is

    <f,g>_F
      = f^T S g,

where

    S=diag(p)-p p^T.

Equivalently, for natural-parameter directions collected in a matrix Z,

    G_F = Z^T S Z.

This is exactly the covariance/Fisher object already present in the FIN
exponential-family analysis.

---

## 2. Certified localized state

Use the center of the interval-certified equal-energy localized root:

    s3 = 1.8199035812800827
    s4 = 1.9139895546687241
    s5 = 1.9145691325468481
    s6 = 1.3672032801955039

at

    g = 3.7183448981203875.

The corresponding reflection-even seven-coordinate dual point is

    theta=(s3,0,s4,0,s5,0,s6).

The probability maximum is approximately

    max_j p_j = 0.836365226651.

Thus the state is strongly localized rather than close to the uniform
distribution.

---

## 3. Phase tangent

Along the continuous effective rotation of the retained k=3,4,5 pairs, the
dual-coordinate phase tangent at the reflection-even representative is,
up to overall sign,

    t_theta =
      (0, 3 s3,
       0, 4 s4,
       0, 5 s5,
       0).

The k=6 parity column has no first-order sine partner on the 12-label carrier.

Its logit tangent is

    v_theta = X7 t_theta.

The Fisher information for the intrinsic phase coordinate is therefore

    I_theta = v_theta^T S v_theta
            = 2.644071152609.

This number is coordinate- and state-specific. It is not a physical inertia.

---

## 4. Hidden H4 Fisher block

Take the Euclidean-orthonormal real Fourier basis

    (k1c,k1s,k2c,k2s).

At the localized state the raw hidden Fisher block is

    G_H =
    [[ 0.043154256877,  0,               0.035425993047,  0              ],
     [ 0,               0.018993442333,  0,              -0.005343982533],
     [ 0.035425993047,  0,               0.054828441137,  0              ],
     [ 0,              -0.005343982533,  0,               0.012487307755]].

Its eigenvalues are

    0.009484130475
    0.013087691423
    0.021996619612
    0.084895006590

and

    cond_2(G_H) = 8.951269366.

Therefore the hidden Fisher geometry is strongly anisotropic after
localization.

At the uniform state, by Fourier orthonormality,

    G_H = I4/12

exactly.

So the localized anisotropy is state-generated.

---

## 5. Phase-hidden mixing

For the five directions

    (phase, k1c, k1s, k2c, k2s),

the Fisher matrix is

    [[ 2.644071152609,  0,              -0.005272485592,  0,              -0.053226541092],
     [ 0,               0.043154256877,  0,               0.035425993047,  0             ],
     [-0.005272485592,  0,               0.018993442333,  0,              -0.005343982533],
     [ 0,               0.035425993047,  0,               0.054828441137,  0             ],
     [-0.053226541092,  0,              -0.005343982533,  0,               0.012487307755]].

The phase-hidden cross-vector has norm

    0.053487043113.

Thus Fisher geometry does not split as

    phase direct-sum k1 direct-sum k2

at the localized state.

---

## 6. Exact reflection decomposition

The localized representative is reflection-even.

Therefore Fisher parity is exact:

- cosine hidden directions are reflection-even;
- phase tangent and sine hidden directions are reflection-odd.

The five-dimensional Fisher matrix therefore decomposes exactly as

    EVEN: (k1c,k2c), dimension 2

and

    ODD:  (phase,k1s,k2s), dimension 3.

The nonzero phase-hidden mixing is allowed entirely inside the odd block.

This is not numerical symmetry breaking.

---

## 7. D12 covariance

Translated localized minima have Fisher matrices related by the induced D12
basis action.

Therefore:
- eigenvalues;
- condition numbers;
- phase-hidden cross norm

are invariant along the 12-state orbit after the corresponding basis is
co-rotated.

The anisotropy is a property of the localized orbit, not of an arbitrary
choice of label origin.

---

## 8. State dependence along the localized stationary branch

Solving the same reflection-even stationary equations on the local branch gives:

    g        I_theta      cond(G_H)   ||cross||     max(p)
    -------------------------------------------------------
    3.65     2.87318438     8.4523     0.0606653    0.80808
    3.70     2.70240175     8.8277     0.0552674    0.82957
    3.71834  2.64407115     8.9513     0.0534870    0.83637
    3.80     2.40649320     9.4375     0.0465337    0.86167
    4.00     1.93463268    10.3615     0.0339940    0.90331
    4.50     1.15776098    11.9528     0.0167612    0.95452

Thus the Fisher ratios change systematically with localization.

There is no state-independent triple of relative coefficients corresponding to

    mu_theta, mu_1, mu_2.

---

## 9. What Fisher DOES provide

For N independent categorical samples from p, the empirical score/mode
fluctuation covariance is proportional to

    G_F/N.

So Fisher geometry is not arbitrary: it is the canonical local statistical
fluctuation geometry of the supplied exponential family.

This gives it a stronger status than an invented kinetic metric.

---

## 10. What Fisher does NOT provide

A statistical metric does not by itself choose an evolution law.

For the same Fisher geometry one may still declare:
- Euclidean gradient flow;
- natural gradient flow using G_F^{-1};
- inertial geodesic-like dynamics;
- heat-bath Markov dynamics;
- non-Markovian storage dynamics.

Therefore a successful Fisher construction cannot by itself close the
kinetic-source problem.

---

## 11. Disposition

`FISHER-KINETIC-METRIC-44` closes as:

    CANONICAL_STATE_DEPENDENT_STATISTICAL_GEOMETRY
    BUT_NOT_UNIVERSAL_KINETIC_METRIC.

It does reduce arbitrariness if a natural-gradient postulate is supplied, but
that postulate remains dynamical input.

---

## 12. Next atom

The raw hidden Fisher block is not the actual stationary hidden residual
covariance away from the uniform state because visible and hidden score
directions are correlated.

The next task is therefore:

    CONDITIONAL-FISHER-RESIDUAL-45

Compare

    G = Y^T S Y

with the stationary Schur-complement residual

    Sigma_z = G-H F^{-1}H^T,

and determine whether the latter gives a more canonical hidden metric after
conditioning on all retained information.
