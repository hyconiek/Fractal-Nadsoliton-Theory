# FISHER-HEAT-BATH-NATURALITY-46
## Hidden whitening plus reversibility does not uniquely select the heat-bath generator

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact local/stationary counterfamily at any interior supplied state;
- numerical positivity interval at the certified localized coexistence state;
- exact reversibility, stationarity and visible/hidden spectral split;
- this is a local stationary naturality no-go, not yet a full nonlinear
  finite-N alternative to every detail of the declared heat-bath model.

This continues `CONDITIONAL-FISHER-RESIDUAL-45`.

---

## 1. Fixed interior state and function space

Fix any interior categorical state

    p_i>0, sum_i p_i=1.

Let L2(p) be the twelve-dimensional real function space with inner product

    <f,g>_p = sum_i p_i f_i g_i.

The constants form one dimension.

The centered space has dimension eleven.

The FIN retained feature space supplies a canonical seven-dimensional centered
subspace V.

Its p-orthogonal complement inside the centered space is a four-dimensional
hidden residual subspace H.

Thus

    L2(p) = constants direct-sum V direct-sum H.

This is the function-space version of the Fisher/Schur decomposition.

---

## 2. Canonical projectors

Let

    P_C = I - 1 p^T

be the projector that subtracts the p-mean.

Let X_c be the retained feature matrix with each column p-centered:

    X_c = X - 1 (p^T X).

Let

    F = X_c^T diag(p) X_c.

Then the p-orthogonal retained projector is

    P_V =
      X_c F^{-1} X_c^T diag(p).

Define

    P_H = P_C - P_V.

They satisfy

    P_V^2=P_V,
    P_H^2=P_H,
    P_V P_H=0,

and both are self-adjoint in L2(p).

---

## 3. Heat-bath/reset generator

The complete-reset generator at stationary law p is

    Q_1 = -P_C
        = 1 p^T - I.

For i != j,

    (Q_1)_ij = p_j >0.

It has:
- eigenvalue 0 on constants;
- eigenvalue -1 on all eleven centered directions.

This is the isotropic label-space heat-bath/reset law.

---

## 4. Exact counterfamily

For any positive constant c define

    boxed:
    Q_c = -P_H - c P_V
        = -P_C -(c-1)P_V.

Then exactly:

### Constants
    Q_c 1 = 0.

### Hidden residuals
    Q_c h = -h       for h in H.

### Retained directions
    Q_c v = -c v     for v in V.

### Reversibility
Because P_V and P_H are p-self-adjoint,

    diag(p) Q_c = Q_c^T diag(p).

Therefore p is stationary.

Thus all four whitened hidden residual directions retain exactly the same
dimensionless relaxation rate 1 as in the reset model, while the seven retained
directions have rate c.

For c != 1 this is not a global clock rescaling.

---

## 5. Markov positivity around c=1

At c=1 every off-diagonal rate is strictly positive because p is interior.

Q_c depends continuously on c.

Therefore an open interval around c=1 remains a valid continuous-time Markov
generator with nonnegative off-diagonal rates.

At the certified localized coexistence state, explicit evaluation gives

    0.978904765460 < c < 1.023692672064

as the maximal interval obtained from all off-diagonal positivity inequalities.

For example c=1.01 is a valid reversible generator.

So the counterfamily is not merely infinitesimal.

---

## 6. D12 covariance

The construction uses only:
- the state p;
- the retained feature subspace;
- the p-inner product.

Under a D12 permutation P,

    p -> Pp,
    X -> PX,

and the projectors transform as

    P_V -> P P_V P^T,
    P_H -> P P_H P^T.

Hence

    Q_c(Pp)=P Q_c(p) P^T.

The same constant c works on every translated localized state.

Thus the counterfamily does not break the internal relabeling symmetry by hand.

---

## 7. Composition / independent-label interpretation

At a fixed p, Q_c is an ordinary one-label continuous-time Markov generator.

For N independent labels, the product generator is the sum of the N one-label
generators and has product stationary law p^{tensor N}; the empirical
composition has the corresponding multinomial stationary distribution.

Therefore the counterfamily is compatible with ordinary independent-label
composition at a fixed interior state.

Important boundary:
the declared FIN heat-bath model has a self-consistent state-dependent target
q(mu). This report proves nonuniqueness of the local stationary generator. It
does not claim that Q_c has already been lifted to an exact finite-N nonlinear
self-consistent law with every original heat-bath property.

---

## 8. Consequence

The conditions

- stationary categorical/exponential-family law;
- reversibility;
- D12 covariance;
- canonical retained/hidden decomposition;
- hidden conditional-Fisher whitening;
- one-label Markov composition

do NOT determine the visible relaxation rate relative to the hidden rate.

At least one dimensionless parameter

    c = retained rate / hidden rate

survives.

Thus `FISHER-HEAT-BATH-NATURALITY-46` gives a counterexample to uniqueness up
to a single overall clock scale at the stationary linearized level.

---

## 9. What extra premise WOULD select heat bath?

If one additionally requires the generator to act as the same scalar on the
entire eleven-dimensional centered space,

    Q f = -gamma f
    for every p-centered f,

then necessarily

    Q = gamma(1 p^T-I).

So full tangent isotropy would select the heat-bath/reset generator up to one
clock rate.

But full tangent isotropy is stronger than D12 symmetry and stronger than
hidden whitening.

Whether FIN sources it is a separate question.

---

## 10. Disposition

    HIDDEN_WHITENING + REVERSIBILITY
    DOES_NOT_SELECT_FULL_HEAT_BATH.

The remaining source issue is now explicit:

    why should retained and hidden sectors share the same microscopic
    relaxation rate?

---

## 11. Next atom

`FULL-TANGENT-ISOTROPY-47`:

1. prove the uniqueness statement rigorously;
2. test whether D12/refinement/Fisher structure can imply full centered-space
   isotropy;
3. if not, classify the symmetry-allowed rate parameters at uniform and
   localized states.
