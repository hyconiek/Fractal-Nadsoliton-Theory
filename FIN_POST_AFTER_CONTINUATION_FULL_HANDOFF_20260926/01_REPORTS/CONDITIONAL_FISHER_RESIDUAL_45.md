# CONDITIONAL-FISHER-RESIDUAL-45
## The stationary hidden OU covariance is the conditional Fisher Schur complement

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact block-covariance theorem;
- independent numerical replay of the localized values reported in
  `EDGEWORTH-STATIONARY-15`;
- exact uniform-state coincidence;
- conditional statistical interpretation only.

---

## 1. Joint visible-hidden score covariance

Let

    X = X7
    Y = orthonormal hidden k=1,2 basis
    S = diag(p)-p p^T.

Then define

    F = X^T S X,
    H = Y^T S X,
    G = Y^T S Y.

The full score covariance is

    C_joint =
      [[F, H^T],
       [H, G]].

When F>0, define

    K = H F^{-1}

and the residual hidden coordinate

    z = y-Kx.

Then exactly

    Cov(x,z)=0

and

    boxed:
    Sigma_z = G-H F^{-1}H^T.

Thus Sigma_z is the Schur complement / conditional covariance left in the
hidden sector after linearly removing all information predictable from the
retained visible sector.

---

## 2. Uniform state

At p_j=1/12:

- different real Fourier sectors are orthogonal;
- H=0;
- Y is Euclidean orthonormal and zero-mean.

Therefore exactly

    G=I4/12,
    Sigma_z=I4/12.

So at the uniform state:
- raw hidden Fisher geometry;
- conditional hidden fluctuation geometry;
- stationary hidden OU covariance

all coincide.

This explains why the distinction is invisible in the symmetric phase.

---

## 3. Localized coexistence state

At the certified localized root center, the raw hidden Fisher eigenvalues are

    eig(G) =
      0.009484130475
      0.013087691423
      0.021996619612
      0.084895006590.

The visible-hidden coupling has

    ||H||_2 = 0.111059154035.

After conditionalization,

    Sigma_z =
    [[ 0.006007281152,  0,              -0.003252045067,  0              ],
     [ 0,               0.015812116764,  0,              -0.006981163288],
     [-0.003252045067,  0,               0.012123631839,  0              ],
     [ 0,              -0.006981163288,  0,               0.008067488259]].

Its eigenvalues are

    0.003956603996
    0.004601351485
    0.013529561507
    0.019923001027.

These reproduce the localized hidden-covariance spectrum already reported in
`EDGEWORTH-STATIONARY-15`.

The condition number is approximately

    5.035379089.

---

## 4. Information interpretation

The difference

    G-Sigma_z = H F^{-1}H^T

is positive semidefinite.

Therefore localization does not merely make hidden fluctuations anisotropic.
A substantial part of the raw hidden score variance becomes linearly
predictable from retained fluctuations.

The Schur complement is the irreducible Gaussian hidden uncertainty relative
to the retained seven-dimensional description.

---

## 5. Exact relation to the stationary OU residual

The admitted leading stationary heat-bath generator is

    L_z
      = -z.grad_z + Sigma_z:Hess_z.

Equivalently,

    dz = -z dt + sqrt(2 Sigma_z) dW

in the declared heat-bath clock.

Its stationary covariance is exactly Sigma_z.

Therefore the same Schur-complement object has two roles inside that declared
model:

1. conditional statistical covariance;
2. stationary covariance of the decoupled hidden OU residual.

This is a genuine bridge between Fisher/statistical geometry and the admitted
heat-bath dynamics.

---

## 6. Whitening theorem

Since Sigma_z>0 at the localized state, define

    q = Sigma_z^{-1/2} z.

Then the stationary hidden OU becomes

    dq = -q dt + sqrt(2) dW_q.

Thus all four hidden residual directions have:
- identical dimensionless decay rate 1;
- isotropic stationary covariance I4

in whitened coordinates.

The localized anisotropy resides in the map between physical/model hidden
coordinates and the whitened residual coordinates, not in four different OU
rates.

This is exact within the declared stationary heat-bath model.

---

## 7. What this does to the kinetic-ratio problem

Inside the specific heat-bath lane, the four hidden residual modes do NOT need
four independent kinetic coefficients: after conditionalization they share one
OU rate.

However this does not solve the broader source problem because:

- the heat-bath update law was supplied;
- its clock normalization was supplied;
- the retained phase/visible sector has a different drift matrix;
- inertial and tree-memory models remain equally compatible with the same
  statics unless independently excluded.

So conditional Fisher geometry explains the hidden OU normalization **after a
kinetic category is chosen**, rather than selecting that category.

---

## 8. Important distinction

Raw Fisher:

    G

answers:

    "How much do hidden score directions fluctuate in the full categorical
    state?"

Conditional Fisher residual:

    Sigma_z

answers:

    "How much hidden fluctuation remains after the retained seven-dimensional
    information is optimally removed at Gaussian order?"

For the physical-emergence programme, Sigma_z is the more relevant hidden
object whenever visible information is already part of the coarse state.

---

## 9. Consequence for the moving-phase memory result

Reports 40-41 projected hidden modes into a phase-locked scalar.

The correct stationary Gaussian amplitudes for that exercise, around a
localized equilibrium, should be calibrated from Sigma_z rather than the raw
G block.

The representation-induced memory kernel shape remains fixed by the moving
phase connection, but its stochastic amplitudes/cross-correlations inherit
Sigma_z.

This gives a cleaner separation:

    memory-kernel SHAPE
      <- representation / phase motion

    stochastic AMPLITUDE
      <- conditional hidden covariance Sigma_z.

---

## 10. Disposition

`CONDITIONAL-FISHER-RESIDUAL-45` closes as:

    EXACT_STATISTICAL-DYNAMIC_MATCH
    INSIDE_DECLARED_STATIONARY_HEAT_BATH_MODEL.

This is stronger than Fisher-44 but still conditional on the supplied heat-bath
kinetics.

---

## 11. Next research atom

### FISHER-HEAT-BATH-NATURALITY-46

Ask whether the heat-bath generator is distinguished among reversible local
Markov generators by requiring simultaneously:

1. stationary multinomial/exponential-family law;
2. D12 covariance;
3. conditional Fisher whitening of hidden residuals;
4. local single-label refresh/composition consistency;
5. no extra state-dependent mobility tensor.

Acceptance:
- uniqueness up to one overall clock rate; or
- explicit alternative reversible local generator satisfying the same
  requirements but giving different visible/phase kinetics.

This directly tests whether the remaining kinetic choice can be reduced from
"postulated heat bath" to a naturality theorem.
