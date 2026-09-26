# DETAILED-BALANCE-STRICT-LOCALIZATION-50
## Canonical half-density deformation exists, but detailed balance leaves an infinite even kinetic freedom

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact classification theorem for pairwise one-label detailed-balance
  deformations of the minimal strict generator;
- exact first-order universality theorem near uniformity;
- numerical replay at the certified localized coexistence state;
- no unique physical update law is selected.

## 1. Uniform strict-side kinetic base

For i != j define the positive symmetric base rate

    b_ij = gamma/12 + eta W_ij,

with gamma>0, eta>=0.  At the uniform distribution this gives

    Q0 = -gamma P_C - eta A.

## 2. Target localized equilibrium

Let

    p_i = exp(h_i)/Z,   h=X7 theta.

Detailed balance requires

    p_i q_ij = p_j q_ji,

hence

    q_ij/q_ji = exp(h_j-h_i).

## 3. General ratio-local theorem

Assume

    q_ij = b_ij f(h_j-h_i).

Then detailed balance is equivalent to

    f(x)/f(-x)=exp(x).

Every positive solution is uniquely

    f(x)=exp(x/2) psi(x),

where psi is positive and even.

Therefore

    q_ij =
      b_ij exp[(h_j-h_i)/2] psi(h_j-h_i)

is the full ratio-local detailed-balance class.

## 4. Canonical half-density member

psi=1 gives

    q_ij=b_ij sqrt(p_j/p_i),

and

    p_i q_ij=b_ij sqrt(p_i p_j)

is manifestly symmetric.

It introduces no additional state-dependent mobility function and reduces
exactly to Q0 at uniform p.

## 5. Infinite nonuniqueness

Examples with the same equilibrium and same uniform base:

- half-density:
      psi=1;

- normalized Barker:
      psi=1/cosh(x/2),
      f=2/(1+exp(-x));

- normalized Metropolis:
      psi=exp(-|x|/2),
      f=min(1,exp(x));

- smooth families:
      psi_a=exp(a x^2).

Hence detailed balance + positivity + D12 covariance do not select a unique
localized kinetics.

## 6. First-order universality near uniformity

For every differentiable even psi with psi(0)=1,

    psi'(0)=0,

so

    f(x)=1+x/2+O(x^2).

Therefore every smooth member has the same first state-dependent correction:

    q_ij =
      b_ij[1+(h_j-h_i)/2+O((Delta h)^2)].

Kinetic ambiguity begins at second order in localization amplitude.

If

    psi(x)=1+a x^2+O(x^4),

then

    f(x)=1+x/2+(1/8+a)x^2+O(x^3).

Detailed balance fixes the coefficient 1/2 but not a.

## 7. D12 covariance

Using one common scalar psi on every pair preserves D12 covariance because
h and W co-transform under label permutations.

## 8. Localized replay

At the certified coexistence center the state is strongly nonuniform.
Half-density, Barker and Metropolis generators all have the same p and all
reduce to Q0 at uniformity, yet their nonzero relaxation spectra differ
substantially.

Thus the freedom is operational.

## 9. Disposition

    CANONICAL_MINIMAL_DEFORMATION_EXISTS
    BUT_IS_NOT_UNIQUE_FROM_DETAILED_BALANCE.

The missing premise is now exactly:

    why should the even acceptance factor psi be constant?
