# STRICT-KINETIC-SPLITTING-49
## Minimal strict Laplacian coupling splits an S12 reset clock into six D12 sector rates

Date: 2026-09-26

Status:
- exact conditional Markov theorem on the uniform/pre-localized strict carrier;
- exact sector-rate formula;
- one nontrivial dimensionless kinetic coupling survives;
- not a kinetic source theorem for the localized phase.

---

## 1. Ingredients

Let

    P_C = I - (1/12)11^T

be the centered projector.

The fully isotropic reset generator is

    Q0 = -gamma P_C.

Let A be the strict symmetric C12 Laplacian,

    A1=0,

with positive nonzero Fourier eigenvalues

    lambda_1,...,lambda_6.

The smallest real linear operator algebra using one strict kinetic correction is

    span{P_C,A}.

---

## 2. Minimal generator

Take

    boxed:
    Q(gamma,eta)
      = -gamma P_C - eta A,

with

    gamma>0,
    eta>=0.

Because A_ij=-W_ij for i!=j,

    Q_ij = gamma/12 + eta W_ij >0

for every i!=j.

Row sums vanish because both P_C and A annihilate constants.

Q is symmetric, so the uniform distribution is stationary and detailed balance
holds.

Therefore Q is a valid irreducible continuous-time Markov generator for every

    gamma>0, eta>=0.

---

## 3. Exact Fourier-sector rates

P_C and A commute and are simultaneously diagonal in the cyclic Fourier basis.

On the constant mode:

    rate_0=0.

On sector k=1,...,6:

    Q|_k = -(gamma+eta lambda_k).

Thus the six relaxation rates are

    boxed:
    r_k = gamma+eta lambda_k.

The strict spectrum directly orders the kinetic splitting.

---

## 4. Strict numerical rates

The frozen strict spectrum is approximately

    lambda1 = 0.754121154207
    lambda2 = 1.577049514428
    lambda3 = 1.961406861976
    lambda4 = 2.199568849333
    lambda5 = 2.298606272079
    lambda6 = 2.342182041146.

For eta>0:

    r1 < r2 < r3 < r4 < r5 < r6.

So the hidden k=1 sector is slowest and the parity k=6 sector fastest in this
minimal strict-side kinetic model.

This ordering is conditional on the sign eta>=0 and the chosen linear grammar.

---

## 5. Clock quotient

Rescaling time can set gamma=1.

The only remaining parameter is

    rho = eta/gamma.

Then

    r_k/gamma = 1+rho lambda_k.

So the minimal strict correction reduces the six arbitrary D12 rates to one
dimensionless kinetic coupling rho.

This is a large reduction compared with unrestricted D12 kinetics, but it is
not parameter-free.

---

## 6. Relation to full heat bath

At

    eta=0,

all six rates coincide and the S12 heat-bath/reset law is recovered.

Any eta>0 breaks kinetic S12 isotropy down to the strict D12 spectral pattern.

Thus the strict operator itself supplies a natural *shape* for kinetic symmetry
breaking if it is allowed to enter the generator.

---

## 7. Why rho is not currently sourced

Statics determines A but does not say that the kinetic generator must contain
A, nor with what coefficient relative to reset.

Both

    eta=0

and

    eta>0

are compatible with the same static strict operator being present elsewhere in
the theory.

Therefore the dimensionless ratio rho is a genuinely additional kinetic
choice unless a microscopic update/variational principle relates it to the
static strict coupling.

---

## 8. Localized-state boundary

Q(gamma,eta) has the uniform stationary distribution because it is symmetric.

It does not by itself generate the certified localized Gibbs/exponential-family
state.

To describe localization one would need:
- state-dependent detailed-balance rates;
- coupling to the FIN potential/order parameter;
- or a larger joint dynamics.

So this theorem belongs to the pre-localized/uniform strict side.

It must not be promoted to the kinetic law at coexistence.

---

## 9. New structural picture

There is now a simple candidate hierarchy:

    pre-strict S12 reset
      Q0=-gamma P_C

        plus strict differentiation A

    ->
    D12 kinetic splitting
      Q=-gamma P_C-eta A

        plus nonlinear/localization dynamics

    ->
    state-dependent Fisher / hidden memory.

Every arrow still needs a typed source law, but the number of free kinetic
parameters can be tracked exactly.

---

## 10. Next atom

### DETAILED-BALANCE-STRICT-LOCALIZATION-50

Ask whether there is a minimal state-dependent deformation of

    Q=-gamma P_C-eta A

that has the supplied FIN exponential-family equilibrium

    p(theta)=softmax(X7 theta)

and preserves:
- D12 covariance;
- detailed balance;
- the strict spectral splitting near uniform;
- one-label locality.

Acceptance:
- a canonical deformation with only gamma and eta;
- or a no-go showing additional state-dependent kinetic freedom is unavoidable.

This is the direct bridge from pre-localized kinetic symmetry to the actual
localized FIN state.
