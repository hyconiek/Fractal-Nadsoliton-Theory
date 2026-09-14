# R7P-113--120 — scoped quantum/spectral comparison

## Scope freeze

The supplied quantum materials are a separate conditional branch.  The explicit
22-product rank-132 separable stationary construction belongs to the
separable-stationarity report.  The later discord-robustness report extends the
scope to a mixed Hamiltonian family and a larger pure-flow-equivalent class;
its marginal assumptions are not interchangeable.  This package therefore
keeps canonical interaction, mixed-family, and pure-family claims separate.

## Spectral semantics and replay provenance

For `C=I/12+W/20` the exact density-gap expressions are
`lambda6/20` and `(lambda6-lambda5)/20`.  The discord theorem, however, is
stated with exact rational lower certificates

- `Delta_L=234218204114629/2000000000000000`,
- `delta_L=27234855667/12500000000000`, and
- `d0=delta_L*Delta_L/7392`.

The currently bundled rounded strict endpoints give the independent lower
quantity `(lambda6_lo-lambda5_hi)/20 = 4357576906719/2000000000000000`,
which is smaller than the historical `delta_L` by exactly
`1/2000000000000000 = 5e-16`.  Consequently the current rounded provider is
**not** promoted as an independent replay of that historical endpoint.  The
historical discord theorem is retained as accepted source-package evidence,
with this provenance gap recorded explicitly.

## No causal bridge

The classical rank-seven variational branch and the quantum branch share
spectral numbers, but the supplied sources provide no map between their state
spaces, generators, preparation resources, gain source, or operational access.
A common factor involving spectral gaps is therefore an algebraic comparison,
not a theorem that discord causes localization or vice versa.

## Controlled loading result

Let `C_gamma=I/12+gamma W` for `gamma in [49/1000,1/20]`.  Writing
`alpha=20 gamma` gives `C_gamma=(1-alpha)I/12+alpha C_0` with
`alpha in [49/50,1]`, so positivity follows by convexity from the accepted
positive endpoint state `C_0`.  For positive gamma the two relevant simple
spectral gaps scale linearly with gamma.  Thus, under the same canonical mixed
family proof hypotheses and retaining the accepted historical rational gap
certificates, the discord floor obeys

`d0(gamma) >= (49/50)^2 d0`.

This is a scoped robustness theorem for a changing supplied marginal.  It is
not a localization theorem and does not derive a gain coefficient.

## Operational resources

The separable construction uses supplied program copies, shared classical cut
choice, local projective measurements, comparison/heralding, and (where joint
properties are tested) joint measurements.  Equal local marginals do not
identify the joint preparation law, and separability does not by itself prove a
universal unknown-input LOCC implementation.

## Final O-lane classification

O contributes one new narrow loading-robustness statement and clarifies
operational nonidentifiability.  It supplies no causal spectral bridge, active
source law, selector, or physical closure.
