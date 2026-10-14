# Conditional contrast and identifiable quotient — RG-E02

## Theorem 1 (M3 mixed log-odds, model identity)
Let a *supplied* symmetric constant-diagonal q12 matrix A, M=3 labelled copies, and Gibbs exponent (g/6)n^T A n define a finite parent. Given the other two labels b,c, the heat-bath conditional distribution of the first label j is

R_bc(j)=exp[(g/3)(A_jb+A_jc)+h_j]/Z_bc.

For every j,l,b,b',c (same c),

D(j,l;b,b';c)=log[R_bc(j)/R_bc(l)]-log[R_b'c(j)/R_b'c(l)]
             =(g/3)(A_jb-A_lb-A_jb'+A_lb').

**Proof.** Take logarithms of each odds ratio; normalization cancels. On subtracting, the same c contribution and any context-independent supplied one-body field h_j-h_l cancel. The remaining alternating four-entry contrast is as shown. The result uses no attempt clock, but the *availability of conditional probabilities* requires a physical preparation/estimation protocol. End.

Consequently, a complete array of mixed contrasts determines the double-centered matrix (g/3) H A H, with H=I-11^T/q, and not the unconstrained matrix A. The double-centering is a gauge choice, not a source for g or A.

## Exact equivalences / no-go
- A -> A + c 11^T changes all Gibbs logweights by a constant; the normalized law is unchanged, even if the field is fixed to zero.
- For mixed contrasts, A -> A + 1u^T + u1^T leaves all D unchanged. In the full Gibbs model this transform adds g sum_a u_{sigma_a}, and a field h -> h-gu gives exactly the same law. Without a field uncertainty, this second transform need not leave unconditional probabilities unchanged.
- (g,A) -> (g/s,sA) for s>0 leaves the model Gibbs laws and R exactly unchanged; gain and interaction scale are not separately recoverable from normalized conditional probabilities.
- With fixed circulant A, centered log R_00 has Fourier nonzero modes proportional to Fourier modes of A. Ratios among nonzero modes cancel g and the common scale.
- A fully uncalibrated readout can have identical output distributions for every source state, destroying all observational identifiability. That is a *counterexample*, not merely insufficient signal-to-noise.

## Distinct static N2 rule
For two copies N=2 with the same Gibbs normalization, the pair difference d has P(d) proportional to exp[(g/2)A_d] when A is circulant. This rule is a separate N2 model at the PHYS-014 g=3, **not** the stationary law or dynamics of the M3 memory campaign. It has its own invisible constant shift/scale degeneracy. PHYS-014's actual proposed observable is difference-only, not M3 R_bc(j).

## Domain
All algebra is exact given the supplied finite model. Numerical Fourier ratios use the supplied exact-decimal formula in binary64 and are not derived nature-selected constants. If one-body fields vary with conditioning, the cancellation premise fails and requires a new declared nuisance model.
