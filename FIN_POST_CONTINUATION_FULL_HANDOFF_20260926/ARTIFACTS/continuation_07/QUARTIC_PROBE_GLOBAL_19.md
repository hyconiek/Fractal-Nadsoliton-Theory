# QUARTIC-PROBE-GLOBAL-19 — global-strength search after the general quartic theorem

Status: **ANALYTIC PHASE REDUCTION + NUMERICAL GLOBAL BRANCH-AND-BOUND**.

The underlying closure coefficient is the proved conditional theorem

    C(phi) = -12 ||P_H(phi^2)||_u^2
             + 2 g <P_H(phi^2),P_H(phi A7 phi)>_u.

The optimization statements below are finite numerical certificates, not new
physical observables.

## Phase reduction

At `g_eq=3.7183448981203875`, expand the hidden k=1 and k=2 channels in the
complex retained Fourier coefficients.  Every phase-sensitive cross coefficient
in the quartic form is strictly positive.  For example the k=1 pair-channel
cross coefficients are approximately `65.58, 69.82, 79.85`; the smallest
relevant k=2 diagonal/cross margin is still about `10.19`.

Hence replacing every retained complex coefficient by its modulus cannot lower
`C`, because each `Re(t_i^* t_j) <= |t_i||t_j|`, and equality is simultaneously
attainable by aligned phases.  The global maximum therefore reduces to four
nonnegative amplitudes for sectors k=3,4,5,6.

## Coefficient-sphere normalization

For `||c||_2=1`, a stationary point is

    (0.2932238705, 0.4863537663, 0.6521040063, 0.5022351449)

with

    C = 0.0932487158813161.

The pure unit-k5 probe has

    C5 = 0.0176875757262002.

A branch-and-bound on `y_i=c_i^2`, using an exact box/simplex monomial upper
problem for each quartic monomial, exhausts the domain at tolerance `5e-5`:

    0.0932487158813161 <= Cmax < 0.0932987158813161.

Thus the best mixed probe is between

    5.27199 and 5.27482

times stronger than pure k5 under this normalization.

The branch-and-bound currently uses ordinary double arithmetic, so this is a
numerical global certificate rather than an outward-rounded interval theorem.

## Uniform-Fisher / fixed-variance normalization

A more natural statistical normalization at uniform equilibrium is

    c^T diag(lambda3,lambda4,lambda5,lambda6)c = 1,

which fixes the leading variance of the probe.  In normalized variables
`w_k=sqrt(lambda_k)c_k`, the optimum found is

    w ~= (0.3150280854, 0.4923890678, 0.6461252850, 0.4907468058),

or in the original aligned retained coefficients

    c ~= (0.2249393462, 0.3320011758, 0.4261715392, 0.3206617498).

It gives

    C = 0.0183907753866966,

whereas equally normalized pure k5 gives

    C5 = 0.00334764299804361.

The same branch-and-bound at absolute tolerance `1e-5` gives a gain interval

    5.49365 <= Cmax/C5 < 5.49664.

Therefore the mixed-probe advantage survives, and slightly strengthens, under
fixed uniform Fisher variance.

## Boundary

No apparatus or experimental normalization is sourced by FIN.  The result only
says that within the declared 7D observable family, k5 is the simplest detector,
not the strongest normalized detector.
