# Exact equilibrium exclusion, independent of a global saddle certificate

Date: 2026-10-14.
Status: analytical finite-model bound with exact rational arithmetic payments,
conditional on the already admitted input spectral enclosures.

Let `u=(1/12,...,1/12)` and use the unchanged all-pairs finite-N Gibbs law
relative to iid-uniform labels. For **every integer N>=1 and g in [5,6]**,

    pi_N,g(||p-u||_2<=1/10)
       <= min(1,(N+1)^12 exp(-11N/20)).

This is an equilibrium mass bound, not an error estimate for the critical
process, a mixing rate or a capacity/mean-exit theorem. Small-N bounds may be
trivial. The critical sequence g_N is in [5,6] for all sufficiently large N
for each fixed bounded lambda.

## Rational payments

The accepted active spectral intervals imply

    A_00=(Lambda3+Lambda4+Lambda5)/6+Lambda6/12 > 5/4,
    ||A||_2=Lambda6<5/2,
    (Lambda3+Lambda4+Lambda5)/Lambda6>11/4.

The verifier compares rational endpoints exactly. The positive partial sum
`sum_(k=0)^12 (5/2)^k/k!>12` proves `exp(5/2)>12`, hence `log12<5/2`.
No floating logarithm, quadrature or optimizer enters these payments.

At the exact critical gain the one-colour state therefore satisfies

    F_gc(e0)=log12-1/2-(Lambda3+Lambda4+Lambda5)/Lambda6 < -3/4.

More broadly, for g in [5,6],

    F_g(e0)=log12-g A_00/2 < 5/2-5*(5/4)/2 = -5/8.

For occupations in the specified ball, use `Au=0`, `D(p||u)>=0` and the
spectral norm bound:

    F_g(p)=D(p||u)-(g/2)(p-u)^T A(p-u)
       >= -(6/2)*(5/2)*(1/10)^2 = -3/40.

These inequalities include boundary occupations with the usual entropy
convention; interior Stirling approximations are unnecessary.

## From type weights to normalized Gibbs mass

The exact type weight relative to iid-uniform labels is

    W_N,g(n)=(N!/prod n_j!)12^(-N) exp[(Ng/2)p^T A p].

The finite-type upper bound gives `W_N,g(n)<=exp[-N F_g(p)]` and there are
at most `(N+1)^12` types. For the one-colour type, multinomial multiplicity
equals1, so **exactly** `W_N,g(N e0)=exp[-N F_g(e0)]>exp(5N/8)`.
The partition function is at least this one weight. Consequently

    pi_N,g(ball)
       <= (N+1)^12 exp(3N/40-5N/8)
       = (N+1)^12 exp(-11N/20).

The [certificate](GIBBS_EXCLUSION_CERTIFICATE.json) and
[stdlib checker](check_gibbs_certificate.py) bind the rational input and all
payments. A deliberately stronger unpaid exponent is rejected. The spectral
bridge retains its inherited assumptions; this does not derive the input A.

## Consequence for UF and for the next question

UF's finite-scaled-horizon process stays near u in occupation coordinates
while its magnified parity has an O(1) non-Gaussian law. Full finite-N
equilibrium assigns the fixed ball exponentially small mass. Thus the
full-occupation limit-order mismatch no longer needs the older coexistence
root as its only supporting witness. This says nothing by itself about **how
long** the actual symmetric-start process stays there: tiny equilibrium mass
does not prove either rapid escape or long trapping.

That missing dynamical time window, rather than another local-normal-form
constant, is the next proposed research question.
