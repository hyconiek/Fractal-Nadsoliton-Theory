# Fixed-population remainder certificate and its limits

Date: 2026-10-20. Same-author numerical/interval stress, not an independent
audit or a new source/process theorem.

For each frozen p define the analytical extension Phi_p(epsilon) by AD's
exact count-row expression, with a_N=a0+epsilon a1 and g=12/Lambda6.
It is evaluated only as a mathematical extension between integer N values;
no fractional copies or new microscopic law are asserted.

AD's symbolic construction gives Phi_p(0)=1 and Phi_p'(0)=0. The interval
automatic differentiation propagates coefficients (value, derivative,
half-second-derivative). Over epsilon in [0,1/264], its last interval J_p
therefore bounds Phi_p''/2. Taylor's integral remainder then gives

    [Phi_p(1/N)-1]N^2 in J_p, N>=264,               (1)

whenever the frozen p is an actual N-type. This uses the analytical zero
identities; merely containing 1 and 0 in two interval outputs would not
prove them. They are independently checked for consistency, not promoted
from containment into exact equality.

Explicitly, Phi_p(0)=tau a0^2/r and the first derivative numerator is
2tau a0a1+(1+kappa)a0^2+a0 b dot grad a0-1. The definitions of a0,a1
make these exactly 1 and 0 algebraically. The derivative enclosure tests
are a different chain-rule evaluation of those same identities.

The certificate is conditional on the inherited exact spectral boxes and
the supplied Fourier construction of A7. Every subsequent arithmetic
operation rounds outward, with no floating implementation premise. It
is uniform in N at each of FOUR p only. It is NOT a certificate uniform
over all p, bounded lambda intervals, actual MT time evolution or the
global exact scaling. It does not independently prove AD-002--004.

## Distinct formula routes used

Let h=(g/2)Ap, u_i=exp(h_i), v_i=p_i exp(-h_i), S+=sum u, S-=sum v.
The implementation uses

    tau=S+S--1, b_i=u_i S--v_i S+,
    kappa=(g/2)[A_diag S+S--v^T A u].

For h_dot=(g/2)Ab, the derivative along b is

    D_b tau=S- sum(u_i h_dot_i)
       +S+[sum(exp(-h_i)b_i)-sum(v_i h_dot_i)].

This gives D_b a0=-a0 D_b tau/(2tau) and the AD first corrector. Separate
ordered-edge sums check tau, b, kappa and D_b tau at the frozen p. The
full row then visits every LEGAL one-copy departure, including entries
at a simplex face without dividing by p_i.

## Arithmetic and adverse result

Endpoints lie on the exact grid 2^-80. Addition, multiplication, reciprocal
and square root use directed integer bounds. Exp uses a positive Taylor
sum at |x|<=1/8, a geometric remainder, range-reduction squaring and
reciprocal for negative arguments. The jet formulas are the chain rule;
interval dependence is retained conservatively, never canceled numerically.
The analytical tau bound 11 exp(-H)<=tau<=11 exp(H) narrows value
intervals safely while retaining their derivative enclosures.

The whole epsilon-box did NOT certify trial positivity at the vertex;
it is deliberately reported as inconclusive, not as proof of a negative
trial. Direct fixed-N values there are positive. The broad vertex J_p
also shows poor numerical sharpness. No point, gain, threshold or cover
was changed to obtain a positive label. The existence proof's eventual
global positivity and practical N0 are not settled by these tests.

At the user's direction work continues without delegation. The
[same-author analytical recheck](SAME_AUTHOR_ANALYTICAL_REVIEW.md) makes the
grid stopping-time term explicit. The old AD independent audit order
remains unexecuted; its independent-author criterion is not discharged.
