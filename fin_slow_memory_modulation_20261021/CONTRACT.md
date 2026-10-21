# OM: the actual memory convolution on a slow test-envelope scale

Date: 2026-10-21. Conditional analytical operator lemma, no new dynamics.

Activate only the operator-estimate component of NM/NEXT_OPEN_LEMMA
steps1–2. Use its exact existing memory k=k_s+k_c, with shape frequency
eta=sqrt3 and mass d^2=27pi^2/8192, and

    k_s(t)=d^2 sin(eta t)/eta,
    k_c(t)=sin(2t+3pi/4)/(16sqrt(pi)t^(3/2))+O(t^(-5/2)).

The source family fixes omega^2=a0+2 in[5/4,7/4]. Its quadratic forcing
uses carrier frequencies in

    V={0} union [sqrt5,sqrt7] union [-sqrt7,-sqrt5].

For S>0, 0<epsilon<=1 and any given complex A in C2([0,S]), study

    C_(nu,epsilon)A(t)=integral_0^t k(u)e^(i nu u)A(epsilon^2(t-u))du,
    0<=t<=S epsilon^(-2), nu in V.

These are analysis test envelopes. No assertion that the actual full-PDE
amplitude already has a uniform C2 slow-envelope representation is added.
The derivative norm below is in original model time t, not slow time.
Test drives characterize the existing response operator; none is injected
as a new physical force or incoming field in the actual RP preparation.

Target: a uniform O(epsilon^2) error with an explicit causal initial-slip
term, including the nondecaying shape part. The scaled drive epsilon^2 A
then has an O(epsilon^4) projected response error. This is a scalar
convolution estimate, not a global spatial H2 approximation of outgoing
radiation or a full nonlinear field theorem at time epsilon^-2.

At the actual continuum threshold nu=2 and shape resonance nu=sqrt3,
retain explicit countercontrols; they are outside V and do not retune a0.
No new packet, incoming field, bath, stochastic law, Gamma value, averaged
envelope equation or physical postulate is selected. All old files,
manifests, plans, AGENTS and previous worktree changes remain preserved.
