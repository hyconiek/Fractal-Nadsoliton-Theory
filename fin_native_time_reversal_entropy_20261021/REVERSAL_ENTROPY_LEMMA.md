# Finite reversible path identities and a sharp entropy transfer bound

Date: 2026-10-21. General elementary tools reused/derived explicitly;
the native FIN application is a separate theorem.

## 1. Do not conflate the two reversed preparations

Let Q be a finite conservative CTMC generator, reversible in positive pi,
and let mu,nu be positive initial probabilities. Write P_mu,T for its
unmarked path law and mu_T=mu exp(TQ). For a path x_0,...,x_J and its
reversed holding intervals, the holding exponentials cancel. Detailed
balance telescopes the jump-rate ratio to pi(x_J)/pi(x_0). Thus

    log[dP_mu,T / d(R#P_nu,T)]
       =log(mu(x_0)/pi(x_0))-log(nu(x_T)/pi(x_T)).           (1)

All paths have common support under these preparations; no null refresh
attempts are inserted into the observed path. Fixed-time jumps have
probability zero, so the cadlag reversal convention causes no endpoint
ambiguity. Positive initial laws include paths with zero genuine jumps.

Because the likelihood in (1) is measurable in (x_0,x_T), those endpoints
retain the entire relative entropy of this forward/reversed experiment.
In particular A equals the KL between Law(X_0,X_T) and Law(X_T,X_0).
This is sufficiency for this SPECIFIC reversal comparison, not for arbitrary
targets, generators or hidden-state inference.

Expectation gives the exact identity

    D(P_mu,T || R#P_nu,T)
       =D(mu||pi)-D(mu_T||pi)+D(mu_T||nu).                 (2)

Consequently, with A=D(P_mu,T||R#P_mu,T),

    E=D(P_mu,T||R#P_mu_T,T)=D(mu||pi)-D(mu_T||pi),
    A=E+D(mu_T||mu),    0<=E<=A.                          (3)

The general boundary-preparation distinction is standard stochastic
thermodynamics, not claimed as a FIN discovery. See Seifert,
[cond-mat/0503686v1](https://arxiv.org/pdf/cond-mat/0503686v1), equations
14–18 and 22–30. We supply the finite reversible derivation above; the
paper's physical heat-bath interpretation is not imported into FIN.

## 2. Reversal-relative entropy from TV and a bounded likelihood

For any probability P and measurable involution R, put P^R=R#P. If
h=log(dP/dP^R) is bounded, then h(Rx)=-h(x) and

    D(P||P^R)=D(P^R||P)
      =(1/2) integral (dP-dP^R) h
      <= TV(P,P^R) ||h||_infinity.                        (4)

TV means half the L1 distance. This factor matters: no extra factor two
is hidden in its definition. For (1) with nu=mu, put f=log(mu/pi). Then
||h||_infinity<=osc(f). If S is an exactly reversal-invariant comparison
with TV(P,S)<=epsilon, contraction under R and the triangle inequality give

    D(P||P^R)<=2 epsilon osc(f).                          (5)

We use only TV closeness to S, NOT D(P||S). The latter can be infinite
when S is reflected and P includes rare outside states or departures.

## 3. Why the extra bound is genuinely necessary

On {0,+,-}, let R fix 0 and swap +,-, and for integer n>=2 put

    P(+)=1/n, P(-)=1/(n 2^n), P(0)=1-P(+)-P(-).

Every probability is positive, but

    TV(P,P^R)=(1-2^(-n))/n ->0,
    D(P||P^R)=(1-2^(-n)) log 2 ->log 2.

Here ||h||_infinity=n log 2 and (4) is an equality. Thus the TV warning
in the old OBS review was correct. It is not removed by a variance test or
an informal appeal to stationarity. The native theorem must pay a useful
likelihood bound AND a sufficiently fast comparison error.

An even simpler support control is P=(1-epsilon,epsilon) and S=(1,0).
TV(P,S)=epsilon, yet D(P||S)=infinity. This is why entropy is not moved
directly from the reflected local Gibbs comparison to the actual parent.

## 4. Given passive records

If a measurable record map F obeys F(Rx)=R_Y F(x), relative-entropy data
processing gives

    D(F#P || R_Y#F#P)<=D(P||P^R).                          (6)

For an integral Y(t)=integral_0^t f(X_s)ds, the correct reversal is
R_Y Y(t)=Y(T)-Y(T-t), not Y(T-t) with an unrepaired origin. This reuses
the actual OBS convention. An arbitrary causal postprocessor, invasive
instrument or postselected experiment is not automatically equivariant.
