# A model-specific physical compatibility test: zero-gap spectator activity

Date: 2026-10-17. New bounded analytical result for the already supplied
q12 native heat bath. It is not an independent activity source or a no-go
against every physical reservoir. The conditional energy interpretation
below is a comparison to thermal physics, not an adopted FIN axiom.

## 1. Exactly which reservoir class is tested

For a tagged copy with cavity counts r, native off-diagonal labelled rates are

    q_ab(r)=nu exp(h_b(r))/Z(r), h=(g/N)A7 r, a!=b.

nu is the same attempt-clock factor. If one interprets g=beta J and
E(n)=-J n^T A7 n/(2N), constant diagonal makes the conditional energy
epsilon_b(r)=-J(A7 r)_b/N up to one common constant. This identifies
dimensionless energy differences only after J,beta have been supplied.

Test the restricted class

    k_ab(r)=v(r) F_ab(epsilon_b(r)-epsilon_a(r)),      (1)

where F_ab may depend on the pair, N and fixed temperature/model parameters,
but not otherwise on r; v(r)>0 is one common context-dependent clock.
This includes population-rate formulas from a weak-coupling stationary
bath with fixed context-independent local transition matrix elements,
when their only context dependence is the Bohr gap, even allowing a
common time rescaling. Collective/context-dependent coupling operators,
memory, nonsecular mixing or a context-changing reservoir lie outside (1).

The microscopic-bath literature motivates this test; it does not assert
that every bath belongs to (1). See PLOS (2022), Section 2, equations (4)--(7),
and the carefully restricted comparison in LITERATURE_DATABASE.md.

## 2. Actual integer FIN contexts and two exactly zero gaps

Let N=12m, m>=2, sigma_j=(-1)^j and fix

    r0=m*1-e3, rt=r0+t sigma, 1<=t<=m-1.

Both are nonnegative integer cavity vectors of total N-1. Reflection of
the symmetric circulant kernel about label 3 gives

    A_03=A_63, A_13=A_53.

Also A1=0 and A sigma=Lambda6 sigma. Hence both transitions 0->6 and
1->5 have **zero conditional energy gap** in both contexts, while

    h(rt)=h(r0)+h sigma, h=g Lambda6 t/N>0.

No approximate equality of transcendental energies is used here. These
identities hold symbolically for the admitted Fourier A7 family.

## 3. A common-clock-invariant activity contrast

The full row normalizer and nu cancel between the two native rates in
one context. Therefore

    R_spec=[q_06(rt)/q_15(rt)]/[q_06(r0)/q_15(r0)]
          =exp(2g Lambda6 t/N)>1.                    (2)

For class (1), both gaps are zero and its common clock cancels, so the
same double ratio is **exactly one**. Thus no member of (1) reproduces
all these native q12 rates. Pair-dependent but fixed coupling amplitudes
are already allowed in (1); they cannot evade the double-context test.
Separate context-dependent channel activities can evade it and are
precisely extra information a provider must supply.

There is also direct normalizer dependence: putting
Ze=sum_even exp(h_j(r0)), Zo=sum_odd exp(h_j(r0)), one has

    q_06(rt)=nu exp(h_6(r0))/(Ze+exp(-2h)Zo),

strictly increasing in h. A binary heat-bath denominator is a function
of its sole energy gap; the twelve-category denominator is not. A physical
derivation of q2 Glauber rates therefore does not by itself derive FIN.

## 4. Binding to the actual local domain, when desired

At exact g_c, take N in12000N and t=N/12000. Then h=1/1000 and
`R_spec=exp(1/500)` for every N in this sequence. Adding the old tagged
label to either cavity gives actual full count types: the unperturbed
types approach u and the perturbed types approach p_theta with theta=1/1000.
Their F converges to `I(theta)-theta²/2`, strictly between zero and
theta^4/[12(1-theta²)]<h_MT. Their one-copy neighbors eventually stay
strictly inside the real MT component, so reflection does not alter these
rates. This is a domain-binding corollary, not an iid sampling theorem or
a claim that these rare conditioned contexts are easy to prepare/observe.

## 5. Physical meaning and nonconclusions

Equation (2) is an **activity-sensitive prediction of assumed native HB**,
not an independently sourced law of nature. It supplies a new mandatory
countertest for a proposed microscopic carrier: its zero-gap transitions
must retain the required spectator sensitivity. It is not the earlier
kinetics-independent mixed-affinity invariant rho36 and is not another
classifier precision result.

Both witness transitions preserve the collective parity. Matching its
drift or even proving its process rigidity does not select these activities.
Passing (2) is also not full kinetic identification: a common partial-reset
factor leaves the double ratio unchanged. Actual count transition rates
include the departure multiplicity n_a and must be divided by it before
using this labelled-rate formula. Existing partial instruments do not
automatically supply the required fully conditioned rates.

The comparison does not exclude all Davies generators: their operators
may be collective/context-sensitive, their projected dynamics may retain
memory, or their chosen model may not equal native HB. Degenerate quantum
coherences must also be controlled before a closed population equation is
used. No physical Born rule, independent A7, temperature scale, clock,
energy conservation of a closed universe or empirical validation follows.
