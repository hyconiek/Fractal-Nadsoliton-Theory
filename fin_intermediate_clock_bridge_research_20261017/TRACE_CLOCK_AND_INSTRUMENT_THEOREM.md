# A conditional reaction bridge: exact FIN trace, paid clock, and an instrument boundary

Date: 2026-10-17. Comparison model ID: **FIN-ONE-HOLE-TRACE-nu-kappa**.
This is a mathematical realization of the supplied unreflected FIN generator,
not a new adopted FIN law, natural activity source or physical vacancy.
The trace/Schur formula is standard; the native-N binding, clock payment and
observation distinctions are made explicit here.

## 1. A reversible intermediate construction for the actual full parent

For fixed N, constant-diagonal symmetric A=A7 and finite g, write

    V(n)=g n^T A n/(2N), h_b(r)=(g/N)(Ar)_b,
    Z(r)=sum_b exp(h_b(r)), s_b(r)=exp(h_b(r))/Z(r).

The visible labelled states are x in {0,...,11}^N. Add intermediate states
(i,x_minus_i): exactly slot i is temporarily empty, while every other label
is retained. At most one intermediate is active globally. From each visible
state activate slot i at rate nu; from its intermediate fill it with b at
rate kappa exp(h_b(r)). All other slots are frozen during that intermediate.
nu,kappa>0 have the same clock dimension. This serial constraint is part of
the **comparison construction**, not derived physical incidence/control.

The unnormalized visible stationary weight is exp(V(n)). For a labelled
intermediate with cavity r, choose

    w_H(i,r)=(nu/kappa) exp[g(r^TAr+A_diag)/(2N)].

On every visible/intermediate edge, exp(V(r+e_a))nu equals
w_H kappa exp(h_a(r)), so the finite augmented chain is reversible.
The stationary law **conditioned on visible states** is exactly native
Gibbs. Its unconditional law also has positive intermediate mass.

The trace process erases time in intermediates. With Q partitioned into
visible V and hidden H, Q_HH is diagonal with entries -kappa Z(r), giving

    Q_trace=Q_VV-Q_VH Q_HH^(-1)Q_HV.

Thus a tagged off-diagonal trace rate is exactly nu s_b(r); summing slots
gives the original labelled LOO heat bath at attempt rate nu. Its count
quotient has activation rate nu n_a to r=n-e_a and refill rate
kappa exp(h_b(r)) to r+e_b. Its trace is exactly the native count generator,
including all multiplicities and null refresh contributions to the diagonal.
No asymptotic approximation is needed for this **trace-clock identity**.

## 2. What finite hidden duration changes

Let lambda=N nu. Starting in a visible state, a marked reset cycle has
duration

    W=Exp(lambda)+Exp(K_i(r)), K_i(r)=kappa Z(r),

and selects b with probability s_b(r), independently of the intermediate
duration. Native marked attempts instead have waiting time Exp(lambda).
At finite kappa these are not the same physical-time record. Under a frozen
cavity, even the actual visible-change waiting law is a geometric sum of
these two-stage cycles, not generally exponential. Mean clock matching is
not equality of full path laws.

For comparison define the **last-visible-state functional**: keep the last
complete visible state during each intermediate interval, update only when
the intermediate ends. This is a specified analysis/readout functional of
history, not an instantaneous Markov projection or a physical memory device.
Actual raw instruments must not silently drop intermediate outcomes.

## 3. A full-N physical-time TV payment for that declared functional

For fixed lambda,K, convexity of TV and the shifted exponential law give

    TV(Exp(lambda),Exp(lambda)+Exp(K))
       <= E[1-exp(-lambda Exp(K))]=lambda/(lambda+K).

Finite contexts permit a measurable maximal coupling conditional on the
current visible state and chosen slot. Use the same refill mark. Before
the first failed waiting-time coupling, the whole last-visible record agrees
with native FIN. The decision to try cycle j before horizon T depends only
on previous native waits; its conditional failure probability is at most
delta=lambda/(lambda+K_min). Native marked attempts are Poisson(lambda T).
The stopped union bound, retaining the next unfinished wait, therefore proves

    TV(Law(last-visible augmented path[0,T]),Law(native path[0,T]))
       <= min{1,(lambda T+1)lambda/(lambda+K_min)}.    (1)

Both chains start from the **same arbitrary visible initial law**, with no
intermediate initially present. Marks here are analysis variables, not
licensed additional observed data. Forgetting them contracts TV.

For M_A=max|A_ab|, Z(r)>=12exp(-|g|M_A)=c_Z>0, uniformly in N and r.
Hence K_min>=kappa c_Z. On a critical horizon T=sqrt(N)S/nu, bounded g
and fixed S, a sufficient serial-realization condition for (1) to vanish is

    kappa_N/nu >> N^(5/2).

This is a deliberately conservative sufficient cost for this particular
single-global-intermediate construction and full coarsened path TV. It is
not necessary, optimal, feasible resource advice or a calibrated clock.
Other simultaneous-intermediate mechanisms need their own proof.
No late-iid mixing/preparation theorem transfers to the changed comparator
for free: charge the entire preparation horizon if that is claimed.

## 4. Raw observation does not inherit that path-TV theorem

Embed native FIN in the larger alphabet with no intermediate outcome.
For every finite kappa and visible start, the raw augmented path visits an
intermediate before T with probability 1-exp(-N nu T). Native paths never do.
Therefore raw-path TV is at least

    1-exp(-N nu T),                                  (2)

independent of kappa. Small stationary intermediate occupancy does not
remove this path event. Fast elimination only licenses a declared censored/
finite-resolution observation, not exact equivalence of ideal raw records.
Equation (1) concerns last-visible records; (2) concerns the raw alphabet.
They are compatible, not contradictory limiting claims.

## 5. Reflection and elimination must not be interchanged silently

The construction above is for the unreflected full parent. Trace first,
then reject native domain exits, if using the MT reflected analysis chain.
Deleting disallowed refills before eliminating the intermediate changes
Z(r) to Z_allowed(r), hence renormalizes accepted transitions.

Example: a three-label frozen reset with s=(1/2,1/3,1/6), nu=1, and allowed
visible labels {0,1}. Native trace-then-reject has rate0->1 equal1/3.
Suppressing label2 refills first gives2/5. This finite negative control is
not a small-N MT example; it demonstrates the order-of-operations error.
An exact reflected augmentation would need separately paid rejection/
memory structure, not an assumed commutation of trace and reflection.

This trace retains the native phase-crossing F edges. It is **not** the
H-restricted diagnostic flow with prescribed PC04 marginals and cannot
be used to declare DH-002 solved. It addresses a separate implementation/
clock question about the already declared parent.

## 6. What remains physical input

The intermediate, its rates, common activation nu, refill kappa exp(h),
initialization and readout functional are supplied comparison ingredients.
They do not derive themselves from A7 or the nadsoliton. This gives a
conditional mechanism by which the twelve-state normalizer arises from
competing exit channels, evading the earlier **direct gap-only** rate class
without contradicting its no-go. It still does not independently source
the activities, absolute time, carrier, locality or physical apparatus.
