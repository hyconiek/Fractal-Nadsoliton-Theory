# Next question: auxiliary representation versus a persistent mediator

Date:2026-10-16. Status: bounded analytical seed and **proposed diagnostic
research**, not execution of the next campaign or adoption of a new FIN law.
Inputs: MP7-032, native LOO HB and the NT scope. Model N,g and A=XX^T stay fixed.

## 1. The existing static joint law supplies an exact update-level bridge

For g>0, write

    Pi(n,theta) proportional to Mult(n)
        exp[-N||theta||^2/(2g)+theta.X^T n],
    theta|n ~ Normal((g/N)X^T n,(g/N)I).

This is the existing Hubbard–Stratonovich identity, not a new source law.
For one selected labelled slot, let r be the N-1 background counts, keeping
the original parent denominator N. The cavity conditional density is

    nu_r(theta) proportional to
        exp[-N||theta||^2/(2g)+theta.X^T r] sum_j exp(X_j.theta).

Drawing theta from nu_r and then the new label j from softmax(X theta)
gives exactly

    P(j|r) proportional to
        exp[g/(2N)||X^T r+X_j||^2]
      proportional to exp[(g/N)(Ar)_j],

because A_jj is constant. Thus this **ephemeral update augmentation**
reproduces the actual leave-one-out heat bath. There is no persistent-field
clock, field propagation or physical mediator in this identity. Notice the
cavity distribution is not theta conditioned on the complete old n.
The finite verifier tests the completed-square weights independently; no
Monte Carlo integral or new dynamics is needed for this algebraic seed.
For local reflection, first make the full twelve-label proposal and reject
outgoing proposals as null updates. Renormalizing only the allowed labels
would change the native reflected generator.

## 2. A natural persistent sampler is a distinct diagnostic model

For comparison only, define AUX-HB-OU-kappa on labelled counts and theta:
each copy refreshes with softmax(X theta), while

    dtheta=kappa((g/N)X^T n-theta)dt+sqrt(2 kappa g/N)dB.

This is a supplied **diagnostic** generator, not a derived replacement of
native FIN. It has the same static joint Pi (and its count-domain-reflected
version has Pi conditioned on D_N), by conditional detailed balance of each
update component. Equal equilibrium does not make the projected count law
the native LOO law. kappa is a new diagnostic parameter, not a sourced constant.

There is a sharp compatibility gate. In a reversible L2(Pi), let P be
conditional expectation onto counts and Q=I-P. For smooth finite-energy
count observables f, a projected semigroup would require

    P L^2 P-(P L P)^2=P L Q L P=(Q L P)^*(Q L P)=0.

The general positive-memory identity is reused from earlier projection work.
For parity mean m, the proposed sampler has

    L m=sigma.softmax(X theta)-m.

Conditionally on n, theta is a full-support Gaussian for every finite N,g>0.
The formula is exact for the unreflected joint law, and for interior count
states with every neighbouring type in D_N. Reflection changes the formula
at departure states and must be charged separately. A sufficiently large
balanced-uniform central type has all its neighbours in the local domain.
The first term at such an interior state is nonconstant: on the parity axis
it is a nonconstant tanh.
Hence ||Q L P m||^2>0. The next campaign must give the required generator
domains/differentiability payment and a bounded witness before admitting a
full continuous-state no-Markov theorem. A merely averaged rate or static
joint density cannot establish dynamic equivalence.

This is not a generic claim that all mediators fail. It concerns the exact
specified sampler. A controlled asymptotic approximation could still hold.

## 3. Exact marginal dynamics still does not determine hidden dynamics

An analysis-only family can retain every native count path. Let rho_n be
the Gaussian conditional above and, at each native transition n->n', refresh
theta from rho_n'. Add an arbitrary rate kappa of refreshing theta from rho_n
while n stays fixed:

    Lhat_kappa F(n,theta)
      =sum_(n'!=n) L(n,n') [integral F(n',y)rho_n'(dy)-F(n,theta)]
       +kappa[integral F(n,y)rho_n(dy)-F(n,theta)].

All members have the same joint invariant Pi, are reversible, and project
exactly to the original count generator. Their hidden histories differ.
For a centered residual coordinate v=theta_l-(g/N)(X^T n)_l,

    E[v(0)v(t)]=(g/N) exp(-kappa t)
                         sum_n pi(n)exp[-r(n)t], r(n)=-L(n,n).

A count jump erases its residual correlation; paths with no count jump
give the formula. This is a concrete specialization of known hidden-lift
nonidentifiability, not a new claim to refute all possible source theories.
It shows why preserving the full native law and static Gaussian identity
still cannot select a physical hidden clock. None of these comparison
lifts is adopted as the fundamental dynamics.

## 4. A falsifiable critical target, not an admitted SDE

The persistent diagnostic model has, on the parity line at criticality,

    m_dot=tanh(h)-m, h_dot=kappa(m-h),
    h=sqrt(Lambda6/12)theta6.

The soft right vector is (1,1), normalized left vector
(kappa/(1+kappa),1/(1+kappa)); the leading native bracket matrix for these
two coordinates is diag(2,2kappa)/N. Thus formal drift/noise projection gives

    r_kappa=kappa/(1+kappa),
    dZ=r_kappa(a Z-Z^3/3)ds+sqrt(2r_kappa)dB.

The independent symbolic fixture verifies the two-mode algebra only. A
process theorem requires the **actual joint jump-diffusion** generator,
seventeen stable directions, stationary moments, reflection/field-tail
payments and joint martingale convergence. No such theorem is admitted here.
Start fixed kappa1 and4, not a broad sweep or a new kappa_N limit.

If proved, the law would show a clock-dependent collective theory sharing
the same equilibrium. It would not select kappa, source nature's mediator,
or identify the transfer Hamiltonian with conserved physical energy.

## 5. Why this is the next physical-meaning gate

The selected question is not another spectral digit: **when does the exact
auxiliary mathematical representation remain the same stochastic theory
after giving its auxiliary variable a trajectory?** The answer distinguishes
an exact update identity, an analysis lift and a new feedback dynamics.
This is a necessary clarification before using theta as a physical field.
Independent source, locality, units and quantum/gravity gates remain unchanged.
