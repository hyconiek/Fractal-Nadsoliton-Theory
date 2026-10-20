# DEL-005: two scoped consequences for the limiting record

Date: 2026-10-20. These are mathematical model predictions, not empirical
apparatus performance or finite-N likelihood/entropy transfer theorems.

## 1. Bounded two-time memory contrast

Let delta_P=p'_P/p_P and compare (1) with the canonical stationary quartic
diffusion L0=partial_z^2+b partial_z. Both have marginal rho, noise sqrt2,
and E[b+delta_P|Z]=b because sum_P p_P delta_P=0. Their compressed
first generator/Dirichlet pairings agree. Their records do not.

Freeze f(z)=exp(-z^2). For the stationary two-time correlation C(t),
the same reversible-domain argument as HM/CPL gives

    C_dagger''(0)-C0''(0)
         =E_rho[ (p'(Z))^2/(p(Z)(1-p(Z))) * f'(Z)^2 ]>0.    (1)

Here f and Lf are bounded, which pays f in the L2 generator domain.
The scalar spectral Taylor formula requires only that domain, not Dom(L^2).
Thus the bounded observable f(Z0)f(Zt) has a strict positive contrast for
all sufficiently small fixed positive t. DEL-004 and the admitted primary
HB stationary critical limit transfer its expectation limit to the two
actual finite count models. No useful numerical lag or finite-N accuracy
is claimed. The generic compression identity is reused, not rediscovered.

## 2. Both limit records have the SAME quadratic variation, not TV-one noise

Unlike CPL, these two diffusions both have path quadratic variation 2t.
They are mutually absolutely continuous on every finite horizon. One exact
likelihood pays this distinction. Under the canonical path law define

    h_P(z)=sqrt(p_P(z)), V_P=(L0 h_P)/h_P,
    Lambda_P=h_P(Z_S)/h_P(Z0) exp[-integral_0^S V_P(Z_s)ds].

Ito gives dLambda_P=Lambda_P delta_P(Z)/sqrt2 dB. Since delta_P is bounded,
stopped second-moment estimates prove mean one and uniform integrability.
An Ito/Bayes calculation changes the generator to L0+delta_P partial_z.
The actual stationary mixture initial law q_P then yields the soft-record
Radon--Nikodym density

    L_S=sum_P sqrt[p_P(Z0)p_P(Z_S)]
                       exp[-integral_0^S V_P(Z_s)ds].       (2)

Each term is strictly positive; h_P is bounded above/below and V_P bounded,
so both path laws are equivalent on finite S. For an augmented reference
with the SAME initial (P,Z0) law and canonical future Z, the initial KL
cost is zero. The usual exponential-martingale identity, or direct Ito
calculation of log Lambda under the changed law, gives by contraction

    KL(Law_dagger(Z_[0:S]) || Law_0(Z_[0:S]))
       <=(S/4)E_rho[(p')^2/(p(1-p))]<=S/30.                (3)

Indeed (p')^2<=1/32 and p(1-p)>=15/64. Pinsker gives the sufficient upper
TV bound min(1,sqrt(S/60)). This is a LIMIT-LAW bound only. No path-TV or
entropy convergence from finite jump paths to continuous diffusion is
inferred, and no physical clock/detector/noise law is introduced.
