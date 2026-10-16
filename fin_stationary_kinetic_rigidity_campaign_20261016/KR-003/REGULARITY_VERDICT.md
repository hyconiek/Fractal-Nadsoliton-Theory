# KR-003 — Can noise concentration NC be removed?

**Outcome of bounded attempt: OPEN_BRACKET_AVERAGING_LEMMA; do not promote unconditional universality.**

## Exact answer to the strength of the existing hypotheses

Stationarity, reversibility, single-copy support, total rate ≤CN and exact statewise scalar drift are sufficient to fix the **mean** of `G_N=N Gamma_N(m,m)` and its weak conditional projection onto the slow coordinate. In this attempt they were **not shown** to imply `E|G_N−E G_N|→0`, nor the time-integrated bracket condition that would replace NC in KR-002.

The rational ten-state control independently constructed in `verify_kr.py` shows that equal Gibbs measure and identical scalar drift do **not** fix `Gamma(m,m)` at a finite count state: a nonzero six-edge reversible conductance circulation changes two statewise brackets by up to 1/8, despite a conserved rate cap. This does not satisfy the frozen q12/A7 critical sequence and is **not** a nonuniversality limit counterexample.

A **genuine q12/A7** balanced-wheel spike is also reconstructed analytically in `LOCAL_Q12_WHEEL_SPIKE_PROOF.md`: the perturbation is supported on ten central count states, including the uniform state, retains exact Gibbs and the full statewise parity drift, changes the pointwise rescaled bracket at the uniform state by order one, and still satisfies NC because the support has vanishing stationary Gibbs mass. It is a *negative control for automatic pointwise equicontinuity*, **not** a different-diffusion example. The earlier source implementation was unavailable; this separate reconstruction relies only on the accepted q12 local-domain facts.

## New weaker law: conditional static noise rigidity

Let `b(z)=chi lambda z−z³/3`, `Z_N=N^(1/4)m` and `G_N=N Gamma~_N(m,m)` under the KR class **without NC**. For each fixed smooth compactly supported `h`:

`lim_N E_pi[h(Z_N)(G_N−2)] =0`.

**Proof:** choose a real `f ∈ C³(R)` with `f''=h`, hence `f'` bounded, `f'''` bounded and `f` at most linear. For the accelerated count generator `A_N=sqrt(N)L~_N`, Taylor's theorem, total rate ≤CN, and `|Delta Z|≤2 N^-3/4` give the uniform statewise remainder

`A_N f(Z) = f'(Z) N^(3/4)L~_N m + (1/2)f''(Z) G_N + R_N`, `|R_N|≤C_f N^-3/4`.

Stationarity entails `E_pi A_N f(Z)=0` for every finite N. The paid stationary drift error gives `E f'(Z_N) N^(3/4)L~m → ∫ f'(z)b(z)rho_lambda(dz)`. Quartic Gibbs weak convergence gives `E f''(Z_N)→∫ h rho_lambda`. The critical invariant density obeys `∫[b f'+f'']rho_lambda=0` (integration by parts; bounded f' and rapidly decaying rho). Therefore `E[h(Z_N)G_N]→2∫h rho_lambda`, which proves the claim.

This shows any weak subsequential limit of `(Z_N,G_N)` has conditional mean `E[G|Z]=2` provided boundedness and the joint limit, **not** that `G=2` a.s. It is not a temporal law: averaging over hidden fibres at one time is not the same as averaging along each trajectory. All finite-N chain variables remain correlated.

## Precisely unpaid inference

The remaining hypotheses do not specify a uniform decorrelation estimate for the bounded centred stationary observable `G_N−E_pi G_N`. One sufficient *weaker-than-NC* replacement is

`(IBA)_S: E sup_{0≤t≤S} |∫_0^t(G_N(X~_N(sqrt(N)v))−E_pi G_N)dv|→0`

for every fixed finite S. NC ⇒ IBA by stationarity, but IBA could hold even with nonconcentrated pointwise G if rapid dynamic averaging is proved. No such conclusion was obtained for all allowed kinetic activities; equally, no admitted q12/A7 family with a rigorously different limiting process was constructed. Do not treat failure to establish a sufficient condition as evidence of a changed limit.

## Conditional proof extension

Every argument in KR-002 remains valid if NC is replaced by IBA, with all other assumptions unchanged: IBA supplies the only place where NC was used, namely predictable bracket convergence. This is a **separate conditional statement**, not a derivation of IBA from the weaker class.
