# NOC-002 — Nonzero first-order separation (amplitude-sensitive)

Let `gamma_0=1-Lambda3/Lambda6` and `gamma_1=1-(Lambda3+1/100)/Lambda6`; then `0<gamma_1<gamma_0<1`. The stationary mean is

`mu_gamma,alpha=e^(-1/2) exp[-alpha (gamma^-1-1)/2]`.

As alpha decreases to zero,

`mu_gamma0,alpha - mu_gamma1,alpha = d alpha + O(alpha^2)`,

`d = e^(-1/2)(gamma_1^-1-gamma_0^-1)/2 > 0`.

An independently reconstructed 70-digit nominal value is `d≈0.0503111362177246220735`; this decimal is **not** an independently outward-certified interval for source transcendental arithmetic. In particular there exist analytic `alpha_*>0` and `c_*=d/2>0` with the uniform separation `mu_0-mu_1>=c_*alpha` for every `0<alpha<=alpha_*`. The exact expression gives positive difference for all `alpha>0`; the restricted lower bound is what is used in NOC-004.

The comparison uses exactly the originally frozen `cos(W)` frequency=1, microscopic gain, and model-time grid spacing. This statistic is **not** invariant under arbitrary multiplicative amplitude gain, unlike the earlier DKO sign-only functional.
