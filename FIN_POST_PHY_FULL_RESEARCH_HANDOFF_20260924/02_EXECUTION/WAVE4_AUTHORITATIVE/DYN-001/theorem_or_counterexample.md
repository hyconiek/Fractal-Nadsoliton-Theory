# DYN-001 theorem / obstruction note

## Theorem (moving-metric exponential tail)

Let `X_t` be a time-inhomogeneous local jump/influence process on labels whose current positions are phases `theta_a(t)` with `|dot theta_a|<=V`. Suppose active transitions are only between circular neighbours and, for some `mu>0`,

`sup_(t,a) sum_b lambda_ab(t)(exp(mu ell_ab(t))-1) <= Gamma_mu < infinity`.

Assume collision updates have zero angular diameter. Then, between collision times, applying the generator to `F_t(b)=exp(mu d_t(a,b))` and using `|partial_t d_t(a,b)|<=2V` gives

`(partial_t+L_t)F_t <= (2 mu V+Gamma_mu)F_t`.

The same inequality passes through zero-diameter collision updates. Gronwall therefore gives

`sum_b U_(s,t)(a,b) exp(mu d_t(a,b)) <= exp[(2 mu V+Gamma_mu)(t-s)]`.

Markov's inequality yields the claimed tail bound. The proof is pathwise in the time-dependent graph and does not count rewiring events.

## Obstruction (raw 1/ell intensity at a crossing)

If a crossing has `ell(t)=v|t-t0|+o(|t-t0|)` and one sets an ordinary jump rate `lambda(t)=c/ell(t)`, then every punctured neighbourhood has diverging compensator:

`integral_(t0-eps)^(t0+eps) lambda(t) dt = infinity`.

Therefore the usual nonexplosive finite-rate jump construction does not cross `t0` without an extra rule. This is not cured by the finite metric-action bound because infinitely many zero-distance transfers can carry negligible angular distance.

## Scope

The theorem is a conditional locality result. The obstruction prevents promotion of the singular conductance itself to a complete microscopic update law. A collision rule or alternative inertial dynamics must be specified before DYN-002 can treat a concrete kinetic model.
