# DYN-001 — Causality under state-dependent local rewiring

Status: **CONDITIONAL_CAUSAL_BOUND_WITH_EXPLICIT_COLLISION_OBLIGATION**.

## Frozen inputs and scope

Repository HEAD checked: `d4c4a0ac7933ff53a958976ce6ddbde5cd514017`; scientific baseline `01c89aa3c24ce5840ef86e950530b6ae8caafd24`.
The relevant prior result is the RG-3 state-sourced circular-order graph: isolated crossings change only local cycle edges, while Haar minimum gaps scale as `M^-2`, so no positive size-uniform dwell time exists. Wave 3 additionally established that fixed-q pinning survives whereas a joint `q,N` limit can wash it out. This task does not source a physical clock or a collision law.

## Result 1 — dwell time is not needed for a metric influence bound

Let `theta_a(t)` be absolutely continuous phases on `S1`, `|dot theta_a|<=V`, and at every non-tie time let the interaction graph connect circular neighbours. Let `ell_ab(t)` be the angular length of an active edge and let a time-dependent local influence generator have nonnegative off-diagonal rates `lambda_ab(t)` only on active edges.

For `mu>0` define the **metric action budget**

`Gamma_mu = sup_(t,a) sum_(b~a) lambda_ab(t) [exp(mu ell_ab(t))-1]`.

If `Gamma_mu<infinity`, and collision updates act only inside the zero-angular-diameter collision block (no teleporting update to a positive-distance label), the time-ordered influence kernel obeys the exponential-moment estimate

`E_a exp(mu d_t(a,X_t)) <= exp[(2 mu V + Gamma_mu)(t-s)]`,

hence

`P_a[d_t(a,X_t)>=R] <= exp[-mu R + (2 mu V + Gamma_mu)(t-s)]`.

The estimate is size-uniform and uses no lower bound on inter-crossing time. Rewiring can therefore be arbitrarily frequent without, by itself, destroying an angular-metric causal cone/tail.

For a degree-two cycle with `lambda_ab <= c/ell_ab` and `ell_ab<=L`,

`Gamma_mu <= 2 c max_(0<ell<=L) (exp(mu ell)-1)/ell <= 2 c (exp(mu L)-1)/L`,

so the apparent `1/ell` singularity cancels in the **metric action** relevant to propagation distance.

## Result 2 — the same singularity still blocks an unqualified raw jump process

For a simple crossing with `ell(t)~v|t-t0|`, the raw integrated jump intensity of `lambda=c/ell` is

`integral lambda dt ~ (c/v) integral dt/|t-t0| = infinity`.

Thus `c/ell` cannot simply be interpreted as an ordinary finite-rate CTMC/update intensity through an exact crossing. A collision prescription is mathematically required: e.g. a capped rate, an instantaneous zero-distance equilibration/exchange rule, a hard-core exclusion, or a different inertial evolution. None is currently sourced by FIN, so selecting one is an added kinetic premise.

This distinction is the main result: **metric locality can survive without uniform dwell, while the raw singular generator still requires a declared collision semantics.**

## Held-out checks

Near-collision regularization `ell_eps(t)=sqrt((vt)^2+eps^2)` makes the raw integrated rate grow like `log(1/eps)`, while the weighted metric-action integral stays bounded as `eps->0`. Mixed fast/slow rewiring schedules do not alter the bound because only instantaneous support, phase velocity and `Gamma_mu` enter.

## Allowed conclusion

A size-uniform conditional influence bound exists in the intrinsic angular metric for a declared nonteleporting collision rule and finite metric-action budget. This does **not** derive a physical light cone, seconds, Lorentz symmetry, or a unique kinetic law.
