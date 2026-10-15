# NOC-001 — Exact boundary transfer scope

The auxiliary reflected chain keeps all actual within-domain q12 count transitions and deletes only outward transitions. With `D_N` the inherited MT local domain and `B_D={n in D_N: an actual count jump leaves D_N}`, the outgoing rate from any configuration is <=N for unit-attempt HB. The Taylor increment of a `C^1` test in `(F,R)` is `<=C||Dv||/sqrt(k)`. Therefore:

`E_pi |(L_D-L)(v o Phi)| <= C (N/sqrt(k)) pi_D(B_D) ||Dv||_infty`.

Exact conditional local Gibbs weights and MT's positive boundary energy yield `pi_D(B_D) <= poly(N) e^(-bN)` for fixed `b>0`; no practical N0 is certified. The unreflected iid-start parent and reflected stationary process may be coupled for the whole required interval, with error at most

`delta_N <= mu0_N(D_N^c) + TV(nu_N P_{t_N}^D,pi_D) + P_iid(tau_D <= t_N + M_N)`.

MT bounds the last term by `mu0_N(D_N^c)+2N(t_N+M_N)(N+1)^12 exp(-bN)` for a polynomial common `t_N`; the preparation TV tends to zero with the inherited MT polynomial mixing guarantee. The **final bounded decision event**, not its unbounded samplewise mean, is transferred with this additive error. Historical MT has been used as an explicit conditional dependency, not silently reproved or altered.
