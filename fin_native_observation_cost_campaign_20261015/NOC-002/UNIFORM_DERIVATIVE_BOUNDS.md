# NOC-002 — Uniform series bounds, including Duhamel functions

Let `gamma_min=min(gamma_A7,gamma_control)>0`. At all `0<alpha<=1/2`, write `a_t=e^(-gamma t)`, `b_t=e^(-t)`. The Gaussian amplitude prefactor differs from its stationary limit by at most `C[alpha e^(-2gamma t)/gamma+(1-alpha)e^(-2t)]`; the phase difference is bounded by `sqrt(alpha) a_t|f|+b_t|r|` because cosine is 1-Lipschitz. Hence

`|P_t g_alpha(f,r)-mu| <= C[exp(-2gamma_min t)+exp(-2t)+exp(-gamma_min t)|f|+exp(-t)|r|]`.

Summing over nonnegative integer t gives the stated uniform **linear growth** `|psi(f,r)|<=C(1+|f|+|r|)`, independent of alpha and of the two frozen gamma values.

For each multi-index beta=(i,j), `1<=i+j<=3`, derivatives hit the cosine phase and contribute factors `(sqrt(alpha)e^(-gamma t))^i(e^(-t))^j`; the positive amplitude prefactor is <=1. Consequently

`sup_{f,r}|D^beta P_t g_alpha(f,r)| <= alpha^(i/2) e^[-(i gamma +j)t] <= exp[-gamma_min t]`.

The geometric series yields `sum_{1<=|beta|<=3} ||D^beta psi||_infty<=C`. For `0<=s<=1`, `P_s psi=sum_{ell>=0}(P_{s+ell}g-mu)` by the OU semigroup, and the same uniform derivative bounds hold. These are exactly the function classes needed in NOC-001 and Duhamel; no uniform bound of the **zeroth** derivative is claimed.
