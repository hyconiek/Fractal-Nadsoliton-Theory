# Exact forward/backward maximal argument

Fix N and let `T=sqrt(N)S`. For a possibly reducible finite reversible generator `L`, set `r=G-E G`, `h=P_ker(L)r`, and let `u` solve `-Lu=r-h` with mean zero on every positive-rate component. Existence is finite-dimensional self-adjoint linear algebra.

Forward martingale `M_t^+=u(X_t)-u(X_0)+int_0^t(r-h)(X_v)dv` has `E|M_T^+|²=E int_0^T Gamma(u,u)(X_v)dv=2TR`. For the stationary reversed chain, `M_v^-=u(X_(T-v))-u(X_T)+int_(T-v)^T(r-h)(X_s)ds` is a martingale with the same second moment. Simple subtraction gives

`M_t^+ + M_T^- - M_(T-t)^- = 2 int_0^t(r-h)(X_s)ds`.

Consequently `sup_t|int_0^t(r-h)|<=1/2 sup_t|M_t^+|+sup_v|M_v^-|` and

`E sup |int_0^t(r-h)|²<= (1/2)*4*2TR+2*4*2TR=20TR`.

For invariant `h`, `h(X_s)=h(X_0)` at all times almost surely, including on disconnected positive-rate graphs; its scaled integral has squared maximal expectation exactly `S²||h||²`. Using `|a+b|²<=2a²+2b²` and the `1/sqrt(N)` scaling yields the coefficient `2S²+40S R/sqrt(N)`.

All operations are legal for each finite N, not a spectral-gap assumption. The proof handles the nonergodic kernel explicitly, rather than falsely inverting `L` on component constants.
