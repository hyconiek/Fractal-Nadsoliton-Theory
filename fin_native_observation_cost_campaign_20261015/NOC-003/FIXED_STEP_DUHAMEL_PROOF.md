# NOC-003 — Exact fixed-step Duhamel comparison

Let `X_t=(b(t),o(t))` be the reflected **finite-state** CTMC, `Phi(X)=(F,R)`, `P^D_t=e^{tL_D}`, and `P_t=e^{tG}` the purely analytical OU semigroup. Denote `psi_N=psi o Phi` and `g_N=g_alpha o Phi=cos(W_N)`.

For `v=psi`, differentiation of `s -> P^D_{1-s}[(P_s v) o Phi]` gives the exact Duhamel identity

`e_N := P^D_1 psi_N-(P_1 psi)oPhi = integral_0^1 P^D_{1-s}{ [L_D((P_s psi)oPhi)-(G P_s psi)oPhi] } ds`.

This differentiation is valid on the finite native state space; the OU test has at most linear growth and derivatives through order three uniformly bounded. The chart `Phi` is **not** assumed to be Markov. Invariance and positivity of `pi_D` give `pi_D|P^D_t h| <= pi_D|h|`. The NOC-001 estimate uniformly for `P_s psi`, `0<=s<=1`, yields

`pi_D|e_N| <= C epsilon_N,k`, with `epsilon_N,k=N^-1/4+k^-1/2 + (N/sqrt(k))*pi_D(B_D)`.

Because `(I-P_1)psi=g_alpha-mu`, it follows **exactly** on each native state that

`g_N-mu = psi_N-P^D_1 psi_N+e_N`.

Note that no factor of the record length `M` appears in the *one-step* error. The final concentration requires the next martingale argument, not the unjustified statement that all dependent measurements are iid.
