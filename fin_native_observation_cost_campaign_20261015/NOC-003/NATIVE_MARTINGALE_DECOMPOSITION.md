# NOC-003 — Dependence-preserving concentration argument

Assume the reflected chain starts in `pi_D`, write `Y_i=X_i=(b_i,o_i)` at integer native times and set `A_M=M^-1 sum_{i=0}^{M-1} g_N(Y_i)`. Let `D_i=psi_N(Y_{i+1})-P^D_1 psi_N(Y_i)`. With respect to the **full** native skeleton filtration, `(D_i)` are martingale differences. By summing the exact Poisson identity:

`sum_{i=0}^{M-1}(g_N(Y_i)-mu) = psi_N(Y_0)-psi_N(Y_M)+ sum_{i=0}^{M-1}D_i + sum_{i=0}^{M-1}e_N(Y_i)`.

The corrector's linear growth and native `pi_D(F^2+R^2)<=C` give `pi_D psi_N^2<=C`. Stationarity then yields `E(sum D_i)^2=sum E D_i^2<=M pi_D psi_N^2<=CM` and `E|psi_N(Y_0)-psi_N(Y_M)|²<=4C`. The defect satisfies `E(1/M)sum|e_N(Y_i)|<=pi_D|e_N|<=C epsilon_N,k`, without an `M` multiplier after averaging.

Chebyshev for endpoint/martingale terms and Markov for the defect, with an absolute threshold `a>0`, give

`P_piD(|A_M-mu|>a) <= C/(M a²)+C epsilon_N,k/a`, for all `M>=1` (the stronger endpoint contribution is O(1/(M²a²)) and is absorbed).

No full-parent polynomial spectral gap is used. `Y_M` is a proof-only endpoint and is **not** observed by the classifier: the statistic uses exactly `M` readings at times `0,...,M-1`. This estimates probabilities of bounded deviation events, not an unbounded moment under the original iid process.
