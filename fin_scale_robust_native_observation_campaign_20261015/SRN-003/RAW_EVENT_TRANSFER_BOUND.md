# SRN-003 — Bounded raw-event concentration with explicitly paid nonsmoothness

Let X_i be the *full* reflected native (b,o) count skeleton at i=0,...,M, stationary with law pi_D. H_i=H(L_alpha phi(X_i)), g_eta=H_eta o L_alpha, psi_eta as above, and e_N=P^D_1 psi_eta(phi)-(P_1 psi_eta)(phi).

Poisson/Duhamel gives the exact full-state identity
 g_eta(phi)-mu_eta=psi_eta(phi)-P^D_1 psi_eta(phi)+e_N.
Summing i=0,...,M-1 produces two corrector endpoints, a martingale difference sum of M terms, and a sum of e_N(X_i). The martingale differences are orthogonal in L2, and the stationary second-moment bound of psi_eta gives E(sum martingale)^2<=C M eta^-2, while the endpoints have E|endpoint difference|²<=C eta^-2. Thus Chebyshev controls mean smooth fluctuation at scale alpha by C eta^-2/(M alpha²), with the endpoints absorbed when M>=1. Stationary L1 Duhamel plus Markov controls average defect by C eta^-3 epsilon_N,k/alpha, **without multiplying by M**.

By the proved native conditional angular anti-concentration, E_piD|H-g_eta|<=C(eta+k^-1/2). Under the nondegenerate independent Gaussian reference, E_G|H-H_eta|<=C eta from Gaussian strip probabilities. Markov controls average raw-minus-smooth by C(eta+k^-1/2)/alpha; replacing mu_eta with raw ideal mean p_gamma(alpha) costs Ceta/alpha. Consequently, for any sufficiently small alpha and all eta<=eta0 (constants fixed independently of N,k,M),

 P_piD(|M^-1 sum_{i=0}^{M-1} H_i - p_gamma(alpha)|>c alpha)
 <= C { eta^-2/(M alpha²)+eta^-3 epsilon_N,k/alpha + (eta+k^-1/2)/alpha },

for any fixed c>0 with corresponding constant C_c. This transfers *raw event probabilities*, not an unproved nonsmooth generator comparison. Boundary term epsilon retains the positive MT departure-layer charge. The proof only uses smooth comparison to estimate the actual raw statistic; it never changes the observer or the parent CTMC.
