# OBS-C02 — Finite hypergeometric readout and critical sampling theorem

**Scientific status:** exact finite-N conditional kernel, asymptotic central Gaussian channel and bounded sign risk under MT local equilibrium; no entropy limit.

## Exact fibre kernel

Given population parity-plus count M=N(1+m)/2 and a fixed previously selected block B of size k, H=# {a in B:S_a=+1} has

 P(H=h|M)=C(M,h) C(N-M,k-h)/C(N,k),
 max(0,k−(N−M))<=h<=min(k,M).

All other h have exactly zero probability. E[m_B|m]=m, Var(m_B|m)=(1-m²)(N−k)/(k(N−1)). With disjoint predetermined B_1,B_2, the full joint law is the multivariate hypergeometric formula

 P(H1=h1,H2=h2|M)=C(M,h1)C(N−M,k1−h1)/C(N,k1) × C(M−h1,h2)C(N−M−(k1−h1),k2−h2)/C(N−k1,k2).

Its conditional cross covariance is Cov(m_B1,m_B2|m)=−(1−m²)/(N−1). Block labels are sampled **without replacement** conditional on counts; they are not independent even given m for finite N.

## Critical channel and scalings

Under H(pi_D,N), Z_N=N^{1/4}m -> Z~rho_lambda. For k_N/sqrt(N)->alpha∈(0,∞), k_N/N->0, conditional finite-population CLT with deterministic bounded ±1 labels and m_N=O_P(N^{-1/4}) gives

 N^{1/4}(m_B−m_N) => alpha^{-1/2}G,
 (Z_N, N^{1/4}m_B) => (Z,Z+alpha^{-1/2}G),

with G standard normal independent of Z, because the conditional Gaussian characteristic function converges uniformly on |Z_N|<=K; the quartic tail deals with |Z_N|>K. One proof uses the martingale difference CLT for successive draws without replacement: maximal standardized increment O(k^{-1/2}); conditional variance sum is (1-m²)(1−k/N)+o(1), uniformly on the central event. The same argument for two disjoint blocks with sizes alpha1/2 sqrt N gives independent G1,G2 **conditionally in the limit** (the finite cross covariance scales as O(N^{-1/2})). Unconditionally two limit readouts share Z and are correlated by Var(Z), not by spatial adjacency or entanglement.

For k_N=o(sqrt N), the conditional entire k-block sign sequence becomes independent of Z in total variation: compare hypergeometric sampling to iid Bernoulli((1+m)/2) with error <=k(k−1)/(2N); and compare that Bernoulli product to unbiased signs by Pinsker/KL O(sqrt(k)*|m|) on the central event. These errors vanish for k=o(sqrt N) and m=O_P(N^{-1/4}). Thus bounded snapshot tests cannot reliably recover Z. This includes fixed k (without falsely applying Gaussian CLT). For sqrt N<<k_N<<N, E[N^{1/2}(m_B−m)^2|m] ≤C sqrt N/k ->0 and so N^{1/4}m_B−Z_N->0 in probability. No claim for k of order N is needed.

## Frozen decision

At lambda=0, rho_0(z)∝exp(−z⁴/12) is even and continuous. For the limiting Gaussian channel, posterior odds of sign(Z) have the same sign as observation y, because for every z>0 and y>0 phi_alpha(y−z)>phi_alpha(y+z). Therefore the Bayes threshold is sign(Y), with ties at Y=0 resolved as INCONCLUSIVE (probability zero in the continuous limit). The error of this declared sign rule is

 R(alpha)=E_{rho_0}[Phi(−sqrt(alpha)*|Z|)].

The finite count event m=0 and block ties m_B=0 are retained as separate outcomes, never silently assigned a desired sign. The limit has neither zero event. Numerical values from scipy quadrature are point diagnostics, not an interval proof of risk or physical predictive performance.

Weak convergence of a mixture does not establish convergence of Shannon/differential entropy or mutual information: these claims remain OPEN.
