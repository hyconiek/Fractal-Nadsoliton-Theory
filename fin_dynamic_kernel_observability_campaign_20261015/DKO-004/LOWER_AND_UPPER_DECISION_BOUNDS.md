# DKO-004 — Decision-theoretic lower and upper bounds

Define E_i=1{W(ih)W((i+1)h)>0} in the stationary Gaussian model and A_M=M^{-1}Σ_{i=0}^{M−1} E_i. The pair sign-agreement probability is

p_γ(α)=1/2+(1/π)arcsin(ρ_γ(α)),
ρ_γ(α)=[(1−α)q+(α/γ)q^γ]/[(1−α)+α/γ], q=e^(-h).

Expand at α=0:

ρ_γ(α)=q+α (q^γ−q)/γ+O(α²),
Δp(α)=p_a(α)−p_b(α)=κ α+O(α²),
κ=[(q^a−q)/a−(q^b−q)/b]/[π√(1−q²)] ≠0.

Nonzero κ follows from the strict monotonicity of (q^γ−q)/γ for 0<γ<1: its derivative is [γq^γ log q−(q^γ−q)]/γ²<0. For α small, |Δp|≥|κ|α/2.

**Dependence paid:** the latent vector (F,U) is a 2D stationary OU Markov diffusion with nonzero spectral gap min(a,b,1). For a bounded function of the path segment [ih,(i+1)h], its covariance with a second such segment separated by (k−1)h is bounded by exp[−γ_min(k−1)h] using the explicit Gaussian OU semigroup L2 contraction (Hermite polynomial expansion; see below). For adjacent segments |Cov|≤1 trivially. Thus

Var_γ(A_M) ≤ C(a,b,h)/M,
C=3+2e^(−γ_min h)/(1−e^(−γ_min h)),

uniformly for 0<α≤α0. Reason: after variance normalization, the independent OU coordinates each contract by e^(-their γ t); every centered square-integrable function admits orthogonal Hermite modes of degree≥1, all contracting at least by e^(-γ_min t). Conditioning on the Markov endpoints yields the stated segment covariance bound. This is not an independence assumption for E_i.

Select the fixed ideal equal-prior midpoint threshold (p_a(α)+p_b(α))/2 for A_M. Chebyshev's inequality then bounds the sum of type-I and type-II errors by a constant times 1/(Mα²), so the ideal Bayes error tends to zero if Mα²→∞ (within the Gaussian experiment). No claim of optimality of the simple sign test is needed.

**Combined resource theorem:** Mα²→0 is impossible for even full Gaussian amplitude records; Mα²→∞ suffices for bounded sign-only records, for the *fixed two-model Gaussian family* and fixed h. At Mα²=O(1) no universal sharp risk constant is asserted here.
