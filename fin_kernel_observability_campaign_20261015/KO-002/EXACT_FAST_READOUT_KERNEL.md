# KO-002 — Exact half-population readout and sampling CLT

**Status:** exact finite hypergeometric observation kernel and conditional triangular-array finite-population CLT proof, with the **whole-population** critical quartic/fast Gaussian law inherited from MT and the comparator-specific local window conditional on KO-001.

## Readout: old physical law, different mathematical observation

N∈12ℕ, fix deterministic B_N⊂{1,...,N} of cardinality k=N/2 *before* examining states. At the common deterministic post-MT origin observe the exact q12 labels in B at one instant. This is not the one-tag parity observer of LI, nor a physical spatial region. Define `f_r(j)=sqrt(2) cos(2π rj/12)` for r=1,3,4. With complete empirical proportions p, let

`F_{r,N}=sqrt(N) Σ_j f_r(j)p_j`,
`V_{r,N}=sqrt(N) k^{-1} Σ_{a∈B} f_r(ξ_a)`,
`U_{r,N}=V_{r,N}−F_{r,N}`.

Write n_j=Np_j and H_j for count of B slots with label j. Conditional on full population counts n, labelled exchangeability gives the **exact** multivariate hypergeometric law:

`Pr(H=h | n)= [∏_{j=0}^{11} binom(n_j,h_j)]/binom(N,k)`, for integers `0≤h_j≤n_j`, `Σ h_j=k`; zero otherwise. No configurations or probability mass are discarded. This uses the fact that, given n, every assignment to slots is equally likely under the unchanged labelled heat bath started iid. Conditional exchangeability is unaffected by the auxiliary reflected *count* comparison when its exchangeable labelled lift is used.

Let `f` be the three-dimensional row vector `(f_1,f_3,f_4)`; set `μ_n=Σ_jp_j f(j)` and `Σ_n=Σ_jp_j(f(j)−μ_n)(f(j)−μ_n)^T`. Exact finite conditional moments:

`E[U_N | n]=0`, `Cov(U_N | n)= N(N−k)/(k(N−1)) Σ_n`.

The finite-population correction `N−k` and `N−1` are required. For the uniform 12-state population, `Σ_n=I_3` by Fourier orthonormality. In the limit α=k/N→1/2, `Σ_n→I_3` in probability and the coefficient tends to `(1−α)/α=1`.

## Conditional CLT, not only a covariance assertion

For **any deterministic triangular sequence** n(N)/N→u with the same q12 alphabet and k/N→α∈(0,1), the following finite-population CLT holds for the bounded f:

`U_N | n(N) ⇒ Normal_3(0,((1−α)/α) I_3)`.

Proof sketch with paid hypotheses: let `b_a=f(ξ_a)−μ_n` be the N deterministic population vectors after fixing the count, so `Σ_a b_a=0`, `max_a ||b_a||≤C` and `N^{-1}Σ_a b_ab_a^T=Σ_n→I_3`. Draw independent Bernoulli(α_N) indicators I_a with α_N=k/N and condition on `Σ I_a=k`: this conditional set of k labels is exactly the uniform without-replacement sample. Its conditional characteristic function is obtained by coefficient extraction (Fourier inversion):

`E[ exp( i t·Σ I_a b_a/√N ) | ΣI_a=k]`
`=[∫_{-π}^{π}e^{-ikθ} ∏_a (1−α_N+α_N exp(iθ+i t·b_a/√N)) dθ] / [∫_{-π}^{π} e^{-ikθ}(1−α_N+α_N e^{iθ})^N dθ]`.

Set θ=x/√N on bounded x. A second-order expansion of the log integrand is uniform because the `b_a` are bounded and α_N→α. The `x·t` cross term cancels **exactly** by `Σ b_a=0`, and the remaining Gaussian exponent has covariance α(1−α)Σ for the unscaled sampled sum `N^{-1/2}Σ I_a b_a`. The θ-integrals outside a shrinking neighbourhood are negligible by the exponential factor `|1−α_N+α_N e^{iθ}|^N`, with the t-perturbation bounded and controlled by the same product expansion. On the central region, dominated convergence and the binomial local central-limit denominator yield the conditional Gaussian characteristic function. Multiply by `1/α_N` to obtain U_N and the claimed covariance. This Cramér-Wold/Fourier argument holds for each fixed t. Its bounds are uniform along all count sequences approaching u by the subsequence principle; no k²/N small-collision approximation is used.

Consequently the conditional Gaussian noise converges **stably** relative to the complete count variables: if `X_N` is any tight statistic measurable with respect to the count n and converges jointly to X, then `(X_N,U_N)⇒(X,G)` with G∼N(0,I_3) independent of X. To see this, condition any bounded product test on n and use the deterministic-sequence conditional characteristic limit plus dominated convergence. No finite-N independence is claimed.

## Join the paid global Gibbs fast law

Under the local Gibbs `π_{ε,D,N}` (and the actual iid-start law through MT exchangeable labelled TV), the inherited joint soft/fast limit is independent: `Z_N=N^{1/4}σ·p⇒Z` with quartic density, and each `F_{r,N}=sqrt(N)Σ f_rp ⇒ N(0,1/γ_r)` with independent modes; `γ_r=1−Λ_r/Λ6` for r=3,4; `γ_1=1` because **the A7 k=1 eigenvalue is zero**. Beware that the full predecessor Laplacian has nonzero harmonic-1 eigenvalue, but A7 does not retain that mode.

To check the normalization: with Euclidean-normalized eigenvector `e_r=f_r/√12`, `F_r=√12 (√N e_r·w)`; MT fast precision on this coordinate is `12−g_c Λ_r=12γ_r`, whence `Var(F_r)=12/(12γ_r)=1/γ_r`. It would be wrong to set it equal to `(12−g_cΛ_r)^{-1}` without this factor 12.

Thus simultaneously, in the stated local stationary window:

`F_r ~ N(0,1/γ_r)`, `U_r~N(0,1)` independent across r and from F and Z, `V_r=F_r+U_r~N(0,1+1/γ_r)`.

As a null check, r=1 has `Var(V_1)=1+1=2`, equal to `1/α`. The two fast modes r=3,4 differ between ε=0,+.01 by their explicitly changed γ_r, while the parity fluctuation remains shared.

## Transfer limits and what is *not* claimed

The tagged labelled law at the common deterministic post-MT origin is exchangeable conditional on counts, and the labelled conditional-uniform orbit lift converts MT's count TV comparison to the relevant labelled snapshot law. Hence bounded event probabilities are transferred, including genuine exits and outside-domain starts. **Unbounded moments are not transferred through o(1) TV**; the claimed moment/characteristic derivations are evaluated under conditional local Gibbs, while convergence in law and bounded event probabilities pass to the actual parent.

No sample-size optimality, usable finite-N N0, physical instrument, spatial block, tag independence at finite N, new observation noise premise, or dynamic fast-path law follows.
