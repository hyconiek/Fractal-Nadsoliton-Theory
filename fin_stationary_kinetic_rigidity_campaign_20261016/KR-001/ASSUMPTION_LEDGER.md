# KR-001 — The frozen diagnostic class and premises

## Fixed objects

* q=12; supplied, unmodified rank-7 A7, exact `g_c=12/Lambda6`, `g_N=g_c+lambda N^{-1/2}` with fixed bounded lambda; N a positive multiple of 12.
* True local MT count domain `D_N`, positive interior, reflected native leave-one-out HB count generator `L^0_N`. Reflection rejects only a jump leaving D_N, retaining all original twelve target probabilities. Count invariant `pi_D,N(n) ∝ (N!/∏ n_i!) exp[(g_N/2N)n^T A7 n] 1_{D_N}`.
* `m(n)=N^{-1} Σ_j sigma_j n_j`, sigma_j=(-1)^j, and `Z_N(s)=N^{1/4}m(X_N(√N s))`. A single-copy jump changes m by 0 or ±2/N.
* `L~_N`: another **count** CTMC on the identical D_N with only permitted single-copy off-diagonal moves. It is reversible with respect to the **same** pi_D,N, has total count-change rate `≤ C N` at every state with C uniform in N, and satisfies `L~_N m(n)=L^0_N m(n)` **for every n∈D_N**, including boundary states.
* Stationary initial preparation exactly `pi_D,N`; no diagnostic arbitrary start, no original iid-start transfer.
* Bracket rate `Gamma~_N(m,m)=L~_N m²−2m L~_N m = Σ_{n'} q~_N(n,n') [m(n')−m(n)]²`, `G_N=N Gamma~_N(m,m)`, and **NC** is `E_pi|G_N−E_pi G_N|→0`. This is an independently stated sufficient analysis assumption, not derived from the other class axioms.

## Reused local theorems — NOT new KR claims

UF-002 provides an exact native scalar drift expansion, after fast-mode selection, in an interior fixed cylinder. MT-003 provides quartic-soft and Gaussian-fast local-Gibbs weak convergence and all fixed polynomial moments. MT local geometry pays exponentially small invariant boundary mass. Their asserted validity is an explicit precondition to KR-002. The extracted versions were read from the original preserved UF and MT handoff ZIPs, not from a current git checkout.

## What is not inferred

Neither a fixed drift nor a stationary Gibbs measure identifies pointwise noise. `G_N` is bounded by `4C` due to jump/rate caps. Reversibility *does* determine its stationary expectation from drift, and under NC this determines its time-integrated limit. Generator pointwise equicontinuity and fast mixing are **not** class assumptions.

**Provenance gap:** The newly described `fin_activity_covariance_research_20261016/THEOREMS.md`, `RESULTS.json` and review completions and the current checkout of `AGENTS.md` are not present as original local files. The user's supplied research-order précis was used for their reported findings; this packet does not claim to have rechecked their original hashes. The visible AGENTS.md available locally is an archival copy only.
