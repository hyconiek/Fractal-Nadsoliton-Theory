# LI-001 — Reference initial law, boundaries, target, and full-horizon transfer

**Disposition: conditional analytic result on admitted MT local window and labelled exchangeability.** Incoming `MT-003` gives uniform conditional-Gibbs moments, and `MT-002/004` plus the OBS tagged lift give the post-burn-in path transfer. These are *inputs*, not reproved/verified as repository intakes here.

Fix `N∈12ℕ`, `g=g_c`, `T_1=sqrt(N)`, `T_obs=sqrt(N)/2`. Write `n_bg=∑_{i=1}^{N-1}e_{xi_i}`; all processes keep the original denominator N. Under the labelled local Gibbs lift of `π_D,N`, up to a factor independent of b,

`π_D(n_bg,b=j) ∝ ( (N-1)!/∏ n_bg,i! ) exp[(g_c/(2N))(n_bg+e_j)^T A(n_bg+e_j)] 1{n_bg+e_j∈D_N}`.

For any background with all 12 completions in D_N, the constant diagonal of A cancels in the conditional distribution, yielding the **exact** conditional

`π_D(b=j|n_bg)=softmax_j[(g_c/N) A n_bg]`.

If a completion is forbidden while the actual completion is admitted, the actual count type is one legal tagged count move from the complement of D_N: the MT departure boundary. By the supplied local barrier and finite-count-type bounds, this exceptional mass `β_N=o(1)` (the MT estimates give exponentially small times a polynomial under the proposed local construction). Do not identify it with zero.

Project b to sign `S=σ_b`. At good backgrounds its conditional sign mean is `q_N=σ^T softmax[(g_c/N)A n_bg]`. Put `α_N=(N-1)/N`, `p_bg=n_bg/N`, `m_bg=σ^T p_bg`, `w_bg=p_bg−α_N u−m_bg σ/12`. Because `Au=0`, `Aσ=Λ_6 σ`, `g_c Λ_6/12=1`, and softmax is uniformly Lipschitz for bounded A,g,

`|q_N| <= C(|m_bg|+||w_bg||)` and `|q_N−m_bg|<=C(||w_bg||+|m_bg|³)`.

By `MT-003`, `E_{π_D}|m_bg|=O(N^(-1/4))` and `E_{π_D}||w_bg||=O(N^(-1/2))`, including the one-tag 1/N shifts. Hence under the local Gibbs lift

`TV(Law(n_bg,S_0), Law(n_bg)⊗Fair{±1}) <= C N^(-1/4) + β_N`.

The iid-start actual parent at deterministic post-MT origin is within inherited `ε_MT,N=o(1)` in *labelled* TV of the local-Gibbs lift, by the accepted exchangeability transfer. Therefore `R_N` with the actual initial joint law is TV-close to a product-start parity process `R_N^prod` with the **same exact R background marginal** and an independent fair sign. The initial-payment error may be taken as `δ_init,N <= 2ε_MT,N+C N^(-1/4)+β_N` (a conservative use of triangle inequalities). The full q12 tag after feedback deletion may be discarded because R-background rates no longer see it, and R-tag parity flip rate depends only on sign and n_bg.

Under this product reference, the background process is autonomous with refresh `softmax[(g_c/N) A(n_bg−e_i)]`; the tag is biased by its background, and the fair independent parity chain Q below is a **change-of-measure reference**, not a FIN law. Do NOT relabel the autonomous background as an unconditionally prepared N−1 Gibbs chain or replace denominator N by N−1.

For the future target `C_full,N=sign(N^(1/4)[m_bg(T_1)+S(T_1)/N])` (three outcomes including ZERO) and `C_bg,N=sign(Z_bg,N(1))`, the difference event obeys

`{C_full,N≠C_bg,N} ⊆ {|Z_bg,N(1)|<=N^(-3/4)}`.

By inherited critical weak convergence and its continuous limiting density, the RHS has probability `o(1)` under R, hence under P up to path TV. This is **not** a finite-N numeric rate and ZERO outcomes are not dropped at finite N.

**Full experiment risk comparison:** the optimal bounded 0–1 prediction risk of the actual P tag-parity path to time 1/2 and actual target at time 1 differs from its product-start R/background-target analogue by at most

`ε_feedback,N + δ_init,N + δ_target,N`, where `ε_feedback,N = |g_c| M_A sqrt((N-1)sqrt(N))/(2N)`, and all three terms vanish. This is an explicit triangle of TV/change-of-target comparisons, not conditioning on a successful run. Exit errors remain in `ε_MT,N` and the admitted path comparisons.
