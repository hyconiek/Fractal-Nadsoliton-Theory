# RI-004 — Actual finite-N parity-history likelihood and sufficiency gate

**Scientific disposition: exact finite hidden-process likelihood, but `FULL_RECORD_SUFFICIENCY_OPEN`.** There is no new physical measurement law.

For one predetermined tag, set finite hidden state `X=(n,b)`, the exact Markov count/tag quotient of RI-001, and observed parity `o(X)=sigma_b` in {−,+}. Let `Q` be its finite off-diagonal generator and define block subgenerators `Q^{uv}` from hidden states observed as u to hidden states observed as v; for u=v keep all invisible transitions *and the original full Q diagonal*, so probability is lost at the rate of visible switches. On a sample history `u0 at time 0`, switch times `0<t1<...<tr<T`, observed signs u0,u1,...,ur, the exact hidden unnormalized filter is the row vector

`alpha_T = alpha_0 exp(t1 Q^{u0,u0}) Q^{u0,u1} exp((t2-t1)Q^{u1,u1}) ... Q^{u(r-1),ur} exp((T-tr)Q^{ur,ur})`.

Here the products act on the appropriate parity-state blocks; `alpha_0` is the *actual* state distribution at deterministic post-burn-in origin t_N conditioned on initial observed parity, retaining the iid-start probability of the initial parity. Its scalar total mass is the exact density of the observed switch-time history relative to the standard product counting/Lebesgue reference; summing/ranging observation paths gives the actual event probabilities. No deterministic external latent signal replaces the interacting (n,b) parent. The exact posterior of target `C_N=sign(m_N(time1))` after a history up to time1/2 is obtained by propagating `alpha_{time1/2}` with the full hidden semigroup for the remaining scaled time, and integrating the sign classes, preserving `ZERO` states. The observation times of this formula are in **model time** `sqrt(N)*s`; g_N remains fixed along the trajectory.

**Sigma-field equivalence at finite N, not in the limit.** For k=1, exact integral Y_N has right derivative `N^(1/4) S_tag(t_N+sqrt(N)s)` between jumps. Thus full continuous exact Y_N and full parity path generate the same observation sigma-field up to standard cadlag modification. In contrast the limiting diffusion Y is a.s. nowhere differentiable in its Brownian part. This illustrates why weak convergence cannot justify convergence of full-record Bayes risks.

**Precise bounded prediction comparison:** target `C_N in {−,ZERO,+}` at critical time1; observing the whole tag parity history to time1/2, compare optimal 0–1 Bayes risk `r_N^full` to the limiting continuous-Y risk `r_*^cont` for the same sign at time1 (or, more modestly, to the frozen nine-cell risk `r_*^9`). Since `R_N` is measurable from the full parity history,

`r_N^full <= r_N^9`, and `r_N^9 -> r_*^9`.

This yields only `limsup r_N^full <= r_*^9`; it is **not a lower bound** and does not imply convergence to the risk of continuous Y. In particular `r_*^cont <= r_*^9` but nothing herein proves `liminf r_N^full >= r_*^cont`.

**Exact missing lemma (`RI-LIKELIHOOD-UI`).** Under the accepted iid-start/MT local window and frozen lambda=0, prove asymptotic uniform integrability and convergence of the likelihoods (or equivalent Bayes deficiency/risk bound) of the full parity switching point-process history relative to the continuous integrated-record diffusion experiment, jointly with the future sign target. It must pay tagged-to-global feedback, hidden fast occupations, mark/timing statistics, boundary exits and order of resolution. Passing the finite 27-cell joint TV theorem or having a formal Gaussian martingale limit does *not* pay this. No independent or fitted observation noise is permitted.

**Not requested/promoted:** full q12 label history is strictly richer than parity path (tagged transitions within a parity class are invisible); the agent must not quietly replace a parity observer with that different data contract. Microscopic Y has no physically persistent storage or irreversible-time assertion. Reversing stationary records under `R Y(s)=Y(S)-Y(S-s)` preserves their admitted reversible statistics, and no entropy-production convergence is asserted.
