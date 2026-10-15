# LI-003 — Joint reference likelihood / future-target limit and unit-mass UI

**Disposition:** conditional analytical process-and-likelihood theorem. The hard *inherited* input is the admitted stationary local critical background process and the path occupation bounds paid in LI-002, not weak convergence alone.

Let `T_obs=√N/2`, `T_target=√N`, `Z_bg,N(s)=N^(1/4) m_bg(√N s)`, and `Y_N(s)=N^(-1/4)∫_0^(√N s) S_t dt`, `0<=s<=1/2`. Q_N is the independent fair-sign reference with autonomous R background. The background's law equals that in R; the finite **full-path** feedback TV bound and admitted UF/MT/OBS local window imply `Z_bg,N ⇒ Z` on `[0,1]`, where Z starts with stationary quartic density `ρ_0∝exp(-z^4/12)` and solves `dZ=−Z³/3 ds+√2 dB`. This requires the originally declared deterministic post-burn-in observation origin and the MT exit/coupling payments, not a freely substituted full global Gibbs state.

Under Q the sign process S is *independent* of the background. The exact martingale `M_S=S−S_0+∫S` has predictable bracket `2t`; `\widetilde M_S,N(s)=N^(-1/4) M_S(√N s)` converges by the bounded-jump martingale FCLT to `√2 W_s`, independent of the limiting background path Z. The exact relation

`Y_N(s)=\widetilde M_S,N(s)−N^(-1/4)(S_(√N s)−S_0)`

makes the uniform difference at most `2N^(-1/4)`, so `(Z_bg,N,Y_N)⇒(Z,√2W)` jointly.

## Joint stochastic-integral score

The LI-002 score equals `(1/2)∫_0^T m_bg(t) dM_S(t)`, equivalently `(1/2)∫_0^(1/2) Z_bg,N(s) d\widetilde M_S,N(s)`. It converges jointly with `(Z_bg,N,Y_N,Z_bg,N(1))` to `(1/2)∫_0^(1/2) Z_s dY_s^0`, where `Y^0=√2W`.

A direct way to pay stochastic-integral convergence is to fix soft radius R and partition `[0,1/2]` into finitely many mesh intervals. Replacing the predictable critical Z integrand by left-endpoint steps produces L² martingale error with bracket `2∫|Z_bg,N−Z_bg,N^step|²ds`; conditional on the independent background the martingale isometry is exact. Tightness of the background process with a continuous limit, followed by the mesh→0 limit and R→∞, makes this error `o_P(1)`. For fixed mesh, the Riemann sum is a continuous finite-coordinate function of background and tagged-martingale increments, so the standard joint FCLT applies. Cross-brackets with background jumps are zero under Q by independence. This is a proof route within the existing fair-tag clock, without imposing measurement noise.

Also, `∫_0^T m_bg(t)² dt = ∫_0^(1/2) Z_bg,N(s)² ds ⇒ ∫_0^(1/2)Z_s² ds` because critical paths converge to continuous paths. The LI-002 Taylor/occupation payment then gives

`(Z_bg,N, Y_N, C_bg,N, log L_N) ⇒ (Z,Y^0,sign(Z_1), log L)`

where `log L= (1/2)∫_0^(1/2) Z_s dY^0_s −(1/4)∫_0^(1/2) Z_s²ds` and `P(Z_1=0)=0`. The target time is 1; likelihood/observation ends at 1/2. All convergence is over a fixed scaled horizon, not over physical seconds or arbitrarily long residence times.

## Unit mass and uniform integrability are **separate** theorems

For each finite N, `L_N>=0` is the Radon–Nikodym density of biased-tag parity history conditional on background against a fair sign chain, so `E_Q L_N=1` (indeed `E_Q[L_N|background]=1`). Under Q, conditional on the limiting Z path, `Y^0` is independent Brownian motion with quadratic variation `2s`. Therefore conditional on Z,

`E[exp( (1/2)∫ Z dY^0 −(1/4)∫Z²ds ) | Z]=exp( [(1/2)²×2/2−1/4]∫Z²ds )=1`.

The conditional Gaussian stochastic integral exists because the Z path is continuous on compact `[0,1/2]` and hence has finite integral of Z². Thus **E L=1**, no global exponential moment of Z is required to show this conditional mass identity.

**Elementary UI lemma.** If `X_N>=0`, `X_N⇒X>=0`, `E X_N=E X=1`, then X_N is uniformly integrable (along the sequence). For any fixed K>0, `min(X_N,K)` is bounded continuous in X_N, so `E min(X_N,K)→E min(X,K)`; consequently `E (X_N−K)_+→E (X−K)_+`, and the latter tends to 0 as K→∞ by monotone/dominated convergence. The tails `E[X_N 1_{X_N>2K}]` are bounded by `2E(X_N−K)_+`; finitely many initial N pose no issue. Therefore `{L_N}` is uniformly integrable. This supplies the missing likelihood mass/rare-history argument without claiming convergence of entropy or all conditional filters.

## Localization tails without a fake uniform N0

For any event A measurable with respect to the autonomous background path, `E_Q[L_N 1_A]=Q_bg,N(A)`. Thus probabilities of bad background tubes under the biased measure equal those of the reference background, and escape exceptions can be propagated separately. Under the original P, use LI-001 path TV through target time. No unproved finite-N total-variation **rate for MT** or likelihood exponential moment is asserted.

**Audit boundary:** this derivation is conditional on the incoming MT/UF/OBS analytical claims genuinely satisfying their stated hypotheses. The present archive contains their producer submissions but not all independently accepted RI/OBS/MT review proof completions requested by L02–L04. It is a manuscript for independent mathematical intake, not machine certification.
