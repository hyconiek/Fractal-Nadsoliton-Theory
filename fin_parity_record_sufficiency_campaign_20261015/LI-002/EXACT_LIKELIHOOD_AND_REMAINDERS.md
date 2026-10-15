# LI-002 — Exact fair-tag likelihood and a vanishing-error expansion

**Disposition:** exact finite-N point-process likelihood and **conditional** `o_P(1)` Taylor remainder, using explicitly identified MT stationary-path moment inputs. No fitted measurement noise, no attempt marks observed.

Under Q_N, the unchanged autonomous R background starts from the exact R background marginal, independent of `S_0∈{±1}` with `P(S_0=±1)=1/2`. The tag is a fair sign-flip chain of rate 1/2, independent of the background. The full observation is its sign history `O_N=S|_[0,T_obs]`; J counts **sign flips only**, not refresh attempts. `M_J=J−t/2` and `M_S(t)=S_t−S_0+∫_0^t S_r dr` satisfy `dM_S=−2S_(t-)dM_J`, `<M_S>_t=2t`, `<M_J>_t=t/2`.

Under the biased-tag product-start R model, given background, the sign-flip intensity is exactly `(1−S_(t-)q_N(t))/2`, `q_N(t)=σ^T softmax[(g_c/N) A n_bg(t)]`. Its likelihood relative to Q for the **observed parity path** (conditional on background) is exactly

`L_N=exp{∫_0^T log(1−S_(t-) q_N(t)) dJ_t + (1/2)∫_0^T S_t q_N(t)dt}`.

Because softmax assigns strictly positive probabilities to both parities on the finite state space, the log is finite. In particular, `E_Q[L_N|background]=1` for every background path, with *no* new initial likelihood term after LI-001's independent-fair start replacement. Initial changes are already separately charged in TV.

The exact score identity is

`−∫_0^T S_(t-) q_N(t) [dJ_t−dt/2] = (1/2)∫_0^T q_N(t)dM_S(t)`.

Define `a_N(t)=q_N(t)−m_bg(t)`. The source parity eigendirection gives the exact parity-line formula `q_N= tanh(m_bg)` when `w_bg=0` at g=g_c; uniform softmax smoothness yields

`|a_N|<=C(||w_bg||+|m_bg|³)`, `|q_N|<=C(|m_bg|+||w_bg||)`.

Here are the **sufficient dynamical controls**, to be obtained from the admitted MT-reflected stationary path plus UF/OBS transfer to P and the full-path feedback-TV transfer to autonomous R:

(A) `sup_{t<=T_obs}|q_N(t)| ->0` in Q-probability;
(B) `∫_0^T_obs q_N(t)^2 dt` is tight in Q-probability;
(C) `∫_0^T_obs (q_N(t)−m_bg(t))^2dt ->0` in Q-probability.

The MT local equilibrium stationary moment statements `E m^2=O(N^-1/2)`, `E m^6=O(N^-3/2)`, `E||w||²=O(N^-1)` imply on the **stationary reflected parent**, for `T_obs=√N/2`,

`E∫m_bg² dt=O(1)`, `E∫m_bg^6 dt=O(N^-1)`, `E∫||w_bg||²dt=O(N^-1/2)`.

Markov inequalities thus pay (B),(C) in probability on the reflected path. Under the admitted near-equilibrium initial coupling and MT exit probability `o(1)`, these events transfer to P **without multiplying TV error by T**. Full-path feedback TV transfers their probabilities to R, and the Q background has *exactly the same law* as R's background. For (A), UF/OBS's compact-horizon path tightness of `N^(1/4)m_N` and fast-tube control show `sup|m_bg|+sup||w_bg||->0`, again transferred by the same path TV. The use of probability transfer, instead of integrating an `o(1)` TV error against an unbounded factor √N, is essential.

On the high-probability event `sup |q_N|<=1/2`, Taylor gives

`log L_N = (1/2)∫_0^T q_N dM_S −(1/4)∫_0^T q_N²dt −(1/2)∫_0^T q_N² dM_J + E_3,N`,

`|E_3,N|<= C sup|q_N| ∫_0^T q_N² dJ_t`.

The stochastic quadratic remainder tends to zero in probability because its predictable quadratic variation is `(1/2)∫q_N^4dt <=(1/2)(sup q_N²)∫q_N²dt=o_P(1)`; stop at bounded bracket if necessary. The cubic remainder is `o_P(1)` by compensating dJ and the same tightness/`sup q_N` bounds; no unmarked attempt is counted as J. By (C) and localization of the martingale bracket, replacing q_N by m_bg in the linear score is `o_P(1)`. By Cauchy–Schwarz and (B),(C), replacing q_N² by m_bg² in the compensator is also `o_P(1)`.

**Conclusion conditional on the admitted dynamic controls:**

`log L_N = (1/2)∫_0^T m_bg dM_S −(1/4)∫_0^T m_bg² dt + o_Q(1)`.

The statement is a likelihood-score expansion *on the original fair parity filtration* under Q; it is not a mere finite-dimensional amplitude CLT. All integrands are predictable (left versions at possible background jumps). Q and R are auxiliary probability measures; neither modifies the original P microscopic law.
