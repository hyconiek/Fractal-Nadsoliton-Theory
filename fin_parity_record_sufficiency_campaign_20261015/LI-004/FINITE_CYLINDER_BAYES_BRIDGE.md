# LI-004 — A bounded full-parity-history Bayes-risk bridge

**Disposition:** CONDITIONAL THEOREM PROVED BY A FINITE-CYLINDER L1 ARGUMENT, provided LI-001–003 and their explicitly inherited source assumptions are correct. This is not universal statistical equivalence, likelihood-process total variation, all-filter convergence or mutual-information convergence.

## Precise claim

At the frozen primary setting `λ=0`, one tag, observation horizon `s=1/2`, target time `s=1`, `N∈12ℕ`, iid-uniform start followed by the declared deterministic MT burn-in `t_N`, let `O_N` be the *whole unmarked tag-parity path* (or its exact integral Y_N, which generates the same finite-N sigma-field). Let `C_full,N∈{MINUS,ZERO,PLUS}` be the future collective sign. Denote the best 0–1 Bayes risk by `r_N^full`. Let `(Z,Y)` be the admitted stationary FIN quartic/observation diffusion `dZ=−Z³/3ds+√2 dB`, `dY=Z ds+√2 dW`, with B and W independent, and let `r_*^Y` be the best 0–1 risk of forecasting `sign(Z_1)` from `σ(Y_s:s<=1/2)`. Then:

`r_N^full → r_*^Y`, **conditionally on the inherited MT stationary local window, critical UF/OBS path transfer and the newly paid LI likelihood comparison**.

The result establishes an asymptotic *single-target bounded decision risk* only. All finite-N ZERO states remain in the risk; the limiting ZERO mass is zero.

## Exact full-record conditional subdensities

Under Q_N (independent fair tag and autonomous background), let `F_N^O=σ(S_t:0<=t<=√N/2)` and `C_bg,N=sign(Z_bg,N(1))`. The biased-tag product-start R measure is `dR^prod/dQ=L_N` on the background/fair-tag observation sigma-field, and the future target is background-measurable. For `c∈{−,0,+}`, define

`γ_{c,N}(O)=E_Q[1_{C_bg,N=c} L_N | F_N^O] = E_{Q_bg}[1_{C_bg,N=c} L_N(background,O)]`.

The last equality uses only the **paid product reference Q_N**. It would be invalid for the original correlated P start. The optimal biased-R risk is exactly

`r_N^{Rprod}=1−E_{Q_tag} max_c γ_{c,N}(O)`.

No posterior division is needed because the γ are unnormalized target subdensities relative to the *fair parity-path reference*. The corresponding limit uses `L=exp( (1/2)∫Z dY^0−(1/4)∫Z²ds)` and `C=sign(Z_1)`:

`γ_c(Y^0)=E_{Z}[1_{C=c} L(Z,Y^0) | Y^0]` (writing the product-reference integral), and `r_*^Y=1−E_{Y^0} max_c γ_c(Y^0)`.

The likelihood change of measure makes Y^0 under L obey the classical observation law `dY=Z ds+√2 dW` as in OBS, while the background marginal remains the stationary quartic process because `E[L|Z]=1`.

## Finite-cylinder L1 theorem (the essential nontrivial bridge)

Let `X_N=Z_bg,N|_[0,1]`, `V_N=Y_N|_[0,1/2]`. By LI-003, `(X_N,V_N,L_N,C_bg,N)⇒(X,V,L,C)`, `L_N` is UI, `(X,V)` has continuous paths and `C=sign(X_1)` has `P(C=0)=0`. For each c and ε>0 choose a *bounded continuous cylinder function* `h_c(X(t_1),...,X(t_a),V(s_1),...,V(s_b))` with fixed times, such that

`E_{Q}|1_{C=c}L−h_c(X_grid,V_grid)|<ε`.

This is legitimate: finite rational-time evaluations generate the sigma-field of the continuous-path product law; conditional expectations onto refining finite cylinders converge in L1, and each bounded measurable finite-coordinate function can be approximated in L1 by bounded continuous functions under its probability law (standard regularity/Lusin approximation). The function may use **background future coordinates only inside the proof integrand**; after the background is integrated out, the resulting rule sees only the observed V coordinates up to 1/2.

By joint weak convergence and uniform integrability of L_N,

`E_{Q_N}|1_{C_bg,N=c}L_N−h_c(X_N,grid,V_N,grid)| → E_Q|1_{C=c}L−h_c(X_grid,V_grid)| < ε`.

Here sign discontinuity is harmless because X_1 has a density. Consequently, L1 contraction of conditional expectation implies

`E_Qtag |γ_{c,N}(O)−H_{c,N}(V_N,grid)| <= ε+o(1)`, where

`H_{c,N}(v)=E_Qbg h_c(X_N,grid,v)`.

The Q product reference is crucial. For every fixed v, background weak convergence gives `H_{c,N}(v)→H_c(v)=E_Qbg h_c(X_grid,v)`. In fact this convergence is uniform on compact sets of v: background tightness confines X_grid to a compact with arbitrarily small error, and bounded continuous h_c is uniformly continuous on the resulting product compact; a finite net in v gives equicontinuity. The observed finite grid V_N is tight, hence `E|H_{c,N}(V_N,grid)−H_c(V_N,grid)|→0`.

Since H_c is bounded continuous and `V_N,grid⇒V_grid`,

`E_Qtag max_c H_c(V_N,grid) → E_Q max_c H_c(V_grid)`.

The map `(x_c)_c→max_c x_c` is 1-Lipschitz in sup norm, so replacing each γ_{c,N} by H_c costs at most `∑_c (ε+o(1))` in expectation; the same approximation holds for limiting γ_c. Let ε→0 to prove

`E_Qtag max_c γ_{c,N} → E_Q max_c γ_c` and hence `r_N^{Rprod}→r_*^Y`.

**This directly pays the feared conditional projection step**: joint weak convergence **alone** is insufficient, but joint likelihood convergence *together with UI, limit-likelihood measurability and product reference* gives the needed L1 approximation at fixed finite observed cylinders. One takes N→∞ at fixed cylinder, then refines; no uniform-in-path supremum, illicit resolution swap, or use of q12 labels/attempt marks occurs.

Finally LI-001 transfers risk from product-start R to actual unmodified P by TV through future time 1 and the small sign-target correction. Hence the stated result.

## What is *not* proved

This does not equate all full-record laws in TV, all posterior filters, KL/entropy, mutual information, full q12 tag histories or arbitrary decision problems. A change in gain/initial law/time or more than one fixed parity tag needs a new proof contract. The argument uses the accepted model's fixed mathematical clock and local Gibbs preparation, not a physical detector or temporal irreversibility.

## Audit obligation

The conditional-cylinder lemma above is mathematically explicit, but an independent intake should audit the LI-002 high-probability occupation bounds as transferred from *stationary reflected MT* and the LI-003 stopped stochastic-integral identification against the actual repository input/reviews. If any inherited MT path theorem is not independently admitted, downgrade the end-to-end claim to `CONDITIONAL` and retain the finite entropy theorem as exact. The archive's absent postcampaign review artifacts must not be invented.
