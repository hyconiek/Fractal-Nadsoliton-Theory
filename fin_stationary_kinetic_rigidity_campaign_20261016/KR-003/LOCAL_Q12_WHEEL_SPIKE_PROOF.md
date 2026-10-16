# An explicit localized q12/A7 wheel stress: equicontinuity can fail while NC holds

This independently reconstructs the **kind** of balanced-wheel negative control reported in the preceding activity-covariance research. It is not claimed as a new direction or a different-limit counterexample. The source file for the earlier implementation was unavailable; this proof uses the frozen count generator and an exact six-edge combinatorial motif.

## 1. Embed the rational wheel in the genuine q12 count domain

Choose four distinct label indices `a=0, b=2, c=1, d=3` of parity signs `(+,+,-,-)`. For N a large multiple of 12, write `u_N=(N/12,...,N/12)` and define the integer base `n_base=u_N−e_b−e_d` (sum N−2). For each q4 vector `s=(s_a,s_b,s_c,s_d)` with nonnegative entries totaling 2, let `n(s)=n_base+s_a e_a+s_b e_b+s_c e_c+s_d e_d`. These ten genuinely distinct count states all satisfy `n(s)−u_N=O(1)` and eventually belong to the interior of the real MT domain D_N (by its admitted central-component property). In particular `s=(0,1,0,1)` corresponds **exactly** to `u_N`.

Order the ten mobile configurations lexicographically (as in `verify_kr.py`). Their indices 0,…,9 refer to the table printed by that script, with six undirected edges:

`(0,6):+1`, `(0,3):-1`, `(1,6):-1`, `(1,4):+1`, `(3,8):-1`, `(4,8):+1`.

Every edge is an allowed single-copy replacement. For the signed edge vector `v`, direct integer calculation verifies for **each** vertex i:

`Σ_{j:{i,j} edge} v_{ij} [m(n(j))−m(n(i))]=0`.

This is a parity-DRIFT divergence-free stress. Crucially the corresponding **squared**-increment sum is nonzero at configurations 3 and 4 (with opposite signs).

## 2. Perturb native symmetric equilibrium conductances, not affinities

Let `c^0_N(x,y)=pi_D,N(x) q^0_N(x,y)=pi_D,N(y)q^0_N(y,x)>0` denote actual reflected native HB conductances for the six interior edges. Set `eps_N=delta min_{e in wheel} c^0_N(e)` with fixed sufficiently small `0<delta<1/4`, and define

`c~_N(e)=c^0_N(e)+eps_N v_e` on precisely these edges; `c~=c^0` elsewhere. Define off-diagonal rates `q~_N(x,y)=c~_N(x,y)/pi_D,N(x)` and corresponding negative row diagonals. Positivity and reversible Gibbs law are automatic. Every altered move remains a single-copy transition. The perturbation of `L m(x)` is `(eps_N/pi_D,N(x)) Σ_e v_e [m(y)−m(x)]=0` by the exact integer wheel identity, at every state including the wheel states; all other states are unchanged. All boundary jumps remain suppressed by native reflection.

At a near-uniform state native HB off-diagonal rates for any selected source/target species are `Theta(N)` since all labels have proportion ≈1/12, the field is uniformly bounded and the target softmax probabilities are bounded away from zero. Because these ten configurations differ by O(1) particles, their Gibbs weights and reversible conductances have **uniformly bounded positive pairwise ratios**. Consequently `eps_N/pi_D,N(n(i))=Theta(N)` for wheel states. The perturbation changes the local scaled bracket `G_N=N Gamma_N(m,m)` at the exact uniform state `n(3)=u_N` by a quantity of fixed nonzero order (in fact negative for the chosen orientation). The total jump rate remains ≤C N uniformly (and with smaller delta may also preserve the rate-one per-copy total attempt capacity using the available self-refresh slack near u).

## 3. Why it does NOT refute critical process universality under NC

The perturbation is supported at only ten count states whose MT-rescaled `(Z,Y)` locations approach `(0,0)`. The accepted continuous joint quartic/Gaussian equilibrium limit has no atom at that point, so

`pi_D,N({n(s):sum s=2})→0`.

The perturbation of G is uniformly bounded, zero outside those states, and its stationary absolute mean therefore vanishes. Since native HB has `E_pi|G_N^0−2|→0` by continuity and MT local concentration, this wheel perturbation also satisfies NC. Reversible generator balance moreover forces `E G~_N=E G^0_N` *exactly* because the parity drifts agree at every state.

Thus the wheel disproves automatic **statewise equicontinuity** of G around u but **does not disprove NC**, IBA, or the KR-002 stationary critical limit. It is a mandatory kill control on confusing the different regularity properties. Its construction requires no new microscopic law adopted as primary FIN.
