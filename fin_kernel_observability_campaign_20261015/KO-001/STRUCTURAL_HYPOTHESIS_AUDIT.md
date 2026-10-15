# KO-001 — Matched operators, structural hypotheses and local-window transfer

**Date:** 2026-10-15. **Status:** exact rational matched-family algebra; a conditional *local-theorem transfer proof* dependent on the admitted MT/UF/OBS/LI mathematical results and the analytic MT canonical-path argument, not independent peer review. This is a new comparator audit, not a new fundamental FIN law.

## 1. Freeze the two models

Let `q=12`, `sigma_j=(-1)^j`, `u_j=1/12`, and let A7 be the previously declared circulant rank-seven `XX^T` operator with positive eigenvalues `Lambda3,Lambda4,Lambda5` (multiplicity two each) and `Lambda6` (multiplicity one). All other nonconstant modes are zero. Let `Pk` be the orthogonal real Fourier projector of harmonic k; explicitly, for k=3,4, `(Pk)_{ij}=cos(2πk(i-j)/12)/6`. Define `A_ε=A7+ε(P3-P4)` and compare **only** ε=0,+1/100 with fixed original labelled q12 LOO-HB unit-attempt dynamics, denominator N and `g_c=12/Lambda6`. No alternate microscopic law is adopted as FIN ontology: ε=1/100 is a diagnostic control.

For harmonic 3, cosine cycle is `(1,0,-1,0)`; for harmonic 4 it is `(1,-1/2,-1/2)`. Hence P3,P4 and their difference have rational entries. The standalone verifier checks using `fractions.Fraction`: `Pk^2=Pk`, `P3P4=P4P3=0`, each projector has rank two, their diagonal entries are both 1/6, row sums vanish, parity vector is annihilated, and their difference is symmetric/circulant. Therefore

- `trace A_ε=trace A7`, `diag A_ε=diag A7`, `A_ε 1=0`, `A_ε sigma=Lambda6 sigma`;
- the eigenvalue lists differ **only** by `Lambda3→Lambda3+ε`, `Lambda4→Lambda4−ε`, with Lambda5,Lambda6 fixed;
- both are PSD, rank 7, real symmetric circulant and dihedrally covariant; the largest eigenvalue is simple parity and the fast gap remains positive (indeed the smallest relative gap remains the k5 value about 0.01860477465).

These statements are algebraic except inequalities using eigenvalues, which are rechecked with `mpmath.iv` from the supplied kernel construction. No new global phase census was assumed.

## 2. Actual generator and detailed balance

For every ε, the unchanged **form** of the model is the slot-wise labelled heat bath: a slot currently in label i refreshes to label j at probability `softmax[(g_c/N) A_ε (n-e_i)]_j` per unit-rate attempt, with null refreshes allowed. For the count chain, its off-diagonal jump rate is `n_i softmax[(g_c/N) A_ε(n-e_i)]_j`. Constant diagonal gives the exact finite-count reversible weights

`π_ε,N(n) ∝ [N!/(12^N ∏ n_j!)] exp[(g_c/(2N)) n^T A_ε n]`.

Check: when n→n-e_i+e_j, the ratio of Gibbs factors matches the reversed heat-bath odds because `(A_ε)_{ii}=(A_ε)_{jj}`, while multinomial ratios and leave-one-out fields supply the remaining terms. This is the same reversible *construction*, not an assertion that full equilibrium is local.

## 3. Exact parity restriction and cancellation

For `p_m=u+(m/12) sigma`, `A_ε p_m=(Lambda6 m/12) sigma=χ m sigma`, `χ=Lambda6/12`. Consequently for both models the projected parity mean-field drift along the parity line is **exactly** `tanh(g_c χ m)-m=tanh(m)-m`, and the local free-energy potential on this line is identically

`I(m)−m²/2`,  `I(m)=((1+m)/2)log(1+m)+((1−m)/2)log(1−m)`.

It has cubic drift coefficient `−1/3`, quartic potential coefficient `+1/12` and the unchanged critical soft jump-bracket coefficient `2`, including the unit clock. Let `R` denote label rotation by two sites. On the normalized tangent, the `R`-fixed subspace is the parity line; transverse modes transform nontrivially. The scalar parity drift is rotation-invariant, so its derivative with respect to **each transverse fast vector** vanishes identically along that parity line. The hidden fast coordinates may affect finite-N brackets; they are not exactly eliminated.

## 4. Uniform local MT window, with necessary scope

The accepted MT source argument is based locally on (i) `F_ε,g(p_m)=I(m)−gχm²/2`, (ii) vanishing first transverse derivative, (iii) the fast Hessian `diag(1/p)−g A_ε` bounded below by a positive scalar for the fixed small cylinder, (iv) bounded smooth derivatives, (v) the combinatorial type-graph path and Gibbs weight inequalities. All five remain valid for the two frozen models.

For the inherited rational cylinder `|m|<.01`, `||w||_2<.0001`, at λ=0 and the original MT parameters `h=10^-10`, `b=5×10^-11`, the source lower fast Hessian bound has the same algebraic margin:

`1/((1.01)/12+0.0001) − 12(1−0.0186) > 0.08`.

Indeed the maximum transverse eigenvalue is unchanged k5, and `g_c Lambda5=12(1−gamma5)`. Along the parity line the unchanged quartic lower bound holds. Consequently the exact same local coercivity inequality `F_ε(p)≥m⁴/24+0.04||w||²` at λ=0 follows on this cylinder, as does a fixed positive free-energy buffer to its outer edges. The local gradient descent/corridor path relies on `||∂_w F||≥.08||w||`, `|∂_m F(p_m+w)−∂_m F(p_m)|≤6.14||w||²`, and the exact parity derivative; these estimates are independent of ε once the fast gap is fixed. Thus the *authorial MT path proof*, if admitted as mathematically valid in the review, transfers to the comparator with model-specific `D_{ε,N}`. The same canonical-path polynomial reflected mixing and exponential departure bound follow as existential asymptotic conclusions. No explicit `N0`, polynomial exponent K, validated stopping date or observed physical time is earned.

**Dependency warning:** the named independent MT and LI proof-completion/review files are not in the currently mounted predecessor archives. The pasted LI intake is an accepted-scope premise, but this package cannot verify the bytes of inaccessible review files. We therefore label the final local-window and consequent observation transfer as *conditional on that accepted input*, while preserving the exact algebra and local Hessian checks as independent results. No full-phase result or B4/capacity is transferred.

## 5. Shared leading parity laws and narrow no-go

With the admitted local prepared window for both models, the UF generator+fast-selection estimates transfer: both have the same limiting `dZ=−Z³/3 ds+√2 dB`, and the same local quartic stationary start `ρ(z)∝exp(−z⁴/12)`. The OBS fixed-k integrated parity record has the same `dY=Z ds+√2 dW` for k=1 with independent W. (The leading source variance and soft drift are equal; fast errors vanish by the paid stable transverse gap and symmetry cancellation.) A fixed finite quantization of two record increments at scaled times 1/4 and 1/2, with finite target such as the later sign at 1, therefore has the same limiting finite-alphabet joint law for ε=0,+.01. Every cell boundary has zero limiting probability because observation increments have continuous Gaussian conditional densities; future ZERO has probability zero in the limit, not at finite N.

For the **different** binary task of distinguishing the operator from that finite record, equal model priors give limiting ideal Bayes error 1/2 by finite-alphabet total-variation convergence. This is a narrow fixed-resolution no-go; it **does not** prove that all microscopic parity timestamps, all losses, all model-index experiments, or full q12 histories are blind. LI's one-target sign-risk theorem cannot be copied to a different model-index target without a new likelihood argument. That stronger claim is explicitly **OPEN_CROSS_MODEL**.

## New vs reused

New: exact rational matched-family checks in this execution, a comparator-specific uniform hypothesis/local-cylinder audit, narrow cross-model fixed-alphabet inference. Reused: source A7, heat bath, MT local path theorem, quartic critical diffusion, OBS integrated record, LI one-target theorem. No new experiment, natural source or inferred physical space.
