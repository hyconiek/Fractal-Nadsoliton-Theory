# DH-001 — Exact dense-complement Hall theorem on actual count fibres

Date 2026-10-17. Classification: NEW RESTRICTED THEOREM, NOT ALL DENSE HALL. The inherited MT/SC/PC assumptions listed here must be independently bound against the actual 2026-10-17 reviews before admission. No new microscopic rule, domain, threshold, or graph is introduced.

## Model and quantifiers

For N in 12ℕ and every actual neighbouring critical parity fibre X_l -> X_{l+1} with both z_l,z_{l+1} in [-1,1], use the *original* MT D_N and departure E_N; parity-conserving 12-element group G_+ acts on counts. H contains every legal within-phase edge and interphase edges incident to E_N. Use the exact PC04 nonnegative corrected marginals O'_l,I'_{l+1} and their fibre-total identity J_l = Σ_{X_l}O'_l = Σ_{X_{l+1}}I'_{l+1}. The subgroup preserves the entire contract, so orbit unions are typed.

The inherited compact input (not proved from the archived inputs in this execution) is: at sufficiently large admissible N, **uniformly** in the band, (i) for x not in E_l, O'_l(x)≥(173/1000) N W_x; (ii) for y in E_{l+1}, I'_{l+1}(y)≤N W_y; (iii) every y outside E_{l+1} has at least nine legal *downward* predecessors in H; (iv) the conditional right-fibre boundary estimate W(E_{l+1})/W(X_{l+1})≤poly(N)e^{-bN} for some b>0; and (v) the compact fibre mass ratios W(X_l)/W(X_{l+1}) stay between two positive constants. The proof below isolates these inputs rather than attributing them to the algebraic test.

## Theorem 1 — Finite exact separated-hole certificate

Let A⊂X_l be a G_+-invariant collection of distinct count states satisfying:

1. every x∈A lies outside the true departure set E_l;
2. ||x−x'||_1>4 for every pair of distinct x,x'∈A;
3. (173/1000) W(A)≥W(E_{l+1}), with W the *actual* count-Gibbs multiplicity weights.

Put T=X_l\A. Then for sufficiently large N as above,

 O'_l(T) ≤ I'_{l+1}(N_H(T)).

**Proof:** Let R=X_{l+1}\N_H(T). Because two distinct downward predecessors of one right count y have the form y+e_a−e_b and y+e_c−e_d, their L1 distance is at most 4. Thus y has at most ONE predecessor in A. If y is not in E_{l+1}, its ≥9 legal predecessors cannot all belong to A, and hence one belongs to T: y∈N_H(T). Consequently R⊂E_{l+1}; this uses the actual full-H neighbourhood, including E-cross-phase edges, with no silently deleted edges. At E, the frozen correction is I'_y=I_y+[-s_y−I_y]_+=max(I_y,−s_y), where 0≤I_y≤.35 N W_y and |s_y|≤N W_y by native count proposal rate ≤N and parity jump ≤2/N. Therefore I'(R)≤N W(E_{l+1}). On A (which avoids E), I'_x≥.18 NW_x and |s_x|≤.007 NW_x, so O'_l(A)=I'_l(A)+s(A)≥.173 N W(A). The total-current identity yields

 I'(N_H(T))−O'(T) = O'(A)−I'(R)≥N[(173/1000)W(A)−W(E_{l+1})]≥0.

All statements refer to exact finite fibres; no approximate statewise drift was substituted. □

## Theorem 2 — An explicit, unbounded-N family of *dense G_+-orbit union cuts* satisfying Hall

Fix any N∈12ℕ and any critical fibre with |z_l|≤1 and z_{l+1}≤1. The even and odd count totals K_e,K_o are integers. Distribute each K_e,K_o as evenly as possible among the six even and six odd labels, obtaining n^0_l with |n^0_{j}−K_{parity(j)}/6|≤1. Define for r=0,...,5 the integer vector t_{2r}=t_{2r+1}=2r−5 and k_N=⌊N^{1/3}⌋. Set n^*_l=n^0_l+k_N t, and A_{N,l}=G_+ n^*_l. Then for sufficiently large N:

- all components of n^*_l are positive, Σ n^*_l=N, and its parity is m_l;
- A_{N,l} lies in the *interior* of actual D_N, far from E_N, and is wholly in phase A with the original inclusive cutoff R≤T(z);
- the orbit consists of 12 distinct members, each pair at L1-distance at least 36 k_N−12>4 (G_+ acts by even rotations and parity-preserving reflections; the fixed pattern t has minimum orbit distance 36);
- by the admitted uniform compact lattice/Stirling law in MT and the central Gaussian transverse profile, W(A_{N,l})≥c N^{-5} W(X_l) for some c>0, uniformly in the band. This bound is inherited from the local Gibbs theorem; it is not established by the finite verifier;
- the MT departure barrier implies W(E_{l+1})≤C N^a e^{-bN} W(X_l) uniformly on this compact band, for some a,b,C>0.

Thus (173/1000)W(A_{N,l})≥W(E_{l+1}) for all sufficiently large admissible N, and Theorem 1 proves Hall for the genuinely dense G_+-orbit-union T_{N,l}=X_l\A_{N,l}. This family was NOT covered by the source theorem requiring the nine selected successors of *each member of T* to overlap at most four times. Nothing here claims that arbitrary dense orbit-unions satisfy Hall.

**Verification of geometry:** Σ_even t=Σ_odd t=0. For the integer k3 modes C(n^*_l)=−6k_N+O(1), S(n^*_l)=−6k_N+O(1); hence R_N=2(C²+S²)/N=O(N^{-1/3}). The frozen threshold T(z) has positive infimum on [-1,1] because gamma3>0 and p_A(z)<1. The soft coordinate is fixed in the allowed band and ||n^*/N−p_{m_l}||=O(N^{-2/3}). The admitted MT free-energy Taylor bound is F=O(N^{-1}) near this representative, below the fixed departure barrier h/2, so it is interior. The six within-parity t-coordinates form a ramp [-5,-3,-1,1,3,5], whose 12 G_+ images have pairwise L1-distance ≥36 before multiplying by k_N. Rounded base differences cost at most 12 in L1. Normalized central Gibbs fibre atoms in dimension ten have mass ≥cN^{-5}; a transverse shift O(N^{1/3}) changes the quadratic exponent by O(N^{-1/3}).

**Status:** An honest new *subclass of dense orbit cuts*, not a proof of complete DH-001. The failure mode for general dense orbit cuts is that their removed complement may have large internally clustered connected regions, so right-hole predecessors can all lie in the complement even away from E. The separated-hole argument then fails, and the weighted Hall condition remains open.
