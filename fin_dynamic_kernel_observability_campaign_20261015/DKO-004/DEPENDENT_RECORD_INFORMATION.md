# DKO-004 — Gaussian observation-fraction cost with dependent records

**Standalone scope:** This theorem concerns the *stationary limiting Gaussian process family* from DKO-002. FIRST take N→∞ for a chosen fixed α∈(0,1) and a fixed finite fast-time window; THEN analyze sequences α→0 and M→∞ among these Gaussian experiments. This is NOT a joint triangular array of FIN parents with α_N→0, T_N→∞.

For the two fixed mode3 rates a=γ3(A7), b=γ3(control), both in (0,1), and fixed sampling interval h>0, define W=√α V (same signs). From the two-OU law its covariance matrix for M stationary records at times 0,h,...,(M−1)h is

K_{γ,α}^{(M)}(i,j)= (α/γ)e^{−γh|i−j|} + (1−α)e^{−h|i−j|}.

This retains the exact temporal dependence between observations. For zero-mean Gaussian laws P_a^M,P_b^M, the KL divergence is EXACTLY

D(P_a^M || P_b^M)=(1/2)[tr((K_b)^{-1}K_a−I)−log det((K_b)^{-1}K_a)].

**Two-sided analytic cost.** Let q=e^{-h}, q_min=e^{-min(a,b)h}; write R_x(i,j)=x^{|i−j|} and D= a^{-1}R_{e^{-ah}}−b^{-1}R_{e^{-bh}}. Then K_a−K_b=αD. Define d_l=a^{-1}e^{-ah l}−b^{-1}e^{-bh l}. The norms satisfy

M d_0² ≤ ||D||_F² ≤ M[d_0²+2Σ_{l≥1}d_l²]=M D_2,

and ||D||_op≤ D_1=|d_0|+2Σ_{l≥1}|d_l|<∞. The strict inequality d_0=1/a−1/b≠0 is crucial.

For α≤1/2, positive private OU component gives a uniform lower eigenvalue ℓ=(1/2)(1−q)/(1+q)>0 for both covariance matrices, since the Toeplitz correlation R_q has eigenvalues in [(1−q)/(1+q),(1+q)/(1−q)] independently of M. An upper eigenvalue bound is L=(1+q)/(1−q)+(1/min(a,b))(1+q_min)/(1−q_min). Restrict α≤α0=min(1/2,ℓ/(2D_1)); then E=(K_b)^{-1/2}(K_a−K_b)(K_b)^{-1/2} has ||E||op≤1/2. From φ(x)=x−log(1+x), φ''(x)=(1+x)^−2, for |x|≤1/2,

(2/9)x² ≤ φ(x) ≤ 2x².

Therefore (1/9)||E||_F² ≤ KL ≤ ||E||_F² and by eigenvalue comparison

[ d_0²/(9L²) ] Mα² ≤ KL ≤ [ D_2/ℓ² ] Mα².

The constants depend only on fixed a,b,h (not M or α in (0,α0]). This proves **KL=Θ(Mα²)** with dependent Toeplitz records, not an iid effective-sample-size heuristic. The reverse KL has analogous bounds.

Pinsker's inequality implies TV(P_a^M,P_b^M)≤√(KL/2) and ideal equal-prior model-ID risk→1/2 when Mα²→0, even if the observer sees all Gaussian amplitudes; by data processing, the same impossibility holds for signs alone. The next document separately proves achievability for Mα²→∞: large KL alone would not suffice.
