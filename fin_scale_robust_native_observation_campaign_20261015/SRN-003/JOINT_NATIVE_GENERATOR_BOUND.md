# SRN-003 — Four-coordinate finite jump generator, including cross-mode brackets

Let q=12, r in {1,3}, phi=(F1,R1,F3,R3) where F_r=N^-1/2 sum_all f_r(x_a), R_r=k^-1/2 sum_a(1_{a in B}-alpha) f_r(x_a). For an old label i and new label j at slot a, the **actual LOO HB** jump vector is
 Delta_a(i,j)=(f1(j)-f1(i),f3(j)-f3(i)) tensor (N^-1/2,(1_B(a)-alpha)k^-1/2),
with rate Q_j(n-e_i)=softmax_j[(g_c/N)A(n-e_i)]. Null refresh increments are zero. Reflection suppresses only outward total-count jumps.

For r,s in {1,3}, set J_rs(i,n)=sum_j Q_j(n-e_i)[f_r(j)-f_r(i)][f_s(j)-f_s(i)]. Native bracket blocks have exact identities
 Gamma(F_r,F_s)=N^-1 sum_a J_rs(x_a,n),
 Gamma(R_r,R_s)=k^-1 sum_a(1_B-alpha)^2 J_rs(x_a,n),
 Gamma(F_r,R_s)=(Nk)^-1/2 sum_a(1_B-alpha)J_rs(x_a,n).
These formulas include all r!=s blocks. At uniform counts J_rs(i)=delta_rs+f_r(i)f_s(i), not a pointwise constant unless old-label feature squares happen to be constant. Orthogonality E_u(f_r f_s)=delta_rs gives limiting 2delta_rs and 2(1-alpha)delta_rs, zero mixed global/private and cross-mode terms. Under conditional hypergeometric block sampling, mixed deviations are O_L1(N^-1/4+k^-1/2) after the inherited E|p-u|<=C N^-1/4 and the block conditional variance bounds. The OU generator is therefore
 G= sum_{r in {1,3}} [-gamma_r F_r partial_Fr + partial_FrFr - R_r partial_Rr +(1-alpha)partial_RrRr],
 gamma1=1, gamma3=1-Lambda3/Lambda6. OU drivers are independent *in the reference law only*; no finite-N independence is asserted.

The exact drift for R_r satisfies L R_r+R_r=k^-1/2 sum_a(1_B-alpha)[a_r(n-e_{x_a})-a_r(n)], where a_r(n)=sum_j f_r(j)Q_j(n). Since sum_a |1_B-alpha|=2k(1-alpha) and a_r is C/N-Lipschitz, |L R_r+R_r|<=C sqrt(k)/N. For F_r, parity-line symmetry cancels pure soft forcing, leaving stationarity-mean error O(N^-1/4) by E m²=O(N^-1/2), E||w||²=O(N^-1), Cauchy-Schwarz, with finite LOO error O(N^-1/2). The 4D absolute third-jump term is O(k^-1/2+N^-1/2), since block events contribute k*k^-3/2 and outside events contribute N*(sqrt(k)/N)^3.

For any v with globally bounded partial derivatives of orders1,2,3, Taylor and the true MT departure layer B_D give
 E_piD |L_D(v o phi)-(Gv)o phi| <= C [N^-1/4+k^-1/2+(N/sqrt(k))pi_D(B_D)] sum_{1<=|beta|<=3} ||D^beta v||_infty.
The implied constants are for the two frozen operators, fixed positive local-cylinder c, and 1<=k<=N/2. This is an L1 **stationary** generator estimate, not a uniform statewise or all-time Gaussian process theorem. It is conditional on the inherited MT uniform local moments, positive domain and departure-layer bound for both models. The rho_eta smoothing used below is *not* a microscopic perturbation.
