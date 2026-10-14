# SF-003 — Prepared compact-stopped stochastic-fold convergence

## Statement (local weak-convergence theorem)

Let the **exact**, supplied q12 rank-seven FIN feature matrix \(X\), its positive fold \((g_f,\theta_f,p_f)\), and the unit source-oriented null vector \(v\) satisfy the isolated-simple-fold contract of SF-001: \(p_f=\operatorname{softmax}(X\theta_f)\), \(\theta_f=g_fX^\top p_f\), and \(H_7=I/g_f-X^\top(\operatorname{diag}p_f-p_fp_f^\top)X\) has exactly one zero eigenvalue with the other six strictly positive. Let the microscopic chain have the *actual leave-one-out generator*
\[L_N f(p)=N\sum_{i,j}p_i\operatorname{softmax}_j\big(g_N A(p-e_i/N)\big)[f(p+(e_j-e_i)/N)-f(p)],\quad A=XX^\top.\]
Fix \(S,R>0\), \(\lambda\in[-\Lambda,\Lambda]\), \(g_N=g_f+\lambda N^{-2/3}\), and a deterministic \(z_0\in(-R,R)\). Initialize with balanced-rounded counts from \(p_f+z_0N^{-1/3}r\), \(r=(\operatorname{diag}p_f-p_fp_f^\top)Xv\). More generally one may use random initial types with \(Z_N(0)\Rightarrow Z_0\), support uniformly contained in \([-R+\epsilon,R-\epsilon]\) with probability tending to one, \(\mathbb E\|Q(p_N(0)-p_f)\|^2\le C/N\), and bad-initial-event probability tending to zero. Define
\[Z_N(s)=N^{1/3}v^\top(g_NX^\top p_N(N^{1/3}s)-\theta_f),\quad Q=I-r\ell^\top,\quad \ell^\top=g_fv^\top X^\top.\]
Stop this process at \(S\) and the first crossing of \(|Z_N|\ge R\) or \(\|Q(p_N-p_f)\|\ge N^{-1/4}\), freezing the value thereafter. Define the analogous frozen process \(Z\) from the solution
\[dZ_s=-g_f\left(a\lambda+\frac b2Z_s^2\right)ds+\sqrt{2g_f}\,dB_s,\quad Z(0)=z_0,\tag{A}\]
where \(a=v^\top\partial_g\nabla\Phi_g(\theta_f)|_{g_f}\), \(b=D^3\Phi_{g_f}(\theta_f)[v,v,v]\), and \(\tau=\inf\{s:|Z_s|\ge R\}\). **Then the stopped \(Z_N\) converges in law to \(Z_{s\wedge\tau\wedge S}\) in \(D([0,S],\mathbb R)\) with the Skorokhod \(J_1\) topology.** The auxiliary fast-tube exit probability is \(O(N^{-1/6})\). In particular the two exit-side indicators and exit time capped at S converge jointly at continuity points; the limiting exit at exactly S has probability zero for interior starts. No global or finite-N uniform distribution-error rate is claimed.

The theorem is conditional on the **certified exact fold/inertia premise** (which the prior intake bound to the q12 exact-decimal feature construction), not on treating a decimal midpoint as exact. It is an analytic theorem about a pre-specified mathematical Markov process, not new physical evidence.

## Proof: no independent noise or fast bath

1. **Full fast chart and preparation.** SF-001 proves \(p-p_f=xr+w\), \(x=\ell^\top(p-p_f)\), \(w=Q(p-p_f)\), \(Jr=0\), \(\ell J=0\), and that the restriction \(B\) to \(\operatorname{Ran}Q\cap\mathcal T\) is Hurwitz. It constructs a derived positive quadratic Lyapunov function and actual balanced-rounded preparations with \(Z_N(0)\to z_0\). The exact coordinate relation implies \(|x|=O(N^{-1/3})\) before soft exit.

2. **Stopped energy and rare auxiliary exit.** From the *same finite-N generator*, SF-002 proves \(L_NV\le-\kappa V+C(x^4+\Delta g^2+1/N)\), including the second-order contribution from every microscopic jump. Its stopped Dynkin estimate gives \(E\int_0^{S\wedge\tau_N/N^{1/3}}\|w_N(N^{1/3}s)\|^2ds=O(N^{-1})\), and \(P(\tau_{\rm fast}<\tau_{\rm soft}\wedge N^{1/3}S)=O(N^{-1/6})\). Hence the microscopic chain remains in a vanishing occupation neighbourhood of \(p_f\), up to the soft/horizon stopping, with probability tending to one.

3. **Actual soft semimartingale.** Write \(d_{ij}=v^\top X^\top(e_j-e_i)\) and \(q_{ij}^{(N)}(p)=p_i\operatorname{softmax}_j[g_NA(p-e_i/N)]\). Up to stopping, \(Z_N\) is a semimartingale with exact compensated jump martingale \(M_N\):
\[Z_N(s)=Z_N(0)+\int_0^{s\wedge\hat\tau_N}\!\beta_N(p_N(N^{1/3}u))du+M_N(s),\]
where \(\hat\tau_N=\tau_N/N^{1/3}\) and
\[\beta_N(p)=N^{2/3}g_N\sum_{ij}q_{ij}^{(N)}(p)d_{ij},\quad
\frac{d\langle M_N\rangle_s}{ds}=g_N^2\sum_{ij}q_{ij}^{(N)}(p)d_{ij}^2\quad(s<\hat\tau_N).\tag{B}\]
Every rescaled jump has absolute size at most \(C N^{-2/3}\). The sum of its absolute third powers per unit scaled time is \(O(N^{-2/3})\), since the total microscopic refresh-attempt rate is N. This is a bound on **actual jump characteristics**, not a diffusion inserted by hand.

4. **Drift expansion with paid fast remainder.** The source soft coordinate is \(\theta_N=g_NX^\top p\). The leading exactly closed retained drift is \(-g_N\nabla\Phi_{g_N}(\theta_N)\), with finite-N LOO error \(O(N^{-1})\). Writing \(\theta_N=\theta_f+N^{-1/3}Z_Nv+\eta_N\), the chart identities give \(v^\top\eta_N=0\) and \(\|\eta_N\|\le C(\|w_N\|+|\Delta g|)\). The fold conditions \(\nabla\Phi_{g_f}(\theta_f)=0\), \(v^\top D^2\Phi_{g_f}(\theta_f)=0\) yield the second-order Taylor formula
\[v^\top\nabla\Phi_{g_N}(\theta_N)=a\Delta g+\tfrac12 bN^{-2/3}Z_N^2+\rho_N,\]
\[|\rho_N|\le C\big(|x|\|w\|+\|w\|^2+|x|^3+|\Delta g||x|+|\Delta g|\|w\|+\Delta g^2\big),\]
with constants uniform on the compact stopped domain. SF-002 provides each of the seven expected integrated bounds after multiplication by \(N^{2/3}\); the largest is \(O(N^{-1/6})\). Therefore
\[\mathbb E\int_0^{S\wedge\hat\tau_N}\!|\beta_N(p_N(N^{1/3}u))+g_f(a\lambda+\tfrac b2Z_N(u)^2)|\,du\longrightarrow0.\tag{C}\]
The gain prefactor correction \(g_N-g_f=O(N^{-2/3})\) also vanishes on the bounded domain. No unbounded path-dependent little-o term remains.

5. **Bracket and third characteristic.** On the stopped domain \(p_N\to p_f\) and \(g_N\to g_f\) uniformly, by the deterministic tube and soft coordinate bound (with jump overshoots). The LOO categorical probabilities converge uniformly to \(p_f\). The bracket integrand (B) therefore converges uniformly to
\[g_f^2\sum_{ij}p_{f,i}p_{f,j}d_{ij}^2=2g_f^2v^\top X^\top\Sigma_fXv=2g_f.\tag{D}\]
Consequently the bracket of \(M_N\) converges uniformly in probability to \(2g_f(s\wedge\hat\tau_N)\) in the stopped representation. There is no covert use of equilibrium Gaussian auxiliary trajectories or an isotropic covariance away from the soft fold channel.

6. **Tightness and identification.** In (C), the drift is a bounded continuous function of \(Z_N\in[-R-o(1),R+o(1)]\) plus a total-variation error tending to zero in L1. In (B) the predictable brackets are uniformly bounded on [0,S], and jumps vanish. The semimartingale tightness criterion (or Aldous' criterion applied to the bounded drift and square-integrable compensated martingale) gives C-tightness in \(D([0,S])\). For \(\varphi\in C_c^3((-R,R))\), compensated Dynkin's formula and a third-order Taylor remainder bounded by the vanishing third characteristic imply that every subsequential limit solves the local martingale problem with
\[G\varphi(z)=-g_f(a\lambda+bz^2/2)\varphi'(z)+g_f\varphi''(z).\tag{E}\]
The drift and diffusion coefficients are smooth and bounded on the closed interval and \(\sqrt{2g_f}>0\); locally extend them Lipschitz outside [-R,R] and stop at its first boundary. Pathwise uniqueness of the extended SDE plus uniqueness of its stopped law identifies every subsequence. Thus the full sequence converges **before removing auxiliary stopping**.

7. **Removal of auxiliary stopping and exit map.** The probability of an auxiliary fast exit before soft exit/S is \(O(N^{-1/6})\), so dropping this extra stop changes the law of any bounded stopped-path functional by at most that probability (couple on the original path). To justify hitting-time convergence without circularly applying a hitting map to an already frozen path, repeat the SF-001/002 estimates with any one fixed larger soft radius \(R^{\prime}>R\), using the same exact initial class (constants may depend on \(R^{\prime}\)). The martingale problem then proves convergence on the interval before the *larger* exit. On those paths the first exit from (-R,R) is an interior first-crossing functional, almost surely continuous at the limiting uniformly elliptic diffusion: Brownian regularity ensures immediate excursions across either level after the first hit, and \(P(\tau_R=S)=0\) for every fixed positive S and interior start. Apply continuous mapping for the inner \(R\) exit to paths localized at \(R^{\prime}\), then freeze at \(\pm R\). This also delivers joint convergence of \(\tau_R\wedge S\) and side; the microscopic overshoot of order \(N^{-2/3}\) vanishes. The extra larger-radius construction is a proof localization and does not alter the declared observable law at R.

Hence the compact-stopped critical process theorem follows. It is **not** a theorem about unbounded-domain explosion, global phase switching, equilibrium starting at the fold without conditioning, finite-size mixing, or real physical time.

## Audit / scope

The key asymptotic proof constants are existential, derived from smoothness and positivity of the inherited fold chart, and are not claimed to be outward numerical bounds. A rigorous independent **check of the inherited numerical Krawczyk certificate and its exact matrix-byte provenance** is not repeated by this campaign; our code checks only a binary64 diagnostic root. Under a false or mismatched inherited simple-fold premise, this result would be a conditional theorem rather than a specific FIN theorem. The application does not use the fixed-time path LDP as a substitute for growing-horizon tightness.
