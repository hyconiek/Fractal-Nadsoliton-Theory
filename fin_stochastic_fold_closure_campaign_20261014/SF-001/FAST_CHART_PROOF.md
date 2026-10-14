# SF-001 — Exact occupation chart, spectral restriction and nonempty preparations

**Status: analytic theorem, conditional only on the inherited isolated simple-fold certificate for the precise q12 feature matrix.** No numeric binary64 fold root is substituted into the asymptotic definitions.

Let \(X\in\mathbb R^{12\times7}\), \(A=XX^\top\), \(p_f=s_{g_f}(p_f)>0\), \(\theta_f=g_fX^\top p_f\), \(\Sigma_f=\operatorname{diag}p_f-p_fp_f^\top\), \(M_f=X^\top\Sigma_fX\), and \(H_7=I_7/g_f-M_f\). The existing local fold certificate states that \(H_7\) has a *simple* zero eigenvector \(v\), \(\|v\|_2=1\), and six positive eigenvalues. Its normalization/sign is retained. In this proof \(p_f,g_f,v\) denote the implicitly isolated **exact** objects; the numerical chart test is only a diagnostic.

On the tangent \(\mathcal T=\{h:1^\top h=0\}\) the empirical-refresh limiting drift \(b(p,g)=s_g(p)-p\) has derivative
\[J=g_f\Sigma_fXX^\top-I_{\mathcal T}.\]
Define \(r=\Sigma_fXv\), \(\ell^\top=g_fv^\top X^\top\). Since \(M_fv=v/g_f\), \(Jr=0\), \(\ell^\top J=0\) *on \(\mathcal T\)* and \(\ell^\top r=g_fv^\top M_fv=1\). Consequently \(Q=I-r\ell^\top\), restricted to \(\mathcal T\), is an idempotent oblique projection, with \(Qr=0\), \(\ell^\top Q=0\). Put \(x=\ell^\top(p-p_f)\), \(w=Q(p-p_f)\); then \(p-p_f=xr+w\) uniquely and \(w\in\mathcal W=\{h:1^\top h=0,\ell^\top h=0\}\) (dimension 10).

The map \(h\mapsto(X^\top h,Rh)\), where \(R\) projects orthogonally off \(\operatorname{span}(1,X)\), is injective on \(\mathcal T\) and dimension-counting makes it an isomorphism to seven retained plus four hidden coordinates. The *leading* nonlinear deterministic system is
\[\dot\theta=-g\nabla\Phi_g(\theta),\quad \dot h=R\,\operatorname{softmax}(X\theta)-h,\qquad \theta=gX^\top p,\ h=Rp.\]
The corresponding full-tangent Jacobian is block triangular, with blocks \(-g_fH_7\) and \(-I_4\). Therefore its spectrum is exactly one zero, six strictly negative numbers \(-g_f\lambda_k(H_7)\), and four instances of \(-1\). This is a *deterministic* tangent result, not a stochastic or physical spectral gap. As \(r\) spans the only null direction and \(\mathcal W\) is invariant under \(J\), the restriction \(B=J|_{\mathcal W}\) is Hurwitz.

Fix any Euclidean-orthonormal basis \(E\in\mathbb R^{12\times10}\) of \(\mathcal W\) and set \(B_E=E^\top JE\). A **derived**, not fitted, quadratic Lyapunov matrix is
\[P=\int_0^\infty e^{B_E^\top t}e^{B_E t}\,dt.\]
This integral converges, is strictly positive definite, and uniquely solves \(B_E^\top P+PB_E=-I_{10}\). Define \(V(w)=(E^\top w)^\top P(E^\top w)\). Exact constants \(m=\lambda_{\min}(P)>0\), \(M=\lambda_{\max}(P)<\infty\) satisfy \(m\|w\|^2\le V\le M\|w\|^2\). Numeric \(P\) and eigenvalues in NUMERICAL_CHART.json are not interval-validated bounds.

### Prepared types and exact gain conversion
Fix finite \(R,S>0\), \(\lambda\) in a compact \([-\Lambda,\Lambda]\), and \(z_0\) in a compact subinterval \([-R+\epsilon,R-\epsilon]\), \(\epsilon>0\). Let \(g_N=g_f+\lambda N^{-2/3}\). Since \(1^\top r=0\) and \(p_f>0\), for large N the vector \(q_N=p_f+N^{-1/3}z_0r\) is strictly positive and sums to one. Balanced (largest-remainder) rounding to nonnegative integer counts \(n(N)\) with \(\sum_i n_i=N\) produces \(p_N(0)=n(N)/N=q_N+\varepsilon_N\), with \(\|\varepsilon_N\|_\infty\le 1/N\) (and \(1^\top\varepsilon_N=0\)). Thus \(w_N(0)=Q\varepsilon_N=O(N^{-1})\), \(E\|w_N(0)\|^2=O(N^{-2})\), and interior/type support is exact. A mixture of bounded such starts extends to a convergent random Z0 under uniform moment bounds; the theorem may simply use deterministic z0.

For arbitrary tangent \(p-p_f=xr+w\), \(\theta_N=g_NX^\top p\), there is the **exact** identity
\[v^\top(\theta_N-\theta_f)=\frac{g_N}{g_f}x+(g_N-g_f)v^\top X^\top p_f.\]
Because \(g_N-g_f=O(N^{-2/3})\), the rounded preparation has \(Z_N(0):=N^{1/3}v^\top(\theta_N(0)-\theta_f)\to z_0\). On \(|Z_N|\le R\), it follows that \(|x|\le C_{R,\Lambda}N^{-1/3}\) for all sufficiently large N, *without any stochastic estimate*.

The two stopping conditions are the first crossing of \(|Z_N|\ge R\) and \(\|w_N\|\ge \delta_N=N^{-1/4}\); also stop at model time \(N^{1/3}S\). Rescaled soft jumps are at most \(C N^{-2/3}\), while fast jumps are at most \(C/N\), so boundary overshoots vanish in their respective scales. The exact generator remains the leave-one-out one; its normalization is N throughout.

**Scope:** exact chart and deterministic stability proof given the inherited certified fold. This does not prove a globally stable phase, a uniform microscopic mixing gap, a unique physical clock or an independent source for A7.
