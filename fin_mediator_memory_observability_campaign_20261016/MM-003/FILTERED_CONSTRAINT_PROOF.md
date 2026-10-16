# MM-003 — Reference-pole-filtered constraint (exact Gaussian algebra)

Let c=Lambda3/Lambda6, gamma=1-c, 0<c<1. The diagnostic global fast drift matrix is J=[[-1,1],[kappa*c,-kappa]], the diffusion bracket diag(2,2*kappa*c). The Lyapunov equation J Sigma+Sigma J^T+Q=0 has stationary Sigma=gamma^-1[[1,c],[c,c]], as independently checked by exact symbolic simplification. The eigenvalues of -J are

    beta_s,f = [(1+kappa) ∓ sqrt((1-kappa)^2+4*kappa*c)]/2.

Evaluating the characteristic polynomial at gamma and 1 shows beta_s<gamma<1<beta_f for positive c and kappa. Stationary F covariance is gamma^-1[w_s exp(-beta_s|t|)+w_f exp(-beta_f|t|)], with w_s=(beta_f-gamma)/(beta_f-beta_s), w_f=1-w_s, both >0. The independently sampled private U contributes c_alpha exp(-|t|), c_alpha=(1-alpha)/alpha. Thus the partial diagnostic readout has THREE component rates (two feedback poles and the private pole), whereas the native readout has the native global rate gamma plus the same private rate1. Full alpha=1 omits the private term.

A common positive model-time rescaling t->nu t changes x=e^-beta_s nu h, y=e^-beta_f nu h, q=e^-nu h without changing the ordering x>q>y. For active unnormalized C_j=A x^j+B y^j+D q^j, define d_j=C_{j+1}-q C_j for j=1,2,3. Algebra gives

    H_ref=d1*d3-d2^2=A*B*x*y*(x-q)*(y-q)*(x-y)^2 < 0.

The coefficient D cancels exactly, regardless of observation fraction. If C_j=A x^j+D q^j (native partial) or C_j=A x^j (native full), H_ref is identically zero. In contrast ordinary C1*C3-C2^2 for native partial equals A*D*x*q*(x-q)^2>0 and is a false-positive memory indicator.

**Bounded-sign law:** At each mode r and lag j=1..4, let s_rj=E[sign V_r(0) sign V_r(j)], with sign(0)=0. For a nonsingular Gaussian pair s=(2/pi)arcsin(rho), so R_rj=sin(pi s_rj/2)=Corr(V_r(0),V_r(j)). Set q=R_null,2/R_null,1 only when R_null,1≠0 and the null r=1 one-rate checks are supported, including R_null,3=R_null,1*q^2. Replace C_j by the normalized R_active,j and evaluate H_ref with this q. Because normalization divides all C_j by the same C_0>0, the sign of H_ref is unchanged and the expression is a functional of bounded ternary observations *only*. It is not an estimate from one five-time observed tuple.

**Gauge restrictions:** Each observed channel may have its own unknown FIXED nonzero real multiplier: sign products cancel its sign. A single common change of model-time units changes q and active rates consistently and is removed by null calibration. No offset, time-dependent sign-changing gain, detector noise, label reassignment or independent channel clock is included. Under an unsupported null reference, the rule returns REFERENCE_UNSUPPORTED, never refits q.

**Margin and scope:** At alpha=1/2,h=1,kappa1/4, nominal H_ref is strictly negative. The adjacent interval file uses mpmath.iv with an expressly ASSUMED input gamma interval [0.16257283698729,0.16257283698731] inferred from the prior nominal spectrum; source spectral enclosure must be independently bound for a certified physical-parameter result. No uniform margin for alpha->0 or q->0,1 has been proved.
