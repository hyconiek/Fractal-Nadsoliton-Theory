# SF-002 — Generator-level Lyapunov inequality and negligible fast-tube exit

**Status: analytic local asymptotic proof under SF-001's certified simple-fold hypothesis and the exact finite-copy heat-bath generator.** All constants below are existential, independent of \(N\), and depend only on the frozen \(X,A,p_f,g_f,E,P,R,S,\Lambda\) and a fixed sufficiently small interior neighbourhood. No numerical finite-N constant or rate is claimed.

For counts \(p\), write \(h_{ij}=(e_j-e_i)/N\), and the actual generator
\[L_N f(p)=N\sum_{i,j}p_i s_j(g_NA(p-e_i/N))[f(p+h_{ij})-f(p)].\]
The exact mean drift is
\[b_N(p,g)=\sum_{ij}p_i s_j(gA(p-e_i/N))(e_j-e_i).\]
Write \(b(p,g)=s(gAp)-p\). Softmax has bounded derivatives of every fixed order on the compact simplex \(p\in\Delta_{11}\) and \(g\) in a fixed compact interval. Its derivative satisfies \(\|D\operatorname{softmax}\|_{1\to1}\le 1/2\). Therefore, for a fixed finite \(A\),
\[\sup_{p\in\Delta}\|b_N(p,g)-b(p,g)\|\le C/N,\]
including boundary states (terms weighted by \(p_i=0\) are zero). This is an honest LOO correction and is not obtained by altering the microscopic generator.

Set \(p=p_f+x r+w\), \(w\in\mathcal W\). The smooth expansion \(b(p,g)=J(xr+w)+O(\|xr+w\|^2+|g-g_f|)\) combined with \(Jr=0\), \(QJw=Jw\) yields, uniformly on the SF-001 stopped domain,
\[Qb_N(p,g_N)=Bw+u_N,\quad\|u_N\|\le C[x^2+|g_N-g_f|+|x|\|w\|+\|w\|^2+N^{-1}].\tag{1}\]
This uniform bound follows from the second-derivative integral remainder on a fixed compact interior neighbourhood and is independent of \(N\). On the stopped domain, \(|x|\le C_RN^{-1/3}\) and \(\|w\|\le N^{-1/4}\), both tending to zero, so the cross terms are absorbable.

Let \(y=E^\top w\), \(V=y^\top Py\), \(P>0\) from SF-001. For the *linear* evolution \(Bw\), \(\nabla V\cdot Bw=-\|w\|^2\). Under a transition \(i\to j\) the soft/fast chart changes by \(\Delta w=Q(e_j-e_i)/N\), hence \(\|\Delta w\|\le C/N\). Since V is quadratic, its generator is **exactly**
\[L_NV=2y^\top P E^\top Qb_N+N\sum_{i,j}p_i s_j(g_NA(p-e_i/N))\,\big\|P^{1/2}E^\top Q(e_j-e_i)/N\big\|^2.\tag{2}\]
The last term is bounded by \(C/N\) because \(\sum_{ij}p_i s_j=1\). For the first term, Young's inequality and (1), choosing the neighbourhood small enough to absorb the \((|x|+\|w\|)\|w\|^2\) terms, give fixed \(\kappa,C>0\) such that
\[L_NV\le-\kappa V+C[x^4+|g_N-g_f|^2+N^{-1}].\tag{3}\]
For example \(-\|w\|^2+C\|w\|(x^2+\Delta g)+C(|x|+\|w\|)\|w\|^2\le-\frac12\|w\|^2+C(x^4+\Delta g^2)\), and \(\kappa=1/(2M)\), after shrinking the neighbourhood if necessary. No specific numerical \(\kappa\) is asserted.

### Correct stopped Dynkin estimate
Let \(\tau_N\) be the FIRST soft or fast exit, bounded above by \(T_N=N^{1/3}S\) (physical model time), and \(\delta_N=N^{-1/4}\). Define \(\bar\tau_N=\tau_N\wedge T_N\). For all large N the entire stopped domain lies in the fixed neighbourhood where (3) holds. Dynkin's formula for the bounded finite-state chain and nonnegative V yields
\[\mathbb E V(p_{t\wedge\bar\tau_N})+\kappa\mathbb E\int_0^{t\wedge\bar\tau_N}V(p_u)du\le\mathbb EV(p_0)+C\big(N^{-4/3}+N^{-1}\big)T_N.\tag{4}\]
Under \(\mathbb E\|w_0\|^2\le C_0/N\), its RHS is \(O(N^{-2/3})\), with the exponent dominated by the actual jump bracket \(N^{-1}\times N^{1/3}\). In particular,
\[\mathbb E\int_0^{\bar\tau_N}\|w_u\|^2du=O(N^{-2/3}),\qquad
\mathbb E\int_0^{\bar\tau_N/N^{1/3}}\|w_{N^{1/3}s}\|^2ds=O(N^{-1}).\tag{5}\]
**Crucially (4) retains the nonnegative EXIT VALUE V at the stopping time; it does not claim an exponentially decaying expectation of a frozen nonzero absorbing state.**

At the event \(F_N\) that the fast tube exits before the soft boundary or horizon, the actual stopped \(\|w_{\bar\tau_N}\|\ge\delta_N\), including a possible overshoot. Hence \(V\ge m\delta_N^2=mN^{-1/2}\) on \(F_N\). Equation (4) immediately gives
\[\mathbb P(F_N)\le \frac{C N^{-2/3}}{m N^{-1/2}}=O(N^{-1/6}).\tag{6}\]
Possible initial bad events are absent for rounded deterministic preparations and otherwise must be assigned probability o(1) separately. The rare multiple-threshold same jump may be resolved by a fixed tie rule (fast first); jumps vanish and this does not alter the limiting soft law. There is NO inference from integrated occupation alone to exit probability; the exit-value term is essential.

### Scaled drift error ledger
Write \(x_s=x(p_{N^{1/3}s})\) and \(w_s=w(p_{N^{1/3}s})\) up to stopping. The accelerated soft generator remainder, after multiplying the physical drift by \(N^{2/3}\), is bounded by a constant times
\[N^{2/3}\big[|x_s|\|w_s\|+\|w_s\|^2+|x_s|^3+\Delta g\,|x_s|+\Delta g\,\|w_s\|+(\Delta g)^2+N^{-1}\big].\]
By (5), Cauchy–Schwarz, \(|x_s|\le C N^{-1/3}\), \(\Delta g=O(N^{-2/3})\) and \(s\le S\), the EXPECTED time-integrated terms are respectively
\[O(N^{-1/6}),\ O(N^{-1/3}),\ O(N^{-1/3}),\ O(N^{-1/3}),\ O(N^{-1/2}),\ O(N^{-2/3}),\ O(N^{-1/3}).\tag{7}\]
Every exponent is negative. There is no unbounded, unexplained little-o remainder.

Moreover the stopped occupation satisfies pathwise
\[\sup_{u\le\bar\tau_N}\|p_u-p_f\|\le C_RN^{-1/3}+N^{-1/4}+O(N^{-1}),\tag{8}\]
where an O(N^-1) boundary overshoot is included. Thus conditional predictable brackets may be evaluated at \(p_f\); no stationarity assertion or Gaussian auxiliary covariance is used.

**Exact limitation:** bounds are asymptotic and existential. They do not provide a certified numerical \(N_0\), an all-operation interval \(P\) enclosure or an empirical rate. The inherited certified simple-fold root/inertia is a separate premise; these estimates do not source the kernel.
