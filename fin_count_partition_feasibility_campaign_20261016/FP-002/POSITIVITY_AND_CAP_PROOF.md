# FP-002 — a new necessary continuum rate-cap bound; **not** a construction

Assume for contradiction-testing only that a family of **exactly** feasible, closed-phase reversible diagnostic generators exists on the true D_N, with full native parity drift, uniform exit-rate cap `CN` and stationary Gibbs initialization. For fixed `S`, use only the original KR stationary full-state drift estimate and the original MT joint static limit. Write `Z=N^(1/4)m`, `rho(z)=K exp(-z^4/12)`, `b(z)=-z³/3`, `q_A(z)=p(z)rho(z)` and `q_B(z)=(1-p(z))rho(z)`, where `p=1/2+exp(-z²)/8`.

The time-rescaled generator is `sqrt(N) L_tilde`. At each closed phase P, apply exact phase-conditioned stationarity to bounded smooth compactly supported `phi(Z)` and the Taylor expansion with jump size <=`2N^(-3/4)`. The per-time absolute third jump characteristic is `O(N^(-3/4))` from the `CN` bound. KR's L1 drift defect vanishes. With `G_N=N Gamma_tilde(m,m)`, 

`0= E_pi[1_P { b(Z) phi'(Z)+ (1/2)G_N phi''(Z)}]+o(1)`.

Crucially `0<=G_N<=4C`, so along subsequences the finite measures `E[1_P G_N 1_{Z in dz}]` have limits `nu_P` dominated by `4C q_P(z) dz`. The MT joint limit and phase-boundary Gaussian zero-atom statement give the marginal phase law `q_P dz`. Distributionally,

`nu_P''=2(b q_P)'`.

As `nu_P` is a finite positive measure dominated by an integrable positive density, the only solution is the absolutely continuous measure `nu_P(dz)=2 J_P(z) dz`, with `J_P(z)=int_{-infty}^z b(u)q_P(u)du`. The otherwise possible affine Lebesgue-density ambiguity cannot be finite on the full real line. Since `b` is odd, `rho` and `p` are even, total `int b q_P=0`. Therefore `J_P(z)>0` for every finite real z, using b>0 on z<0 and the equivalent right-tail representation for z>=0. This is a **necessary** candidate bracket measure, not proof that feasible conductances exist.

Consequently every successful uniformly capped family satisfies the pointwise necessary bound

`4C >= g_P(z):=2J_P(z)/q_P(z)` for P=A and P=B and every real z.

In particular at z=0, `p(0)=5/8`, `1-p(0)=3/8`, and

`C >= (1/4) g_B(0) = J_B(0) / [2*(3/8)*rho(0)]`, 

`J_B(0)=int_0^infty (u³/3)*(1/2-exp(-u²)/8)*rho(u)du`.

The inequality is **strictly stronger than `C>1/2`**, because `1-p(u)>3/8` for all u>0 and `int_0^infty (u³/3)rho(u)du=rho(0)`. A high-precision **noncertified numerical diagnostic** yields `C>=0.6470219874711511267...`; this is not an outward interval theorem nor a practical sufficient C. Native and diagnostic per-label attempt clocks were not changed; `CN` is simply the permitted uniform total rate cap.

The conditional mean bracket at fixed z is also forced to be `g_A(z)` or `g_B(z)`, and its mixture recovers exactly `p g_A + (1-p)g_B=2`; this is consistent with BA exact parity compression and no contradiction. At z=0, unverified high-precision values are g_A~1.64714723007, g_B~2.58808794988. This necessity result rules out every cap `C<=1/2` even if all microscopic phases were feasible at higher C. **It does not solve the required all-N positivity/flow problem.**
