# MM-002 — Generator-derived fast characteristics

Write alpha=k/N and f_r(j)=sqrt(2)cos(2*pi*r*j/12), r=1,3. For a predetermined block B define a_a=N/k-1 if a in B, and -1 otherwise. Then F_r=N^-1/2 sum_a f_r(x_a), U_r=N^-1/2 sum_a a_a f_r(x_a), V_r=F_r+U_r. For the AUX-HB-OU labelled generator, q_j(theta)=softmax_j(X theta), whose proposal does *not* depend on the refreshed slot's old label. On any unrestricted count configuration (or reflected interior),

    L F_r = sqrt(N) sum_j f_r(j) q_j(theta) - F_r,
    L U_r = -U_r,                           (EXACT)

because sum_a a_a=0. Reflection rejects whole proposals whose total n leaves D_N; it does NOT renormalize q. All boundary terms must be retained as O(sqrt(N)*Pi_D(boundary)) in stationary L1 (fixed alpha), assuming the inherited exponential local-Gibbs departure bound.

Take an orthonormal Fourier X representation of A7: in the active cosine r=3 direction X_{j,3}=sqrt(Lambda3/12) f3(j). Let c=g_c*Lambda3/12=Lambda3/Lambda6 in (0,1), gamma=1-c, and H_3=sqrt(N)*sqrt(Lambda3/12)*theta_3. The exact field drift is

    L H_3 = kappa (c_N F_3-H_3),    c_N=g_N Lambda3/12,
    d<M_H3>/dt = 2 kappa c_N.

At critical lambda=0, c_N=c. Cross jump/continuous brackets between H and F,U vanish identically. In the null r=1 direction A7 has Lambda1=0, so no H1 drives F1.

The exact predictable count-jump brackets are, with w_a(F)=1 and w_a(U)=a_a,

    d<M_{X_r},M_{Y_s}>/dt
      =N^-1 sum_a w_a(X)w_a(Y) sum_j q_j(theta)
        [f_r(j)-f_r(x_a)][f_s(j)-f_s(x_a)].

Here X,Y are F or U. For r=s the inner square is E_q f_r(new)^2 + f_r(old)^2 - 2 f_r(old) E_q f_r(new), NOT identically 2 statewise. The N=12 independent jump fixture with block indices all even, uniform counts and theta=0 gives FF=2, UU=2 and FU=1: mixed brackets vanish **only in the stationary/conditional-block limit**, not as exact finite-state identities.

On the admitted MT local positive cylinder, the parity coordinate has theta6=O_Lp(N^-1/4) while all transverse field coordinates have O_Lp(N^-1/2). Parity-axis symmetry implies the f1/f3 mean proposals vanish at every pure soft configuration. A Taylor expansion around that axis, bounded softmax derivatives, and stationary Hölder estimates give for fixed T

    E_Pi integral_0^T |sqrt(N) f3·q(theta)-H3| dt
       <= C_T [N^-1/4 + N^-1/2 + |c_N-c|],

and the corresponding f1 term tends to 0. Conditional hypergeometric block moments and full stationary local moments yield, in integrated L1,

    <M_F3> -> 2t,    <M_U3> -> 2 c_alpha t,
    <M_F3,M_U3> -> 0,    <M_H3> -> 2 kappa c t,

with zero cross-mode limits. The same local estimates pay bracket errors O(N^-1/4)+O(N^-1/2) (not a claimed universal finite-N bound). Gaussian theta tails are controlled by its known conditional normal rho_n and stationarity, not dropped. For each fixed alpha away from zero, jumps of F/U are O(N^-1/2); the OU field has continuous Brownian martingale part. These are generator calculations plus stated inherited stationary moment hypotheses, not a guessed two-OU covariance.
