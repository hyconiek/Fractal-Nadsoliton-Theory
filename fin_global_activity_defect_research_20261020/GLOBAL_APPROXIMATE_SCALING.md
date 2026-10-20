# AD-001: derive a globally uniform O(N^-2) margin residual

Date: 2026-10-20. Conditional smooth finite-model analysis, not a numerical
matrix scaling or an assumed replacement law.

## 1. The exact equation to approximate

Let epsilon=1/N, r=q-1 and d_N=N r+1. For a!=b define

    delta_ab=e_b-e_a,
    E_ab(p,g)=(g/2)[(Ap)_b-(Ap)_a],
    c_ab(p,g)=p_a exp(E_ab),
    ell_ab(g)=(g/2)(A_aa-A_ba).

The true SQRT comparison count rate is exactly
Q_N(n,n+delta_ab)=N c_ab(p,g)exp(epsilon ell_ab). The last exponential
is the LOO correction; it is not optional. For a positive trial function
a(p,g), set z_app(n)=sqrt(pi_N(n))a(p,g)/sqrt(d_N). Its row margin is

    z_app(n)(K z_app)(n)=pi_N(n) F_N[a](p),
    F_N[a](p)=a(p)/(r+epsilon)
       *[epsilon a(p)+sum_(a!=b)c_ab exp(epsilon ell_ab)a(p+epsilon delta_ab)]. (1)

This identity includes the diagonal 1 and count multiplicities. Terms with
n_a=0 vanish; every remaining p+epsilon delta is a legal simplex point.

## 2. A derived first correction

Write

    tau=sum c_ab, b_vec=sum c_ab delta_ab, kappa=sum c_ab ell_ab,
    a0=sqrt(r/tau),
    a1=[1-(1+kappa)a0^2-a0 b_vec dot grad a0]/(2 tau a0),
    a_N=a0+epsilon a1.                              (2)

All functions use the supplied A,g; no target diffusion coefficients enter.
The derivative is tangent to sum p=1; b_vec is tangent, so the expression
is independent of a smooth normal extension.

For |g|<=G, H=G max|A_ij| and p in the closed simplex,
r exp(-H)<=tau<=r exp(H). All functions in (2), with their fixed derivatives,
are uniformly smooth on a fixed neighbourhood of that compact simplex
and gain interval. In particular a0 is bounded above/below and a_N stays
positive with uniformly bounded fixed derivatives for all sufficiently
large N. No division by p_a or min pi_N occurs.

Taylor's formula on the legal edge segment gives uniformly

    sum c_ab exp(epsilon ell_ab)a_N(p+epsilon delta_ab)
      =tau a_N+epsilon[kappa a_N+b_vec dot grad a_N]+O(epsilon^2).

Substitution into (1) gives numerator

    tau a0^2
      +epsilon[2 tau a0 a1+(1+kappa)a0^2+a0 b_vec dot grad a0]
      +O(epsilon^2)
    =r+epsilon+O(epsilon^2).

Thus the ACTUAL full-domain trial margins obey

    delta_N(n):=F_N[a_N](p)-1,
    sup_all_counts |delta_N(n)|<=C/N^2.              (3)

The uniform remainder uses bounded second derivatives of a0, first/further
fixed derivatives of a1 and exp(epsilon ell); finite q bounds the sum.
Counts on simplex faces cause no missing departure term or singularity.

At u, A*u=0 and circulant symmetry give tau(u)=r,
D_T tau(u)=0, b_vec(u)=0 and kappa(u)=g tr(A)/2. Hence
a0(u)=1 and a1(u)=-kappa(u)/(2r). The leading trial a0 by itself has
a nonzero O(1/N) row residual there for g>0 and nonzero PSD A. Omitting
the cavity kappa or transport-gradient term in (2) would not pay (3).
This is the single motivated refinement of the leading approximation.

## 3. Conservative analysis operator and local expansion

Let u_xy=K_xy z_app(x)z_app(y), a positive symmetric auxiliary flow.
Its margins are pi_N(x)(1+delta_N(x)), NOT pi_N exactly. In particular one
must not apply the same-margin entropy-deficit theorem to it as if it were
an admissible stationary optimizer.

Its off-diagonal transition table u_xy/pi_x is usable for comparison.
If needed, make it conservative by choosing self probability

    P_app(x,x)=1-sum_(y!=x)u_xy/pi_x
              =a_N(x)^2/d_N-delta_N(x)>0

eventually: the first term is uniformly at least c/N and (3) is O(N^-2).
This is a reversible ANALYSIS kernel with pi_N, not an adopted physical
law or the exact entropy source. Only its off-diagonal entries will enter
the characteristic comparison.

For the actual row-zero A7,

    tau(p)=(sum_a p_a exp(-g(Ap)_a/2))
              (sum_b exp(g(Ap)_b/2))-1.

Its tangent derivative at u vanishes because A*1=0. Therefore, uniformly
on the simplex/gain interval,

    |a_N(p)-1|<=C(||p-u||^2+N^-1),
    |a_N(p+delta_ab/N)-a_N(p)|<=C/N.                 (4)

These are properties of the trial function, not pointwise convergence of
the unknown exact v_N. Its actual comparison is paid by the next lemma.
