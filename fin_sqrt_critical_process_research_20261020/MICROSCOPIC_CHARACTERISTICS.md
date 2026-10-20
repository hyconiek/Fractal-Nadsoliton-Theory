# SQ-001: exact characteristics and the paid fast selection

Date: 2026-10-20. Use [CONTRACT.md](CONTRACT.md), keeping the original
all-pairs/LOO convention and all ten fast count directions.

## 1. Actual jumps and finite-N cavity error

A move a->b changes m by (sigma_b-sigma_a)/N. For the full SQRT parent,

    Q_N m=sum_ab p_a R_ab^N(p)(sigma_b-sigma_a)=b_N(p),
    N Gamma_Q(m,m)=sum_ab p_a R_ab^N(p)(sigma_b-sigma_a)^2=B_N(p).

Terms a=b vanish. Put R_ab(p)=exp((g/2)[(Ap)_b-(Ap)_a]) and define
b(p,g),B(p,g) by replacing R^N with R. With H=G M_A,

    R_ab^N=R_ab exp[(g/(2N))(A_aa-A_ba)],
    |R_ab^N-R_ab|<=H exp(2H)/N.

Indeed |g(A_aa-A_ba)/(2N)|<=H/N and R_ab<=exp(H) on the simplex.
Summing the actual per-departure weights p_a gives uniform bounds

    |b_N-b|<=2(q-1)H exp(2H)/N,
    |B_N-B|<=4(q-1)H exp(2H)/N.                      (1)

No native HB rate or covariance formula has been substituted.

## 2. Exact limiting parity line and the new clock factor

On p_m, (Ap_m)_a=chi m sigma_a. For k=g chi, a copy of sign s has
q/2 possible opposite labels, each at rate exp(-s k m). Hence exactly

    b(p_m,g)=q[sinh(km)-m cosh(km)],                  (2)
    B(p_m,g)=2q[cosh(km)-m sinh(km)].                 (3)

In particular B(u,g)=2q. The line formula is a limiting drift, not exact
finite-N parity lumpability; (1) retains the actual cavity correction.
For k near 1 on the fixed MT cylinder, Taylor's formula with bounded
exponential derivatives gives

    b(p_m,g)=q[(k-1)m-m^3/3]
                +O(|k-1||m|^3+|m|^5).              (4)

The linear coefficient is exact; the cubic coefficient
k^3/6-k^2/2 differs from -1/3 by O(|k-1|). No kinetic clock is refitted.

## 3. Smooth symmetry removes the entire first fast variation

Let T shift the q labels by two. A commutes with T, sigma is fixed by T,
and p_m is fixed by T. Thus b(Tp,g)=b(p,g). For w perpendicular to both
1 and sigma,

    (1/6)sum_(j=0)^5 T^j w=0.

Since b is smooth, its derivative at p_m is a linear invariant covector.
Differentiating the symmetry and averaging yields exactly

    D_w b(p_m,g)[w]=0.                               (5)

This is a consequence of the actual SQRT smooth function, not an
assumption of fast independence, parity Markovity or a general statement
for nonsmooth kinetics. Second-order Taylor then pays

    |b(p_m+w,g)-b(p_m,g)|<=C_w ||w||^2.               (6)

For example, with J=G sqrt(q) M_A, a valid conservative constant is
C_w=q(q-1)exp(H)(2J+J^2). To see this, set
a_ab=(g/2)(A_b-A_a) as a row covector. Each scalar edge p_a exp(a_ab p)
has second variation exp(a_ab p)[2w_a(a_ab w)+p_a(a_ab w)^2].
Use |a_ab w|<=J||w||, p_a<=1, |sigma_b-sigma_a|<=2 and the one-half
Taylor remainder on the segment p_m+theta w inside the simplex.
There can be a NONZERO quadratic fast term; the verifier retains one.

Uniform first derivatives of B on the fixed compact region also give

    |B(p,g)-2q|<=C(|m|+||w||),                       (7)

because A*u=0 makes B(u,g)=2q for every bounded g. Equations (1),(4),(6)
and (7), not a shared Gibbs law alone, supply the drift/noise limit.

## 4. Correct scaled and reflected characteristics

The critical space/time scaling gives

    beta_N=N^(3/4)b_N,
    gamma_N=N Gamma_Q(m,m)=B_N,
    |Delta Z_N|<=2N^(-3/4).

The accelerated absolute third jump characteristic is bounded by
2N^(-3/4)gamma_N<=8R N^(-3/4); the fourth one is bounded by
4N^(-3/2)gamma_N. These are actual jump sums, not fitted Gaussian noise.

Let B_D be the true MT departure layer. Suppressing outgoing edges changes
beta by at most 2R N^(3/4)1_B_D and gamma by at most 4R 1_B_D.
The reflected bracket is nonnegative and bounded by 4R everywhere.
These corrections are paid in expectation below; the full-generator
martingale is not falsely declared exact for a reflected chain.
