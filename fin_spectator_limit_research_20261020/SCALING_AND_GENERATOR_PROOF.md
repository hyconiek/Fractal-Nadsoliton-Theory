# SL-001: finite-source scaling and the derived generator

Date: 2026-10-20. Conditional self-audited mathematics, not peer review,
mechanization or a new physical law. Existence/uniqueness of positive
scaling is inherited from the hypothesis and coercive-dual proof.

## 1. Primary finite source and a covariance guard

The supplied spectator family has R_c=K+cI. To track a reference-gauge
change correctly, also consider the mathematical pair (K,D), where
D=diag(d_x), d_x>0 is fixed. The primary source is still D=I. Its scaling
and kernel satisfy

    pi_x=z_c,x[(Kz_c)_x+c d_x z_c,x],
    P_c,xy=(K_xy+c d_x 1_(x=y))z_c,x z_c,y/pi_x.        (1)

Let H=D^(-1/2)KD^(-1/2). Reference-gauge invariance gives exactly the same
kernel from H+cI with scaling sqrt(D)z_c. Therefore it suffices to prove
the primary theorem with K replaced by H. This is coordinate bookkeeping,
not permission to choose another D to obtain a preferred law.

Write s_x=sqrt(pi_x), v_c,x=sqrt(c d_x)z_c,x and

    A=max_x (Hs)_x/s_x < infinity.

Every entry is fixed and pi is positive; no numerical gap or partition
function estimate is being assumed. Equation (1) becomes

    pi_x=v_c,x^2+c^(-1)v_c,x(Hv_c)_x.                  (2)

## 2. Paid uniform scaling error

All terms in (2) are nonnegative, so 0<v_c,x<=s_x. Define r_c,x=v_c,x/s_x.
Then

    0<=1-r_c,x
       =v_c,x(Hv_c)_x/[c pi_x(1+r_c,x)] <=A/c.         (3)

The upper bound uses Hv_c<=Hs and v_c<=s; it is valid even when A/c>1
and then merely loose. For c>=2A it also gives r_c,x>=1/2. No state is
dropped because its Gibbs weight is small; positivity and finite A are
actual fixed-carrier premises.

## 3. Derived conservative rate law

Set Q_c=c(P_c-I). Its off-diagonal entries are

    Q_c,xy=H_xy sqrt(pi_y/pi_x) r_c,x r_c,y, x!=y.

Q_c is an auxiliary Poissonized embedding for estimates. The actual finite
observed parent in this packet still has kernel P_c on its intrinsic grid;
its path law is not the Q_c continuous-time path law.

Hence the **derived**, not inserted, limit is

    Q_xy=H_xy sqrt(pi_y/pi_x), x!=y,
    Q_xx=-sum_(y!=x)Q_xy.                              (4)

The diagonals of Q_c are also minus their outgoing sums, by exact source
stochasticity. Keeping only c d_x z_c,x^2≈pi_x in the source diagonal
would lose this escape term; it is not a valid diagonal limit.

Both Q_c and Q are conservative finite-state generators, with outgoing
rates <=A. Detailed balance follows from
pi_x Q_xy=H_xy sqrt(pi_x pi_y), and likewise for Q_c. Thus pi is invariant.
For the induced row-sum norm ||B||_infinity=max_x sum_y |B_xy|,

    ||Q_c||_infinity, ||Q||_infinity <=2A,
    ||Q_c-Q||_infinity <=4A^2/c.                       (5)

Indeed 0<=1-r_x r_y<=2A/c; the off-diagonal row error is at most
(2A/c)sum_(y!=x)Q_xy<=2A^2/c, and the diagonal adds the same mass.
No commutativity of Q_c and Q is assumed.

For the actual completion D=I, H=K. The square-root-rate shape is a known
conditional form recorded in the old MaxCal/source audit. The new result
is its paid derivation from the frozen finite spectator family, with errors
and path law, not its general mathematical discovery.

## 4. A complete joint-gauge convention

If K'_xy=b_x K_xy b_y, the same finite completion requires
D'_xx=b_x^2 d_x. Then H'=H, P'_c=P_c for every c, and Q'=Q. Changing only
K while leaving D unchanged is **not** that same source family, even when
the selected law at c=0 was gauge equivalent.

For example, use a three-state complete reference K=ones and uniform pi.
Let b=(1,2,3). With the correctly transformed D'=diag(1,4,9), all limiting
off-diagonal rates are still 1. Holding D=I instead gives rates 2,3,6 on
the three undirected edges: not one common clock multiple. This is a
non-adopted completion countercontrol, not a new model fit.

Thus the conditional source of activity is the supplied pair (K,D), or
equivalently H. A reference class at one isolated level does not by itself
source a completion. The primary combinatorial spectator construction
already supplies D=I; it is not independently selected as nature's law.
