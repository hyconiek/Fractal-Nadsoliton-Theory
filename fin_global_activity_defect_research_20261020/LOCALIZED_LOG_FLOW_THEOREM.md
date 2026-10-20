# AD-002: local control of exact flows without dividing a global error

Date: 2026-10-20. Self-contained finite positive-flow identity and its
application. The trial's margins are unequal and the boundary is retained.

## 1. Exact local identity

Let c_xy=K_xy z_x z_y be the exact source flow with margins pi_x.
Let u_xy=K_xy z_app,x z_app,y have margins
tilde_pi_x=pi_x(1+delta_x), and put R_x=z_x/z_app,x.
Restrict sums to allowed edges, including self edges. Then
log(c_xy/u_xy)=log R_x+log R_y. For any nonempty subset D,

    S_D=(1/2)sum_(x,y in D)(c_xy-u_xy)log(c_xy/u_xy)
       =sum_(x in D)(pi_x-tilde_pi_x)log R_x
          -sum_(x in D,y outside D)(c_xy-u_xy)log R_x. (1)

This follows by symmetry of the internal flows and their ACTUAL row sums.
It is not a conditional optimizer identity. No equality of margins or
global-to-local normalization shortcut is assumed.

For positive t, (t-1)log t>=4(sqrt(t)-1)^2. One proof differentiates
log s-2(s-1)/(s+1), whose derivative is (s-1)^2/[s(s+1)^2], and then
sets t=s^2. Hence

    S_D>=2 sum_(x,y in D)(sqrt(c_xy)-sqrt(u_xy))^2.  (2)

## 2. Pay the true boundary and small local mass

IS's actual GLOBAL bounds give c1 N^-1/2<=v_x<=c2 N^1/2. Since a_N
is uniformly bounded above/below, R_x=v_x/a_N(x) satisfies

    sup_x |log R_x|<=C L_N, L_N=1+log(N+1).

From AD-001, |delta_x|<=C N^-2 and tilde_pi_x<=2pi_x eventually.
Let B_D be the states in D with a legal departure. The cross-flow term
in (1) has absolute value at most

    C L_N sum_(x in B_D,y outside)(c_xy+u_xy)
       <=C L_N pi(B_D),

because exact rows sum to pi_x and auxiliary rows sum to tilde_pi_x.
The margin term is at most C L_N N^-2 pi(D). Therefore

    H_D^2:=sum_(x,y in D)(sqrt(c_xy)-sqrt(u_xy))^2/pi(D)
       <=C L_N[N^-2+pi(B_D)/pi(D)].                 (3)

For the actual MT domain the second ratio is eta_N=pi_D(B_D), bounded
by polynomial*exp(-Nb). This proof never divides a GLOBAL O(N^-2)
error by the exponentially small pi(D). It establishes the local estimate
directly, using the real departure layer. A tiny complete-graph negative
control retains why the omitted cross term would be invalid.

## 3. Weighted flow variations, not pointwise scaling

For any nonnegative weight w_x, Cauchy--Schwarz and the row bounds give

    sum_(x,y in D) w_x |c_xy-u_xy|/pi(D)
       <=sqrt(H_D^2)
          *sqrt(sum_(x,y in D)w_x^2(sqrt(c_xy)+sqrt(u_xy))^2/pi(D))
       <=C sqrt(E_piD w^2) sqrt(L_N[N^-2+eta_N]).    (4)

This is valid also for w=1+|Z_N| because actual MT moments bound its
second moment. An unweighted L1 estimate multiplied by the worst-case
N^(1/4) amplitude would be insufficient; (4) pays that weight directly.

The result controls averaged INTERNAL edge flows of the exact global
source relative to the derived trial. It does not assert v_x->1 uniformly,
close complete path laws, a constant edge lower bound or equal posterior
information. Both exact and auxiliary objects remain mathematically typed.
