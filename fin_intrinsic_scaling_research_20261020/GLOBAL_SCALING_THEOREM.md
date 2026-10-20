# IS-001: a global polynomial scaling bound on the Hamming carrier

Date: 2026-10-20. Self-contained finite conditional mathematics, not a
claim of literature priority. Matrix scaling/source existence is inherited.

## 1. General theorem and its exact hypotheses

Let q>=2, N>=1 and K be the unmarked self-plus-single-coordinate reference
on E^N. Let pi>0 be supplied, and z>0 solve pi_x=z_x(Kz)_x. Suppose on
every oriented legal off-diagonal edge

    pi_x/pi_y<=C, C>=1.

Put

    D=N(q-1), D_out=(N-1)(q-1), d=1+D,
    Mbar=C+(C-1)D_out,
    B=1+D Mbar.

Then for every legal edge, both orientations included,

    Mbar^(-1)<=z_x/z_y<=Mbar,                        (1)
    P(x,x)=z_x^2/pi_x>=1/B,                         (2)
    v_x^2>=d/B, v_x=sqrt(d)z_x/sqrt(pi_x),
    v_x v_y>=d/B,                                  (3)
    P(x,y)>=B^(-1)sqrt(pi_y/pi_x), x~y.              (4)

No uniform independence or product form of pi is assumed. All coordinates
and labels are those of the supplied finite hypothesis, not physical cells.

## 2. Maximum-neighbour-ratio argument using actual cliques

Let M=max_(x~y) z_x/z_y. It is finite, attained and >=1, since both
orientations occur on a finite positive graph. Choose a maximizing edge
x=(a,r), y=(b,r), a!=b, differing in one fixed coordinate.

Let S=sum_c z_(c,r), the WHOLE q-state clique, including the self state.
Let H_x,H_y sum z over neighbours changing any OTHER coordinate. Exactly

    pi_x=z_x(S+H_x), pi_y=z_y(S+H_y).

For every other-coordinate background r', (a,r') and (b,r') are themselves
adjacent. Thus z_(a,r')>=z_(b,r')/M, whence H_x>=H_y/M. Moreover every
y' in H_y is adjacent to y, so H_y<=D_out M z_y. Since S>=z_x=M z_y,
H_y/S<=D_out. Therefore

    C>=pi_x/pi_y
       =M(S+H_x)/(S+H_y)
       >=1+(M-1)S/(S+H_y)
       >=1+(M-1)/(1+D_out).

This proves M<=C+(C-1)D_out. The same proof covers N=1 with H_x=H_y=0.
The self edge is counted once: replacing the source by a sum of N fully
marked update kernels would change this exact equation and is not allowed.

For each state x, its D off-neighbours satisfy z_y<=Mbar z_x. Consequently

    pi_x=z_x(Kz)_x<=z_x^2(1+D Mbar)=B z_x^2,

which proves (2). Equations (3),(4) follow by multiplying positive square
roots and using P_xy=sqrt(pi_y/pi_x)v_x v_y/d. The trivial global upper
v_x<=sqrt(d) remains true from the positive diagonal, but a constant upper
edge bound is neither proved nor needed for the local-mixing consequence.

For fixed q,C, Mbar=O(N), B=O(N^2), d/B=Omega(1/N). For example
B<=q[1+(q-1)(C-1)]N^2. When C=1 the sharper exact formulas give
Mbar=1, B=d, hence v_x=1 and P=K/d, as required.

## 3. Actual Gibbs application and the count warning

For the supplied A7 constant-diagonal all-pairs exponent, a legal labelled
move has exponent difference h_b-h_a using the same LOO cavity. Uniform
|h_j|<=H proves the required C=exp(2H), independent of N. Thus (1)--(4)
hold for its actual GLOBAL entropy optimizer, not a locally reselected law.

At count level pi_n contains N!/prod n_a!, and adjacent probability ratios
include n_a/(n_b+1); they are not uniformly bounded by exp(2H). Instead use
the paid RT orbit transport. If m_n=N!/prod n_a!, then
z_count(n)=sqrt(m_n)z_label(x), pi_count(n)=m_n pi_label(x), so the v is
unchanged on each orbit. Each count edge can be lifted to a legal labelled
edge. With Q the already derived SQRT count generator,

    P_count(n,n')=Q(n,n')v_n v_n'/d>=Q(n,n')/B.       (5)

All departure/reverse multiplicities stay inside Q and the Gibbs weights.
Count self probabilities equal labelled self probabilities because an
actual one-slot label change never preserves its count vector.

The new bound is polynomial, not the old order's stronger constant-a target.
That distinction remains in the claim register. The following local-window
theorem needs only this paid bound, not an unsupported improvement to it.
