# CH: incompatible next-event mechanisms despite causal forward laws

Date: 2026-10-21. Analytical model proof, same-author self-audit; not
formalization, independent review or physical validation.

## 1. Reused laws and the exact new question

The HYP source and RT transport have already proved, on E^N with q labels,

    K_xy=1 for x=y or one changed coordinate, otherwise 0,
    d=1+N(q-1), pi_g(x)=exp(g V(x))/sum_y exp(g V(y)),
    pi_g(x)=z_g(x)(K z_g)(x),
    P_g(x,y)=K_xy z_g(x)z_g(y)/pi_g(x).

The positive scaling is unique. Its count quotient is strongly lumpable;
one count move a->b has n_a microscopic continuations. These facts are
inherited, not a new maximum-entropy axiom, projection theorem or Gibbs
source. One unmarked self edge remains fixed.

For the supplied symmetric row-zero constant-diagonal A, remove only the
constant self contribution to the Gibbs potential and put

    V(x)=(1/N) sum_(i<j) A_(x_i,x_j).

This is exactly the same Gibbs law as exp[g n^T A n/(2N)], not a changed
interaction. The multinomial factor remains in the count law.

The separately authorized CC classification gives

    R_ab^s(h)=(1/q) exp[(h_b-h_a)/2+s(h_a+h_b-2 mean(h))],
    h_b=(g/N) sum_(j!=i) A_(b,x_j), a=x_i.

Every fixed real s is allowed. We will eliminate s rather than fit it.
We compare next actual changes, not HYP's discrete self probabilities to
CC's continuous holding times. No Poissonization is adopted.

## 2. A finite scaling susceptibility, with its invertibility payment

For any finite g the margin Jacobian in log z is diag(pi_g)+c_g, where
c_g(x,y)=K_xy z_g(x)z_g(y). For a nonzero real vector v,

    v^T[diag(pi_g)+c_g]v
       =(1/2) sum_xy c_g(x,y)(v_x+v_y)^2
       >=2 sum_x c_g(x,x)v_x^2>0.

Thus the finite smooth implicit-function theorem applies locally at each
g; uniqueness joins the local solutions. In particular z_g is analytic
near zero. This is the already paid positive-diagonal dual Hessian, not
an assumed spectral gap for a growing stochastic system.

At g=0, pi_0=q^(-N), z_0=(q^N d)^(-1/2). Set

    u_x=partial_g log z_g(x)|_0, Vbar=E_uniform V.

Differentiating the ACTUAL margin equation gives

    (d I+K)u=d(V-Vbar).                                  (1)

There is no missing derivative of the partition function. Its contribution
is Vbar; for our row-zero pair potential Vbar=0.

For completeness, consider products of centred one-coordinate functions
depending on exactly r distinct coordinates. A full sum over a dependent
coordinate is zero, so summing its q-1 different replacements gives minus
the original product. An unused coordinate contributes q-1 copies.
Consequently K has eigenvalue d-rq on this interaction-order subspace.
Products of constants and any basis of centred label functions span the
whole finite carrier. Hence d I+K has eigenvalues 2d-rq, all positive:
their smallest value is N(q-2)+2. This independently pays (1).

Each row-zero pair term lies in order r=2, so

    u=c_N V, c_N=d/[2(d-q)], N>=2.                       (2)

In particular c_2=23/22 and c_3=17/22 at q12. These are finite-state
gain susceptibilities, NOT a critical drift or a physical constant.
At a fixed present x and two legal different destinations y,y',

    partial_g log[P_g(x,y)/P_g(x,y')]|_0
       =c_N[V(y)-V(y')].                               (3)

The current-state denominator cancels exactly. A common clock, even one
depending on x and g, cannot change (3).

## 3. A parameter-free count identity for EVERY CC member

Now N=3. By symmetric circulant structure write

    a0=A_00, a1=A_01, a2=A_02, Delta=a0-a1.

At x=001, distinguish the following count destinations:

    A: 000, from the one 1-slot changing to 0;
    B: 002, from the one 1-slot changing to 2;
    C: 012, from either of the TWO 0-slots changing to 2.

Let J_A,J_B,J_C be their conditional probabilities at the next count
change; all other destinations are retained in the same normalized row.
For any positive rates define r_A=Q_A, r_B=Q_B, r_C=Q_C/2.
Then r_C/sqrt(r_A r_B)=J_C/[2 sqrt(J_A J_B)].

The LOO fields for the 1-slot are h_b=2g A_b0/3. For a 0-slot they are
h_b=g(A_b0+A_b1)/3. Their means vanish by the row-zero property. Direct
substitution into the classified CC rates gives

    log r_A= -log q + (g/3)[a0-a1+2s(a0+a1)],
    log r_B= -log q + (g/3)[a2-a1+2s(a2+a1)],
    log r_C= -log q + (g/6)(a2-a0)
                          +(gs/3)(a0+a2+2a1).

The neutral normalization, all s terms and any common row clock cancel:

    F_CC(g):=log[J_C/(2 sqrt(J_A J_B))]+g Delta/3=0,     (4)
    J_C^2=4 exp(-2g Delta/3) J_A J_B.

This is exact for every g and every real s, even if a different s were
chosen at each g. No selector value of s has been inferred. The identity
also retains the real count multiplicity 2, not a marked null attempt.

## 4. The ACTUAL HYP source violates that identity near zero

Let z_A,z_B,z_C be the labelled scaling at representatives 000,002,012.
Permutation symmetry makes the two microscopic C continuations equal.
The quotient formula therefore gives exactly

    F_HYP(g)=log z_C-(log z_A+log z_B)/2+g Delta/3.

At zero F_HYP(0)=0. The three supplied pair potentials are

    V_A=a0, V_B=(a0+2a2)/3, V_C=(2a1+a2)/3,
    V_C-(V_A+V_B)/2=-2 Delta/3.

Equation (2) now gives at q12,N3

    F'_HYP(0)=(1-2c_3)Delta/3=-(2/11)Delta.             (5)

The true A7 is PSD and has the parity eigenvector sigma with eigenvalue
Lambda6>0. For v=e0-e1 its squared parity projection is 4/12=1/3.
All other spectral contributions are nonnegative. Therefore

    2 Delta=v^T A7 v >=Lambda6/3,
    Delta>=Lambda6/6>0.                               (6)

This sign uses the supplied exact PSD/parity premises, not the rounded
safe_seed matrix. Analyticity and (5) prove that some delta>0 exists such
that F_HYP(g)<0 for EVERY 0<g<delta. No useful numerical delta is claimed.
Combining with (4), for those gains NO CC parameter, and NO common clock
at this very state, reproduces the HYP next-count-change law.

This is a failure of simultaneous mechanism compatibility, not a proof
that either independent hypothesis is mathematically inconsistent. It
does not establish the sign at g_c or at every positive g. A first-order
expansion has not been silently extrapolated to critical coupling.

For a general symmetric circulant row-zero q>=3 matrix with Delta>0,
the same calculation gives F'_HYP(0)=-q Delta/[6(q-1)].
The actual q12 result above is its stated specialization.

## 5. A bounded decision contrast and uniform separation over s

Avoid singular log estimates by defining, on a next-event probability row J,

    D_g(J)=J_C^2-4 exp(-2g Delta/3) J_A J_B.             (7)

For g>=0 and Delta>0, D is bounded between -1 and 1. Every CC row has
D=0 exactly. At g=0 the HYP row has

    J_A=J_B=1/33, J_C=2/33.

Using (5) in (7), rather than discarding normalization derivatives, yields

    partial_g D_g(J_HYP(g))|_0=-16 Delta/11979.          (8)

Here is a complete TV payment. In two independent mathematical copies
of the SAME next-event experiment, define

    f_g(U,V)=1_(U=C,V=C)
             -2 exp(-2g Delta/3)[1_(U=A,V=B)+1_(U=B,V=A)].

Its expectation is D_g(J), and its range has length at most 3. For any
two event laws J,J', expectations differ by at most
3 TV(J tensor J,J' tensor J')<=6 TV(J,J'). The last inequality follows
by inserting J' tensor J and coupling one factor at a time. This is a
mathematical product experiment, not an assumed physical resetting device.
Consequently

    inf_(s in R) TV(J_HYP(g),J_CC,s(g)) >= |D_g(J_HYP(g))|/6,
    liminf_(g->0+) inf_s TV(J_HYP(g),J_CC,s(g))/g
           >=8 Delta/35937>0.                         (9)

Only one next-count-change law is being compared; hidden histories or a
new apparatus are not needed for the MODEL statement. This is neither a
critical-process theorem nor a finite-N resource or empirical bound.

## 6. What this says, and does NOT say, about causality

Both hypotheses already define causal forward laws in their stated model
categories. For HYP, prefix probabilities are pi(x0) times successive
P(x_i,x_(i+1)); summing a future endpoint uses P1=1. The conditional next
state given a finite past is P evaluated at its current state. The global
static optimization selecting P does NOT consult realized future events.
For CC, a finite state's holding law and next-event probabilities are
determined by its current supplied rates. Its Gibbs conditional law and
preparation remain inputs; there is no future-event argument in a rate.

This is model nonanticipation, not an empirical intervention theorem.
Varying g here compares different autonomous laws; it is not a future
control protocol and cannot establish backward causation. Reversibility
of either stationary law also does not mean an absence of causal order.

Thus causal order alone cannot select which mechanism is nature's law.
The new result instead excludes their automatic conjunction: a whole-law
entropy selection and an exact local additive-channel response need not
describe the same causal branching. Correlations/equilibrium alone would
not reveal that distinction. No third principle is adopted as a repair.

F01 still asks for an independently grounded restriction on the actual
evolution. The nadsoliton remains primordial informational whole, with
observer emergent; no lower substrate, physical clock/spacetime, apparatus,
natural A7, source closure, QM/GR or ToE has been manufactured here.
