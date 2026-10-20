# FH-003: a paid native diagonal limit, with its enormous history cost

Date: 2026-10-20. Conditional on the unchanged HYP/IS/MT/AD mathematical
premises, including the AD self-audited critical theorem. This is not an
independent audit of AD or a change to primary FIN.

## 1. Quantitative continuity on the actual labelled Hamming carrier

Let m=q^N, d_N=1+N(q-1), and suppose adjacent pi ratios lie in [1/C,C],
C independent of N. Set delta=1/(L pi_min), and assume delta<=1/2. Largest-
remainder occupation rounding gives |pi^L/pi-1|<=delta at every state.
Along pi(t)=(1-t)pi+t pi^L the adjacent ratios are at most C'=3C. Use the
already proved IS bound at that very pi(t):

    B=1+N(q-1)[C'+(C'-1)(N-1)(q-1)],
    P_t(x,x)>=1/B,
    c^L_xy>=mu_0=pi_min/(2B sqrt(C')) on all allowed entries. (1)

The last claim uses IS P_t(x,y)>=sqrt(pi_t(y)/pi_t(x))/B on off-edges;
diagonal entries have the stronger lower bound pi_t(x)/B.

Write u=log z(t). Differentiating pi=z(Kz) exactly gives

    (I+P_t) u_dot=pi_dot/pi(t).

With s=1/B and R_t=(P_t-s I)/(1-s) stochastic, the convergent alternating
Neumann series yields

    ||(I+P_t)^(-1)||_(infinity->infinity)<=1/(2s)=B/2.

Since ||pi_dot/pi(t)||_infinity<=delta/(1-delta),

    max_allowed |log(P^L_xy/P_xy)|
       <=theta=(B+1)delta/(1-delta).                  (2)

For theta<=1, each row TV is at most (exp(theta)-1)/2<=theta. Initial
conditional laws pi^L(.|D), pi(.|D) differ by at most delta/(1-delta).
Sequential coupling therefore gives the safe bound

    TV(Q^L_[0:k](.|D),Q_[0:k](.|D))
       <=2delta+2k theta.                            (3)

No stationary mean rate or unspecified Markov closure replaces this full
history calculation. Rational rounding may break exact finite-N slot
symmetry; (2),(3), not a falsely claimed exact symmetry, control that error.

## 2. A concrete exact-integer cycle length

For the actual labelled Gibbs input with |g_N|<=G and max|A7_ij|<=M,

    pi_N,min>=q^(-N) exp(-G M N),
    C=exp(2GM)

are conservative global bounds: the exponent (g_N/(2N))n^T A7 n has
range at most GMN, and a legal labelled LOO move changes it by at most
2GM. The gain bound and A7 are supplied, not derived here.

Freeze **L_N=q^(N^2)**, an exact integer, and round occupations on the
labelled carrier as above. This is an existence construction, not an
enumeration, feasible resource recommendation or new choice of A7.
For all sufficiently large N, delta_N->0, theta_N->0 and L_N mu_0>=1.
Using mu_0 in FH-002 gives

    A_N <= (m^2+m-1)log(L_N+1)
              +m^2(m-1)/(mu_0 L_N)+log m.

For any nonempty initial count event D_N, its labelled probability is at
least pi_N,min. Combining the exact preparation cost and (3),

    TV(U_n^[k](.|D_N),Q_N,[0:k](.|D_N))
      <=sqrt(A_N/[2 pi^L_N(D_N) floor(L_N/k)])
          +2delta_N+2k theta_N,                     (4)

where pi^L_N(D_N)>=(1-delta_N)pi_N(D_N). For every fixed polynomial
k_N<=N^r the right side tends to zero. To check all growth terms, its
leading squared term is bounded by a polynomial times

    exp[(2 log q+log q+GM)N-(log q)N^2],

and the second A_N term has exponent at most
(3 log q+2(log q+GM))N-2(log q)N^2, again times a polynomial.
Both vanish; (2),(3) also do. Constants depend on fixed q,G,M,r, not on
an unproved uniform convergence of finite-dimensional implicit maps.

## 3. An actual critical-path corollary, not just iterated limits

Set q=12, N in12N, exact critical gain g_c+lambda/sqrt(N), fixed bounded
lambda, and D_N the actual MT local count domain. Condition the ROOT of
the finite global occupation ensemble to lie in D_N. This is preparation,
not a claim that the global Gibbs law is already local critical equilibrium.

Take k_N=ceil(d_N sqrt(N) S) for fixed S. Equation (4) compares the ENTIRE
labelled discrete record through that critical horizon with the original
unreflected HYP source started from pi_D,N. Under that source, comparison
with its stationary reflection costs at most k_N pi_D,N(B_D), which
vanishes by the accepted exponential departure bound. The inherited AD
theorem therefore yields, under this specified finite ensemble,

    N^(1/4)m(X_floor(d_N sqrt(N)s)) => Z(s),
    dZ=12(chi lambda Z-Z^3/3)ds+sqrt(24)dB,
    Z0~rho_lambda, 0<=s<=S, in D([0,S]),J1.           (5)

Thus the same collective limit has a finite admissible-history realization
without a new noise law or additional spectator coordinates. This closes
the growing-carrier limit-order payment for THIS exceedingly expensive
sequence, not its optimal cost. Uniform counting of count words, a short
cycle, an arbitrary preparation, native HB compatibility, real-time units
or a physical interpretation of L_N do not follow.

## 4. Infinite-time warning and F01

For every finite L, extend a random cycle periodically forever. The shift
is a finite deterministic permutation, has entropy rate zero, and repeats
exactly. For a nontrivial connected selected P, every row has at least two
positive entries and max_xy P_xy<1. Any specified infinite periodic word
then has probability zero under P. The countable union of all periodic
words still has probability zero. Hence the finite-cycle extension and
the infinite Markov law have **infinite-history TV exactly one**, for
every finite L, despite the paid finite-window limit.

No finite-entropy eternal source obstruction is evaded. Limits in history
length, observation horizon, N and entropy rate cannot be interchanged.
The construction explains conditional finite-record typicality from an
already supplied counting reference and Gibbs occupation input. It does
NOT explain why nature chooses those inputs or typicality itself. Strict
natural F01/F03/F06 remain OPEN_SOURCE; no law is promoted to QM/GR/ToE.
