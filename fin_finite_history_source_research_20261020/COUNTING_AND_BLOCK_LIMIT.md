# FH-002: a direct counting proof and whole-record bound

Date: 2026-10-20. Natural logarithms. This is a specialized elementary proof
of a known Markov-conditioning mechanism, not a claim of its discovery.

## 1. A balanced integer table supplies enough actual words

First suppose pi=n/L exactly, with every n_x>0. Put

    mu=min_(K_xy=1)c*_xy>0,
    a_m=m^2+m-1, b_m=m^2(m-1),
    A_L=a_m log(L+1)+b_m/(mu L)+H(pi).                (1)

Assume L mu>=1. For x<y on legal off-diagonal edges choose
M_xy=M_yx=floor(L c*_xy), and set

    M_xx=n_x-sum_(y!=x) M_xy.

This has exact row AND column sums n, nonnegative entries and legal support.
All allowed entries are positive: off-diagonal floors are at least one,
and M_xx>=L c*_xx>=1. It is connected. Write c^M=M/L. Off-diagonal errors
have absolute value at most 1/L, diagonal errors at most (m-1)/L, whence

    sum_xy(c^M_xy-c*_xy)^2<=b_m/L^2,
    KL(c^M||c*)<=sum_xy (c^M_xy-c*_xy)^2/c*_xy
                    <=b_m/(mu L^2).                 (2)

The KL inequality follows from log u<=u-1 and equal total masses. The
inherited exact entropy identity gives F(c*)-F(c^M)=KL(c^M||c*).

Fix a root a and an in-arborescence t_x toward a using positive M edges.
At each vertex arrange its outgoing destinations as a list with the given
multiset M_x. At x!=a require its last destination to be t_x; root a has
no last-destination restriction. The number of list choices is exactly

    product_x [n_x!/product_y M_xy!] *
      product_(x!=a) M_(x,t_x)/n_x
      >=L^(-(m-1)) product_x [n_x!/product_y M_xy!].   (3)

Following these lists starting at a consumes every entry. If it first
stalls, balance forces its endpoint to be a. Its used walk is then closed,
so residual incoming/outgoing counts balance at every vertex. If any list
is unused, its last tree edge is unused. Following these residual tree
edges forces unused incoming and hence unused outgoing all the way to a,
contradicting its exhausted list. Thus there is no unused list. Reading
successive departures recovers every row list from the produced word, so
different lists give different rooted words. Parallel equal destinations
are not given artificial edge labels.

For completeness, the elementary multinomial lower bound is

    n!/product_j k_j! >=(n+1)^(-m) exp[n H(k/n)].      (4)

Indeed under the multinomial law with probabilities k_j/n, its own count
vector k is a mode: if another positive-probability vector has a surplus
in i and deficit in j, transferring one unit i->j raises its probability
by l_i k_j/[(l_j+1)k_i]>1. Iterate to k. At most (n+1)^m count vectors
sum to probability one, so the mode has probability at least (n+1)^(-m).
Zero k_j prohibit positive l_j and cause no exception.

Equations (2)--(4) prove, for the **actual** word count Z_n,

    log Z_n >=L h* - a_m log(L+1)-b_m/(mu L),         (5)
    h*=F(c*)=-sum_x pi_x log s_x.

No large graph enumeration or floating Stirling approximation is used.

## 2. Full finite-history divergence, not only pair convergence

Let Q_[0:L] be the ordinary stationary Markov law with initial pi and
transition P. View U_n on L+1 states by appending X_L=X_0. All its paths
have Q density pi_(x_0) exp(-L h*), and its root marginal is pi. Hence

    D_L=KL(U_n||Q_[0:L])=L h*-log Z_n+H(pi)<=A_L.    (6)

For 1<=k<=L let U_n^[k] be its k-transition block. Partition the L
transitions into floor(L/k) blocks plus a remainder. KL chain rule against
the Markov law and conditional convexity imply that each full block costs
at least KL(U_n^[k]||Q_[0:k]): forgetting the earlier past retains the
starting state and only decreases KL. Every block has the same marginal
by circular stationarity; the initial marginal KL is zero. Thus

    KL(U_n^[k]||Q_[0:k])<=A_L/floor(L/k),
    TV(U_n^[k],Q_[0:k])<=sqrt(A_L/[2 floor(L/k)]).    (7)

This pays whole histories. Matching single-step flows alone would not do
so. U_n itself is not asserted Markov at finite L, nor do we apply an
entropy-RATE estimate to its zero-entropy periodic extension.

For any nonempty preparation event D depending on X_0, both initial
marginals give it probability pi(D). Partitioning the KL by this event,

    TV(U_n^[k](.|D),Q_[0:k](.|D))
      <=sqrt(A_L/[2 pi(D) floor(L/k)]).               (8)

The factor for a rare preparation is retained. Coarse-graining to count,
parity or any bounded record contracts these same full-block bounds.

## 3. Irrational input and continuity are paid, not rounded away

For a fixed positive pi choose integer n(L) by flooring L pi and assigning
the remaining units to its largest fractional remainders, using a fixed
index rule for ties. Then pi^L=n(L)/L has |pi^L_x-pi_x|<=1/L. For sufficiently
large L it is positive. Apply (1)--(8) with its OWN optimizer P^L. These
are rational occupations, not a claim that Gibbs pi is rational.

The map pi->z is smooth on positive vectors. In logarithmic coordinates
the derivative of z(Kz) is diag(pi)(I+P), invertible: equivalently its
symmetric Hessian is the positive definite matrix paid by the earlier
coercive-dual proof. The ordinary finite-dimensional implicit-function
theorem and uniqueness give a smooth map near pi. Consequently on a fixed
carrier P^L->P, mu^L stays positive, and for fixed k (7) gives the actual
stationary P block limit. This is a fixed-carrier statement until the
explicit N-dependent payments in [the native bridge](NATIVE_LIMIT_ORDER.md).

This finite-ensemble representation retains a counting measure, support
and occupation constraint. It does not derive any of them from ontology.
