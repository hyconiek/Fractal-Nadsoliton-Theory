# AD-004: actual original-source critical weak path theorem

Date: 2026-10-20. Self-audited conditional analytical theorem, not formal
proof, independent review, physical source or primary FIN replacement.

## 1. Statement and exact preparation

Let Y_j evolve the reflected ORIGINAL entropy-selected count kernel P_D,N,
selected from full A7 Gibbs before reflection. Start Y_0 exactly pi_D,N.
For fixed finite S and bounded fixed lambda, N through multiples of 12,
put h_N=1/(d_N sqrt(N)) and

    Z_N(s)=N^(1/4)m(Y_floor(s/h_N)), 0<=s<=S.

Then in D([0,S]),J1,

    Z_N => Z,
    dZ=q(chi lambda Z-Z^3/3)ds+sqrt(2q)dB,
    Z0~rho_lambda, q=12,                            (1)

with future Brownian increments independent of Z0 in the joint natural
past filtration. This is not a statement about the full global Gibbs
stationary preparation, which has different mass, or an arbitrary start.

## 2. Actual discrete martingale, containment and UI

On the critical grid the exact increment martingale is

    Delta M_j=Delta Z_j-h_N beta_N(Y_j),
    E[(Delta M_j)^2|F_j]=h_N[G_N-h_N beta_N^2](Y_j).

Its jumps are <=4N^-3/4. Raw G_N<=4d_N/N<=4q gives a deterministic
bracket bound on every fixed horizon, up to one grid step. A fourth-power
discrete expansion, conditional mean zero and
E[(Delta M)^4|F]<=16N^-3/2 E[(Delta M)^2|F] give a uniform fourth moment
of the martingale maximum. The terms 6M^2 E(Delta M)^2 and
4M E(Delta M)^3 are bounded using E M_j^2<=C j h_N and the jump bound,
as in SQ; this pays square-martingale UI.

For the soft radius stop retain the nonnegative terminal square in

    E[Delta(Z^2)|F_j]=h_N[2Z beta_N+G_N].

Since 2Z f(Z) is bounded above, AD-003's weighted drift payment and
stationarity bound the stopped expected square by E Z0^2+C S+o(1).
Consequently the probability of hitting radius R is <=C_S/R^2, including
initial exceptions and overshoots. This is maximum containment, not an
unsupported bound at random times from one-time stationarity.

On compact stops f is bounded, the integrated residual tends to zero in
L1 and martingale increments have bounded bracket. The drift Riemann sum
and piecewise-constant integral differ by at most
h_N max|beta_N|<=2N^-3/4. These estimates and vanishing jumps yield
C-tightness. AD-003(5),(6) identify the predictable bracket as 2q s;
its grid endpoint error is bounded by 4q h_N.

Pass both martingale identities against bounded continuous cylinders of
the initial state and past (Z_N,M_N). Fourth-moment UI and continuous
limits justify passage to expectations. The limit M is continuous with
M^2-2q s a martingale in the JOINT natural past filtration, so it is
sqrt(2q) times Brownian motion in that filtration, with increments
independent of Z0. The same locally Lipschitz cubic SDE as SQ is uniquely
nonexplosive: its generator on z^2 is bounded above by
2q+(3q/2)(max(chi lambda,0))^2. Therefore the full sequence converges.
Its zero-current quartic stationary law is the actual MT initial limit.

No full count-path TV comparison with an auxiliary CTMC was assumed at
any step. A bounded-jump discrete process, not inserted continuous noise,
supplies the limiting Brownian martingale.

## 3. Actual unreflected late iid parent

The IS theorem supplies the local window for THIS original source. Let
k_N be a deterministic intrinsic-step origin with k_N>=k_loc and
k_N+floor(d_N sqrt(N)S)<=floor(d_N exp(zeta N)), fixed zeta<b.
Its unreflected iid count law at k_N is TV-close to pi_D by e_N^prep->0.

Maximally couple those initial count laws and propagate the same full P_N.
The complete future-path TV costs e_N^prep once. Couple the pi_D start
to its stationary reflected comparison until a genuine departure. Over
m_N=floor(d_N sqrt(N)S) steps the additional error is at most
m_N pi_D(B_D), which vanishes exponentially by MT. Thus the actual
unpostselected parent obeys

    N^(1/4)m(Y_(k_N+floor(d_N sqrt(N)s))) => Z(s)      (2)

with limit (1). Initial iid is a stipulated statistical diagnostic of the
given hypothesis kernel, not selected by the stationary entropy axiom.
No unbounded iid moment/entropy transfer is inferred from this TV argument.

The leading critical amplitude law matches the earlier SQRT process in
the already frozen normalization. Finite-N branching/activity remains
different from HB/SQRT, and richer histories need not be statistically
equivalent. This does not identify nature's kinetic law, physical clock,
quantum noise or a fundamental source of A7/Gibbs/entropy maximization.
