# Completing the stopped-law identification, and a bounded exit-side corollary

Date: 2026-10-14.
Status: **ANALYTIC_INTAKE_COMPLETION**, conditional on the admitted fold and
SF-001/002 estimates. This is not code-certified stochastic proof or peer review.

## 1. What needs making explicit

SF-003 correctly pays drift, bracket, shrinking jumps and fast-energy errors.
Testing only `C_c^3((-R,R))`, however, does not by itself specify absorption:
interior tests also fail to see reflecting or waiting behaviour at the boundary.
Its larger-radius argument is useful for transfer of inner exit events, but
first one must identify the absorption law at that larger radius too.

The following joint-clock argument supplies the omitted intermediate detail.
It does not change the microscopic generator, preparation or stop.

## 2. Fixed deterministic prepared start

Use the actual rounded type sequence from SF-001, fixed R,S and lambda, and
the joint soft/fast stopping time in **scaled time** `sigma_N<=S`. Set

    A_N(t)=t wedge sigma_N,
    Z_N(t)=Z_N(0)+integral_0^t f(Z_N(u)) dA_N(u)+M_N(t)+e_N(t),
    f(z)=-g_f(a lambda+b z^2/2).

The SF-002 remainder ledger gives `E sup|e_N|->0`. The stopped martingale has
bounded predictable bracket density, jump size `a_N=O(N^(-2/3))`, and

    sup_t |<M_N>_t-2g_f A_N(t)| -> 0

in probability (in fact a deterministic vanishing bound on the good local
domain). The bracket density is globally bounded for the finite-alphabet,
bounded-g family, including when a bad initial event is separately removed.

Second moments of M_N are bounded. Fourth moments are bounded uniformly as
well: for compensated jump martingales, the fourth-power formula has terms
`6 M_-^2 q2`, `4 M_- q3` and `q4`, where `q2<=C`, `|q3|<=Ca_N`,
`q4<=Ca_N^2`. Integrating, using `E M_t^2<=Ct` and Cauchy–Schwarz, bounds
`E M_S^4` by a fixed C(S). Doob's inequality also bounds the fourth moment
of the supremum. This pays uniform integrability both for M_N and for
`M_N^2-<M_N>`, not just formal convergence of their characteristics.

Aldous' increment bound for M_N follows from its bracket: conditional stopping
increments of length h have second moment at most Ch. Jumps vanish, so the
martingale is C-tight. The bounded drift integral and vanishing total-variation
error give C-tightness of Z_N. The clock A_N is bounded and 1-Lipschitz, and
sigma_N belongs to compact [0,S]. Thus jointly along a subsequence,

    (Z_N,M_N,A_N,sigma_N) => (Z,M,A,sigma),
    A(t)=t wedge sigma,

with continuous path components. Work in their joint **natural past
filtration**, completed and right-continuous; do not reveal the future value
sigma at time zero.

The uniformly integrable martingale and square-martingale identities pass to
the joint natural filtration, by [Whitt's Theorem 5.3 and its proof](https://arxiv.org/html/0712.1929v2#S5.Thmtheorem3).
Here the auxiliary process in that theorem includes Z_N and A_N. Theorem 2.1's
deterministic covariance clock is **not** invoked directly for the random
stopped bracket. Consequently

    Z(t)=z0+integral_0^t f(Z(u)) dA(u)+M(t),
    <M>_t=2g_f A(t).

The time sigma is a stopping time for the right-continuous joint filtration:
`{sigma<t}={A(t)<t}`, and right limits give the non-strict version. Brownian
representation, adding independent increments after sigma if necessary,
gives `M(t)=sqrt(2g_f) B(t wedge sigma)` in that filtration. This is a
representation of the limiting noise, not an extra microscopic bath.

## 3. Identify sigma with the first soft hit

Fast-first exit has probability tending to zero by SF-002. On its complement,
prelimit paths lie in `[-R-o(1),R+o(1)]`; if `sigma_N<S`, their terminal value
is a soft-boundary crossing with vanishing overshoot. Joint convergence and
continuous limit evaluation therefore imply

    |Z(t)|<=R for all t<=S,
    Z is frozen after sigma,
    sigma<S => |Z(sigma)|=R.

Let `tau_R` be the first time the corresponding unfrozen SDE hits ±R. The
drift can be extended boundedly and Lipschitz outside the interval. Pathwise
uniqueness identifies Z with this SDE before sigma. If `tau_R<sigma`, the
nondegenerate Brownian noise makes outward excursions immediately after
tau_R, contradicting that Z remains in [-R,R] before sigma. The immediate
crossing property follows, for example, from Brownian oscillations and the
bounded drift's O(h) increment near the boundary. Conversely, when sigma<S,
its boundary value forces `tau_R<=sigma`; an interior premature freeze is
impossible. Thus

    sigma=tau_R wedge S.

The law is the absorbed/frozen diffusion, not a reflected or sticky extension.
Uniqueness of the Lipschitz-extended SDE implies uniqueness of this stopped
joint law and identifies every subsequence. SF-003's process theorem follows.

For random prepared starts, retain convergence of the initial law, the stated
fast second moment and explicit exceptional-event probabilities. Localize away
the exceptional initial event first; do not claim a uniform numeric rate for
an unspecified o(1) event. The deterministic rounded class is the primary
admitted application and does not need that extension.

## 4. Exit variables without circularity

Now repeat the already completed stopped-law argument at one radius R'>R.
Before its larger stop, the limiting diffusion actually crosses the inner
boundary. The inner exit-time/side/freezing map is almost surely continuous
there. On a finite horizon the hitting time has no atom at its endpoint:
bounded-drift Girsanov equivalence reduces that zero-probability event to the
corresponding Brownian event. Thus capped inner exit time and the three outcomes
`UPPER_BY_S`, `LOWER_BY_S`, `SURVIVE_S` converge jointly. Boundary overshoots
vanish. Fast auxiliary stopping changes bounded-event probabilities by o(1).

This is stronger than finite-dimensional convergence and does not assume that
the hitting map of an already frozen path is automatically continuous.

## 5. Eventual exit side can be transferred, but not unbounded moments

The compact limiting diffusion exits almost surely. To see this without a
global cubic potential, extend its drift with a finite bound K. Over any
fixed block length h, a Brownian increment satisfying

    sqrt(2g_f) Delta B > 2R+Kh

forces a boundary exit from every initial point in [-R,R]. This event has a
strictly positive, uniform block probability. Independent Brownian future
increments give an exponentially vanishing limiting survival tail.

Let p_N^+(S) be the full parent chain's upper-first-exit probability by scaled
time S, without the auxiliary fast stop; let p_N^+(infinity) denote eventual
upper-first-exit probability. The microscopic chain is finite and irreducible
with positive refresh rates. For sufficiently large N both boundary sides are
accessible: the nonconstant `Xv` has maximum and minimum strictly around its
positive-p_f mean. Consequently this compact soft-boundary exit occurs almost
surely also at each such N.

For every fixed S, the completed theorem gives

    p_N^+(S)->p^+(S),   P_N(tau_N>S)->P(tau>S).

Since eventual-side indicators are bounded,

    |p_N^+(infinity)-p^+(infinity)|
       <= |p_N^+(S)-p^+(S)|+P_N(tau_N>S)+P(tau>S).

Take `limsup_N`, then `S->infinity`: the RHS vanishes. This adds a valid
**eventual exit-side convergence corollary**, with that explicit order of
limits. It requires the theorem for every fixed horizon, not just a single
stored finite-time array. No finite-N tail rate or mean-time convergence is
obtained. Rare finite-N long lifetimes could still affect unbounded moments.

For lambda=0, rounded p_f starts, exact `R=c=(4/b)^(1/3)` and canonical
boundaries ±1, the SF-004 rational interval therefore encloses the limiting
microscopic **local first-exit side** probability. It is not the probability
of exit by canonical time one, global phase switching or a physical event.
