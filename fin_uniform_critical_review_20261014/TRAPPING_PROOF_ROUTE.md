# A proposed finite-N local-equilibrium/trapping route

Date: 2026-10-14.
Status: **PROOF_STRATEGY**, not an executed metastability theorem.
The equilibrium exclusion in GIBBS_EXCLUSION_PROOF is established separately;
none of the dynamic time-window targets below follows from that bound alone.

## 1. Local free-energy coercivity from the same model

Use `p=p_m+w`, `p_m=u+m sigma/12`, `w` orthogonal to 1 and sigma. Exactly

    F_g(p_m)=((1+m)log(1+m)+(1-m)log(1-m))/2
             -(g Lambda6/24)m^2.

At g_c the quadratic term cancels and the remaining even series starts
`m^4/12+m^6/30+...`. The fast first derivative at p_m is zero by shift-two
symmetry. The positive fast Hessian at (u,g_c) remains uniformly positive in
a sufficiently small **fixed** neighbourhood. Taylor integration in w then
targets, for bounded lambda and all large N,

    F_gN(p)>=c_m m^4+c_w||w||^2-C Lambda^2/N,

where the displayed capital Lambda bounds |lambda|, not a Fourier eigenvalue.
For example the scalar Young bound preserves m^4/24 and pays
`(3/2)(chi*(g_N-g_c))^2`. Prove the full domain/constant statement; do not
borrow a merely pointwise Hessian minimum.

## 2. Fix the stopping domain and graph, not a physical wall

One possible construction is a small cylinder

    U_epsilon={|m|<epsilon, ||w||<epsilon^2},

inside the simplex interior. Choose h>0 below its entire boundary's lower
energy, uniformly for all large N. Define D_N as the communicating component
containing a balanced-uniform type of

    {p=n/N in U_epsilon : F_gN(p)<=h}.

This component is a mathematical observation domain. No jump is suppressed
in the **actual** process. A reflected comparison will be introduced only if
and when needed for an explicit coupling proof.

Prove nonemptiness, iid-start coverage and boundary properties. In particular,
a sufficiently small fixed central type ball should connect to the balanced
type by genuine count paths of low energy; it is not enough that a continuous
set looks convex. Do not enumerate astronomical count lattices.

For the inner exit layer

    B_N={n in D_N : q_N(n,D_N^c)>0},

bounded jumps and a paid boundary Lipschitz bound should imply
`F_gN(n/N)>=h-C/N`, hence eventually a fixed positive lower bound b<h.
The artificial outer cylinder must not create an unpaid cheaper exit.

## 3. A useful initial-law domination identity

Let mu0_N be the **actual iid-uniform initial type distribution**, and let

    W_N(n)=mu0_N(n) exp[(Ng_N/2)p^T A p]

be the unnormalized Gibbs weight. Since A>=0 and eventually g_N>0,
`mu0_N(n)<=W_N(n)` for every type. W_N is invariant for each declared
finite-N reversible generator. Positivity of its Markov semigroup therefore
gives the exact measure domination

    mu0_N P_t <= W_N P_t = W_N.

This is a comparison with an **unnormalized** invariant measure; it is not
the false assertion that iid-uniform is the Gibbs initial state.

If all rates out of one count state sum to at most C_rate N, expected outward
jump count gives the candidate bound

    P_mu0(tau_D<=T)
       <=mu0_N(D_N^c)+C_rate N T sum_(n in B_N) W_N(n)
       <=mu0_N(D_N^c)+C_rate N T (N+1)^12 exp(-N b).

Derive the count compensator and domain assumptions before admitting it.
This does not need normalization by an exponentially large global partition
function; that factor cancels in the invariant-measure argument.

For each fixed c<b, T_N=exp(cN) would then yield a vanishing exit probability
once iid coverage is paid. This is a **one-sided trapping result**, not
`tau~exp(Nb)`, a true communication height or a global capacity/prefactor.
For HB/MET/BAR the maximum rate constants follow from their fixed microscopic
rules, not physical time calibration.

## 4. Match a conditional local equilibrium, not the full Gibbs law

The auxiliary equilibrium is

    pi_D,N(n)=W_N(n)1_(n in D_N) / sum_D W_N.

At large N, the target joint fluctuation law is quartic in `Z=N^(1/4)m` and
Gaussian in `Y=sqrt(N)E_f^T w`, with the actual fast Hessian. Symmetry kills
the dangerous linear fast energy term; higher mixed terms and count-lattice
Riemann/Stirling errors still need proof. Show tightness and required moment
uniform integrability, not just pointwise Taylor expansion. A local partition
prefactor is optional unless a subsequent bound actually uses it.

This should match UF's limiting diffusion equilibrium. It remains conditional
local Gibbs, while the admitted global Gibbs exclusion keeps the full parent
late equilibrium distinct.

## 5. The actual missing relaxation gate

Suppress only jumps out of D_N in an auxiliary reflected generator; establish
its graph connectivity and invariant pi_D,N. Couple it with the actual chain
using unchanged clocks until first exit. A reflected gap/mixing bound is new
work: finite-horizon UF weak convergence does not supply it.

A useful target is any proved polynomial relaxation time, not necessarily
the optimal critical prefactor. If `gap_D,N>=c0 N^(-a)` and the initial density
relative to pi_D,N is polynomially bounded, an explicit L2/TV argument yields
a polynomial t_loc(N). The latter initial-density bound can use local
coercivity/type weights, but must be paid for the iid start, not asserted.

Conversely, a Rayleigh test with m and the conditional quartic variance should
give an upper gap bound of order N^(-1/2). A claimed faster uniform gap must
respect this obstruction. No exact Markov parity quotient may be assumed.

With a **proved** relaxation gate and the above exit bound, one can target

    sup_(t_loc(N)<=t<=exp(cN))
      TV(Law_mu0(p_N(t)),pi_D,N)->0.

State each law, conditioning and coupling error precisely. If that gate
cannot be paid, export LOCAL_RELAXATION_OPEN and stop that promotion; do not
fit a decay pole or insert a reflecting wall into the physical parent.

## 6. Interpretation

The potential payoff is an actual finite-N effective statistical regime with
separated local relaxation and sector escape times. It is not a proof of an
autonomous small-N memory bit: sign-information decay and sector lifetime
are different quantities. It does not source spatial propagation, physical
energy/time, natural initialization/gain or a fundamental theory.
