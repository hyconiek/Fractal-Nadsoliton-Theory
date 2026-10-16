# Tests of a kinetic-source bridge after MM

Date: 2026-10-16. Status: **ANALYTICAL_RESEARCH_SUBMISSION**, not an accepted
physical bridge or an execution of KB-001–003. The finite checker supplements
the arguments below; it is not a proof assistant. Statements concern the
declared model and classes, not the impossibility of a future FIN theory.

## 1. The existing occupation action cannot acquire a nonzero pre-GENERIC current by rewriting

**Inputs.** Reuse the admitted TH-002 interacting LOO path action, rather
than the independent-particle examples of Mielke–Peletier–Renger. Fix a finite
symmetric matrix A, finite gain g, and p in the interior of the q=12 simplex.
Let T={v: sum v_i=0}; covectors are understood modulo constants. Define

    s_i(p) = exp(g(Ap)_i) / sum_j exp(g(Ap)_j),
    F(p) = sum_i p_i log(12 p_i) - (g/2) p^T A p,
    a_i = log p_i - g(Ap)_i                   (modulo constants),
    H(p,xi) = sum_ij p_i s_j [exp(xi_j-xi_i)-1].

TH-002 already proves H(p,xi)=H(p,a-xi) and the corresponding
L(p,v)-L(p,-v)=a.v. These are inherited identities, not discoveries here.

**Candidate tested.** [Kraaij–Lazarescu–Maes–Peletier,
2017, §III.2, equations (69)–(70)](https://arxiv.org/html/1706.10115v1#S3.SS2)
derive pre-GENERIC from a symmetry recentered around a reversible current.
At the occupation-velocity level, retaining the FIN entropy S=-F and its
whole action would require some W(p) in T satisfying

    L(p,W+v) - L(p,W-v) = a.v                 for every v in T.       (1)

**Scoped result.** Such W is necessarily zero.

Proof. The Legendre dual of (1) is

    H(p,a/2-xi) - H(p,a/2+xi) = -2 xi.W.                           (2)

For completeness, substituting v -> 2W-v in the supremum defining H gives
H(p,a-eta)=H(p,eta)+(a-2eta).W; now put eta=a/2+xi. The inherited
TH-002 symmetry makes the left side of (2) identically zero. All tangent
covectors are allowed, so W=0. No numerical limiting argument is used.

One can see the same fact directly from the exact pair factorization

    H(p,a/2+xi) - H(p,a/2)
       = 2 sum_(i<j) sqrt(p_i s_j p_j s_i) [cosh(xi_j-xi_i)-1].     (3)

Every weight is positive. The expression is even in xi and strictly convex
on T. Its Hessian quadratic form is

    2 sum_(i<j) sqrt(p_i s_j p_j s_i)
           cosh(xi_j-xi_i) (u_j-u_i)^2 > 0  for nonzero u in T.

Thus the missing moving-frame current is not hidden in a second choice of
center for this same action. For any proposed W!=0, choosing xi=W makes
the failure of (2) have absolute size 2||W||^2 in these declared coordinates.
This is an exact mathematical falsifier, not an observable confidence test.

**Boundaries.** This does not forbid GENERIC with zero reversible component,
arbitrary algebraic drift decompositions that do not preserve the full action,
a different entropy, another time-reversal involution, an extended state
space, hidden divergence-free edge currents, momentum variables, or controlled
limits of a newly supplied microscopic theory. It concerns full occupation
velocity; occupation coordinates do not retain labelled identities. At the
simplex boundary some rates vanish; no boundary version is asserted.

The object H is a large-deviation Hamiltonian, not mechanical energy. In
particular, its costate xi must not silently become a physical momentum.
The result identifies the missing *nonzero independently sourced reversible
current/carrier* if this particular Hamiltonian-source route is pursued.
It does not select HB over MET/BAR, derive A7, or supply K07.

**Novelty.** The input action and ordinary reversal are old. The mathematical
argument is a direct application of published generalized-symmetry theory.
This submission exports its explicit FIN obstruction and falsifier; it makes
no claim of an original general theorem or a comprehensive priority search.

## 2. Conservation and Gibbs equilibrium alone still do not source kinetics

First reuse a basic finite-chain fact. For a finite irreducible reversible
generator Q with pi>0, if Q E=0 then

    0 = -<E,Q E>_pi
      = (1/2) sum_xy pi_x Q_xy (E_y-E_x)^2.

Connectivity forces E to be constant. Thus the native finite labelled
irreducible heat-bath chain has no nonconstant state function conserved for
every initial condition. A stationary expectation of an energy-like
function is a different, much weaker statement. On a reflected disconnected
domain, the conclusion is only componentwise. This elementary theorem is
reused, not claimed as a new general result.

Adding an energy ledger does not solve the source problem. Here is a fully
specified exact counterexample to the stronger proposed shortcut
"microcanonical total energy + Gibbs marginal + reversibility fixes rates".

Take system states x=0,1,2 with E=(2,1,0). At reservoir energy e=x there are
d=(1,2,4) distinct microstates b. The seven allowed pairs (x,b) all lie on
the declared total-energy shell E_x+e=2. A uniform shell law projects to

    pi=(1,2,4)/7 proportional to exp[-(log 2) E].

For each symmetric positive matrix a_xy (x!=y), put transition rate a_xy
between every pair of microstates belonging to distinct x,y, and no
within-x jumps. Use the negative row sum on the diagonal. This seven-state
generator G is symmetric, conserves the declared total energy, and has
uniform stationarity. Its projection is **strongly lumpable**, with

    Q_xy=d_y a_xy,  x!=y.

Indeed, G J=J Q for the indicator embedding J of macrostates; aggregate
exit rates do not depend on b. Hence this is equality of full projected
path laws, not merely projected mean drift. Detailed balance is exact:
pi_x Q_xy=d_x d_y a_xy/7=pi_y Q_yx.

Two choices give

    a^A_xy=1:       Q_A = [[-6, 2, 4], [1,-5,4], [1,2,-3]],
    a^B_01=2,
    other a^B=1:    Q_B = [[-8, 4, 4], [2,-6,4], [1,2,-3]].

They have identical energies, reservoir degeneracies and equilibrium.
They are not constant time rescalings: Q_B01/Q_A01=2 but Q_B02/Q_A02=1.

This is a countermodel in a *declared stochastic* class. It is not a
deterministic Hamiltonian source, an apparatus, an adopted FIN law or a
proof of thermalization. The reservoir degeneracies and transition clocks
were supplied. The fixed-shell energy is intentionally a bookkeeping
constraint; precisely this exposes why such a constraint alone cannot be
advertised as a new physical source. Multiple disjoint shells could be
included without removing the within-shell activity freedom.

The 2017 pre-GENERIC paper (§II.2) explicitly discusses adding an auxiliary
scalar to obtain full GENERIC. In the 2024 Hamiltonian derivation, energy
bookkeeping instead has an antecedent microscopic construction. These two
levels of justification must not be confused in a FIN proposal.

## 3. A quantitative regularity test for an autonomous compression candidate

The old repository Program 15 already supplies a Cauchy-measure unitary
dilation and a finite-dimensional recurrence obstruction. Do not reissue
either as a newly discovered source. [Mielke–Peletier–Zimmer, 2024,
Appendix C, Remark C.2](https://arxiv.org/html/2404.09284v2#A3)
also explicitly explains the generator-domain obstruction.

The following finite-window bound makes that regularity issue testable for
a proposed FIN completion. Let A=-Q be the nonnegative self-adjoint finite
native generator on L2(pi). For a real normalized centered observable f set

    C_FIN(t)=<f,exp(-t A)f>,
    a=<f,A f> > 0,       b=||A f||^2.

Suppose an asserted conservative representation uses an isometry I into a
Hilbert space, a strongly continuous unitary group U_t=exp(-it B), and
h=I f in Dom(B). This includes a stationary measure-preserving classical
flow through its Koopman representation when the observable has the stated
regularity. Write sigma^2=||B h||^2 and C_B(t)=Re<h,U_t h>.

For t>=0, the spectral theorem and the scalar inequalities
cos u>=1-u^2/2 and exp(-u)<=1-u+u^2/2 (u>=0) give

    C_B(t)-C_FIN(t) >= a t - (b+sigma^2)t^2/2.                    (4)

Consequently on any window [0,T], with
t_*=min(T,a/(b+sigma^2)),

    sup_(0<=t<=T)|C_B(t)-C_FIN(t)|
        >= a t_* - (b+sigma^2)t_*^2/2 > 0.                       (5)

If T>=a/(b+sigma^2), the lower bound is a^2/[2(b+sigma^2)].
Thus exact compression of this relaxing observable down to t=0 is
incompatible with this regularity. Approximation can remain possible as
the allowed second spectral moment increases, or after an explicit initial
time layer is excluded. This is a necessary test, not a sufficient source
or complete convergence theorem.

In the inherited native N2 heat-bath difference process every centered
mode decays at rate 2. Therefore a=2 and b=4. For the purely illustrative
bound sigma^2<=16 and T>=1/10, the correlation error is at least 1/10.
A two-frequency conservative control has C_B(t)=cos(4t); at t=1/10 its
defect from exp(-2t) is about 0.10233. The checker encloses this defect by
exact rational alternating-series bounds. No physical bandwidth, finite-N
confidence level or empirical time resolution is asserted.

**Crucial scope.** B is a spectral/Koopman generator; its second moment is
not automatically a mechanical energy or energy variance. This does not
rule out a positive-energy classical Hamiltonian bath, singular observables,
infinite baths, discrete/collision dynamics, or limits with a vanishing time
layer. The 2024 paper actually uses compression vectors outside the bath
generator's domain and a suitable enlarged space for positive-temperature
Gaussian initial data. It is not refuted by (4).

## 4. Composition-type check

[Pollard, 2017, Chapter 4](https://arxiv.org/html/1709.09743v2#Ch4)
builds open Markov processes by specifying boundary maps and gluing states;
its monoidal operation is disjoint union. The boundary black-box result is
about steady populations and currents. This is a useful correct composition
law, but it does not by itself give the joint state law of simultaneous FIN
cells: two q12 cells have 144 joint configurations, whereas the disjoint
union of two q12 state sets has 24 states. A further physical and probabilistic
typing map would be needed. Existing FIN incidence/role gates remain open.

## 5. Scientific disposition

Results 1–4 reject specific shortcuts and refine the contract for a possible
source. They do not pay the source itself. A noncircular internal conservative
completion would have to supply its carrier, evolution, coupling, preparation
and projection *before* solving backward from the desired FIN generator.
Full labelled cylinder laws or a declared bounded error are required where
block observations are claimed. Same equilibrium, a common covariance or
one projected semigroup matrix element is insufficient.

The outcome is **SCOPED_OBSTRUCTIONS_AND_SOURCE_CANDIDATE_SPECIFICATION**.
The model-to-physics bridge remains **OPEN_INPUT**; KB-002/003 are not activated.
