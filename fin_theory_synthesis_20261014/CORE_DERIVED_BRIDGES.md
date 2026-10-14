# Derived structures: what can be obtained without adding a physical law

Date: 2026-10-14. Scope: the admitted dimensionless finite-copy model,
with its already declared leave-one-out heat bath. This law is an existing
modelling premise, not a source for physical time or nature's A7.

Notation: q=12, A=XX^T, X=X7, p=n/N, mu=X^T p,
s(p)=softmax(g A p). One attempt per copy is the supplied model clock.
Scalar rates below change if that clock changes. The proofs do not export
units, locality, fields of spacetime, quantum dynamics or gravity.

## 1. Interacting jump Hamiltonian — derived, not postulated

The exact count generator has rates

    n -> n-e_i+e_j at rate n_i softmax_j[g A(p-e_i/N)].

For a smooth test function f on the simplex define the nonlinear generator

    H_N f(p) = (1/N) exp[-N f(p)] L_N exp[N f](p).

For linear f(p)=xi·p this is exactly

    sum_(i,j) p_i softmax_j[g A(p-e_i/N)] [exp(xi_j-xi_i)-1].

At fixed finite g and bounded A its leading limit is

    H(p,xi) = sum_(i,j) p_i s_j(p) [exp(xi_j-xi_i)-1]
            = (sum_i p_i exp(-xi_i))(sum_j s_j(p) exp(xi_j))-1.

This follows directly from uniform continuity of softmax and the O(1/N)
self-exclusion. Terms with p_i=0 carry zero rate. A full path LDP additionally
requires a topology, initial-law LDP and comparison/tightness hypotheses;
TH-002 must check them rather than equating this generator calculation with
a completed global theorem. The relevant framework is
[Kraaij's finite-state interacting-jump theorem](https://arxiv.org/html/1502.05185v3).

The Hamiltonian here is a cumulant/rate-function object, **not** a conserved
physical energy or a quantum Hamiltonian. Its dual Lagrangian is a candidate
path cost derived from this stochastic process; it is not a new microscopic
variational postulate. No finite-time independent-Poisson-count assumption
has been inserted.

### Reversal identity

Let F(p)=D(p||u)-(g/2)p^T A p. For interior p, constants in grad F drop out
because H is invariant under xi->xi+c1. With f_i=log p_i-g(Ap)_i,

    p_i s_j exp(f_j-f_i) = p_j s_i.

Relabelling the indices proves

    H(p,xi)=H(p,grad F-xi),   H(p,grad F)=0.

These are exact identities for the limiting Hamiltonian. They do not pay
the global Gamma=B4/capacity/prefactor gate or establish a unique global
quasipotential without the appropriate basin/reachability theorem.

## 2. Closed mean does not imply closed trajectory fluctuations

Differentiating H at xi=0 gives

    b(p)=s(p)-p,
    D(p)=diag(p)+diag(s)-p s^T-s p^T.

D is the leading count quadratic-variation matrix; the covariance rate of
p_N is D/N+O(N^-2). These leading formulas are already present in report131
and the process convention dictionary138; they are reused, not claimed as
new discoveries.

The retained mean drift closes:

    dot mu = X^T softmax(g X mu)-mu.

But the tilted retained Hamiltonian is

    H_mu(p,eta)=H(p,X eta),

which contains sum_i p_i exp(-X_i·eta). This is not generally determined
by mu. Even its quadratic term can distinguish the hidden fibre.

### Exact finite-N witness, no simulation campaign

Set epsilon=1/60 and h_j=cos(2*pi*2j/12). Define

    p_+=u+epsilon h,   p_-=u-epsilon h.

At N=120 these are actual count states: n_j=10 +/- 2h_j, where
2h=(2,1,-1,-2,-1,1,2,1,-1,-2,-1,1). Both have exactly mu=0 by Fourier
orthogonality. Use the already declared g=0 control, so refresh is uniform.

For the retained cosine-k5 coordinate,

    x_(j,c5)^2 = (Lambda5/12)(1+h_j).

Since sum h_j=0 and sum h_j^2=6, the generators applied to mu_c5^2 at these
states differ by exactly

    epsilon Lambda5/N = Lambda5/(60*120) > 0.

Thus the mean observation is not an exact autonomous Markov state for all
preparations, even though its deterministic drift closes. A mean-only local
path Hamiltonian would already fail to reproduce this second derivative.
This is a scoped constructive consequence, consistent with earlier hidden-
noise/rank-four results; it is not a no-go for every limiting reduction.

### Observable lift, not new material degrees of freedom

Let S(p)=X^T diag(p)X. Given mu and the leading retained quadratic variation,
the known s(mu) terms recover S. Unscaled Fourier products satisfy

    c3*c4=(c1+c5)/2,    s3*c4=(-s1-s5)/2,
    c3*c5=(c2+c4)/2,    s3*c5=(-s2-s4)/2.

Thus first and second feature moments recover all four missing k1/k2
components; together with normalization they span all12 label functions.
The lifted state determines the full occupation vector p. It is not a magic
seven-variable stochastic compression or recovery of every labelled copy's
identity. These are functions of existing labels, not inserted particles,
physical geometry or a new information substrate.

## 3. Actual collective drift versus Gaussian auxiliary field

For constant g>0 define theta_ret=g mu and

    Phi_g(theta)=||theta||^2/(2g)-log(mean_j exp(X_j·theta)).

The existing leading microscopic drift gives exactly

    dot theta_ret = -g grad Phi_g(theta_ret).

This identifies the **model** scalar mobility in the retained deterministic
coordinates; it does not derive the heat-bath premise itself. Leading noise is

    D_theta(p)=g^2 X^T D(p) X,

not a freely chosen isotropic covariance. At an equilibrium p=s,
D_theta=2g^2 M(p), M=X^T(diag p-pp^T)X.

The Hubbard-Stratonovich variable theta_aux is different. The existing
Gaussian-parent identity gives theta_aux conditional on labels as a Gaussian
with mean g mu and covariance(g/N)I. Consequently its stationary marginal
is a Gaussian convolution of the discrete theta_ret distribution. Matching
that marginal to exp(-N Phi) does not identify the microscopic trajectories.
Stable-phase covariance compatibility is already covered by MP7-033.

## 4. Intrinsic stochastic fold: promising target, not yet a theorem

At a certified simple fold (theta_f,g_f), let v be the unit null vector of
Hess Phi. Then M_f v=v/g_f, hence the leading retained noise along v is

    v^T D_theta(p_f)v = 2g_f.

No new noise amplitude is needed in that channel. Other directions generally
do not have the isotropic MP7 auxiliary-Langevin covariance.

Let a=v·partial_g grad Phi and b=D^3 Phi[v,v,v], with the orientation fixed
by the existing fold certificate. The balance of the actual drift/noise
suggests the scaling

    g_N=g_f+lambda N^(-2/3),
    theta_ret-theta_f along v = N^(-1/3) Z,
    t=N^(1/3) s,

and the candidate stopped limiting law

    dZ = -g_f [a lambda+(b/2)Z^2] ds + sqrt(2g_f) dB_s.

TH-004/005 must derive coordinate corrections, verify coefficients, control
the fast modes/initial layer and prove convergence up to a declared exit.
This cubic local model is not a normalizable global equilibrium potential;
its stopped/exit scope matters. It is not a proved FIN global mixing law.
The analogous Potts result is a useful benchmark, not a transferred theorem:
[Cuff et al.](https://arxiv.org/abs/1204.4503) obtain a spinodal window of
order n^(-2/3) and n^(4/3) mixing in single-update steps. FIN uses N total
attempts per unit time, so clocks must be converted before comparing powers.

This candidate is intrinsic mesoscopic statistical physics of the supplied
model. It supplies neither propagating spatial fields nor a physical second.
Its Brownian driver has Gaussian increments; “non-Gaussian collective law”
refers to the nonlinear state/exit distributions, not non-Gaussian white noise.

## 5. Literal existing relational costs: two different obstructions

### Change-only cost: identity freezing

The cost displayed in GEN-676/702 is of the form

    C(R,R';M,M') = alpha ||R'-R||^2
                 + beta D(P_R||P_R') + gamma ||M'-M||^2.

If weights are nonnegative, D is KL between defined probability laws and
the identity successor is admissible, every term is nonnegative and the
identity has cost0. With alpha,gamma>0, minimization forces R'=R,M'=M.
Even without uniqueness, the identity is a minimizer. Missing P_R definitions
do not justify asserting a nontrivial evolution. An imposed compulsory memory
update or exclusion of identity would change the premises and must be sourced.

Here R denotes the relation variables actually included in the squared
penalty. If a version of the proposal also contains unpenalized x_i variables,
the proof freezes the penalized relations/memory, not every unspecified
variable. Zero-cost freedom in such variables does not itself select a
generative law; it is another missing transition/selector definition.

This exact conditional obstruction shows why penalizing change alone cannot
derive generative dynamics. It does not refute other already-defined costs.

### State cost plus inertia penalty: separability

GEN-931/932 instead displays a scalar-entry structure cost plus squared
change penalty. On a product edge domain,

    C(R,R_t)=alpha sum_e f(r_e)+beta sum_e(r_e-r_(t,e))^2

has a product argmin set. With unique scalar minimizers the update of each
edge depends only on itself: zero cross-edge influence. Ties, nonproduct
constraints or other undefined C_I/C_M terms must not be silently selected
to create interactions. For alpha,beta>0 and feasible R_t,

    F(R_(t+1)) <= F(R_t)-(beta/alpha)||R_(t+1)-R_t||^2.

This is a Lyapunov inequality, not conserved physical energy. Nontrivial
cycles and changes consisting only of a constant-F permutation cannot occur
under that exact minimization. It does not certify any claimed numerical
SOC/soliton result or exclude all nonseparable relational models.

These obstructions use existing formulas; no triangle reward, graph locality,
controller, new force, f choice or scale has been appended to repair them.
