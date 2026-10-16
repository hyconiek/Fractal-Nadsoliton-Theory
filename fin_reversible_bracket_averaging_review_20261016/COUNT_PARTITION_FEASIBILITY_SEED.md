# A concrete next object: intrinsic count partition and an exact flow gate

Date:2026-10-16. **Analytical candidate preparation, NOT an admissible
counterfamily or completed IBA refutation.** No physical law/state is added.

## 1. Freeze one partition using only existing counts

Use lambda0, exact g_N=g_c, N in12N and the actual MT D_N. Let

    Z_N=N^(1/4)m,
    C(n)=sum_j n_j cos(pi j/2), S(n)=sum_j n_j sin(pi j/2),
    R_N=2[C(n)^2+S(n)^2]/N.

C and S are integers: the harmonic entries are0,+1,-1. R_N is the squared
radius of the existing normalized cosine/sine k3 fast coordinates; it is
neither a spatial radius nor the earlier BA Poisson cost. Put
gamma3=1-Lambda3/Lambda6>0 and freeze

    epsilon=1/8, p(z)=1/2+epsilon exp(-z^2),
    T(z)=-(2/gamma3) log(1-p(z)),
    A_N={n in D_N:R_N(n)<=T(Z_N(n))}, B_N=D_N\A_N.

The symbols A_N,B_N here denote a partition, not a changed interaction
matrix or a new preparation. A future candidate generator would use only
existing count states. Its phase label is a function of those states,
not an appended frozen hidden bit. The threshold is fixed before any
feasibility result; no fitting or threshold search is proposed.
This chosen analysis partition is not a uniquely derived fundamental FIN law.

## 2. Static properties already follow from the admitted model

MT's joint local equilibrium gives Z_N->Z with
rho0(z) proportional to exp(-z^4/12), independently of the fast Gaussian
pair. The normalized k3 features have variance1/gamma3 each: their
Euclidean normalization gives the factor12 in converting the MT precision
12I-g_c A. Thus R_N converges to an independent radius with CDF

    P(R<=t)=1-exp(-gamma3 t/2), t>=0.

The boundary R=T(Z) has zero probability. Therefore for bounded continuous
f, and also f(z)=z^4 by the actual MT moment UI,

    E_piD[f(Z_N)1_A_N] -> E_rho0[f(Z)p(Z)].

Write mu=E exp(-Z^2), M=E[Z^2 exp(-Z^2)]>0. The phase masses tend to

    a=1/2+epsilon mu, b=1/2-epsilon mu,

both bounded away from zero. Nonempty positive-mass phases are supplied
asymptotically; no small-N MT-domain claim is made.

The partition passes one necessary **exact** finite-N balance test. The
circulant Gibbs energy, outer cylinder, sublevel and balanced seed are D12
invariant, hence the actual connected component D_N is invariant. A one-step
label rotation flips m and rotates the k3 pair; R_N and Z_N^2 are invariant.
Thus both A_N and B_N are rotation invariant. Native reflected HB covariance
gives d(Rn)=-d(n). Consequently

    sum_(n in A_N) pi_D(n)d(n)=sum_(n in B_N) pi_D(n)d(n)=0.

This necessary balance is NOT sufficient to realize closed diagnostic
phases; finer internal cuts and the rate cap can obstruct them.

## 3. What would follow IF admissible closed-phase rates exist

Suppose an actual sequence of reversible generators in the unchanged KR/BA
class has no edges between A_N and B_N, retains all statewise native parity
drifts and one uniform CN cap. Internal reducibility is allowed. Each phase
is stationary under the corresponding conditional Gibbs law, so

    E[G_N|A_N]=-2N E[m d|A_N],
    E[G_N|B_N]=-2N E[m d|B_N].

The already paid weighted native drift approximation gives N m d=Z b0(Z)
in stationary L1, where b0(z)=-z^3/3. Hence, with r_N=G_N-EG_N,

    E[r_N|A_N] ->-4epsilon M/a,
    E[r_N|B_N] ->+4epsilon M/b.

To see the signs and coefficients without quadrature, integration of
(z exp(-z^2)rho0(z))' yields

    E[(Z^4/3-1)exp(-Z^2)]=-2E[Z^2 exp(-Z^2)]=-2M.

Closed-phase invariance and Jensen at time S now give

    liminf E sup_(t<=S)|integral_0^t r_N(X_sqrt(N)s)ds|
       >=8S epsilon M=S M>0.

Thus a feasible closed-partition family would refute universal IBA without
having to first identify a whole different limiting SDE. This is an
IF statement. No such family has been constructed in this intake. The
global parity-cut/conditional-mean theorem remains valid and is compatible
with these different phase means because phase membership depends on Z.
The obstruction uses the full unpostselected stationary pi_D start;
conditional phase means are proof tools, not an altered preparation.

Any such family must change the stationary bracket on nonvanishing mass:
the same conditional-mean calculation and native L1 concentration imply

    liminf E_piD|G_tilde-G_native|>=8epsilon M=M>0.

Thus deleting crossing edges and repairing rates only in a vanishing-mass
interface band cannot work under a uniform rate cap. A local wheel/spike
or a boundary-only rejection repair is not the desired construction.

## 4. Conditional extension to strictly positive irreducible rates

If a closed-partition candidate exists, set eta_N=N^-2 and

    Lhat_N=(1-eta_N)Ltilde_N+eta_N L0_N.

It has the same Gibbs law and exact drift, the same single-copy support
and a uniform max(C,1)N cap. Every native permitted edge is now positive;
the connected actual D_N makes Lhat irreducible. This is a diagnostic
construction, not a new primary law or a fitted time rescaling.

The original candidate and this mixture can be coupled by common
(1-eta_N)Ltilde jumps. Different extra jumps have total rate at most
eta_N(C+1)N. Over model time sqrt(N)S, full count-path TV is at most
(C+1)S eta_N N^(3/2)=O(N^-1/2). Their centered G functions differ uniformly
by O(eta_N), because both means coincide exactly from matched drift.
The additive bracket functionals are uniformly bounded on fixed S, so
the strictly positive family inherits the same nonzero IBA obstruction.

This IF implication does not construct the missing closed family. It
shows why irreducibility alone would not rescue the universal claim;
an N-independent lower-rate comparison is a different hypothesis.

## 5. The only construction gate: exact finite flow feasibility

For every permitted undirected count edge e={x,y} with both endpoints in
the same phase, introduce c_e>=0. Let

    (Ac)_x=sum_y c_xy(m_y-m_x), b_x=pi_D(x)d(x),
    (Bc)_x=sum_y c_xy, v_x=C N pi_D(x).

Seek **Ac=b, Bc<=v, c>=0**, with one finite C independent of N. This pays
the actual drift, detailed balance, support, cap and reflection. No edge
may cross the partition or leave D_N. Same-parity edges have zero A column;
they can be omitted without loss for nonnegative-rate feasibility. This
does not assert connectivity or positivity on every permitted edge.
All three of c,b,v can consistently use unnormalized Gibbs weights;
the common positive partition-function factor cancels.

The following exact dual test is necessary AND sufficient. For all real
phi_x and beta_x>=0 satisfying on each permitted within-phase edge

    (m_y-m_x)(phi_x-phi_y)+beta_x+beta_y>=0,

one must have

    sum_x phi_x b_x+sum_x beta_x v_x>=0.             (F)

Necessity follows by multiplying the primal equations/inequalities.
For sufficiency add nonnegative slack w with Bc+w=v. The target (b,v)
must lie in the finitely generated cone with columns (A_e,B_e) and (0,e_x).
Such a cone is closed: reduce conic representations to linearly independent
supports; a convergent sequence has a subsequence in one of finitely many
supports, where the inverse coefficient map is continuous. If the target
is outside, project it onto the closed cone. The vector from the target to
its projection has nonnegative dot product with every column and negative
dot product with the target. Its last coordinates are beta>=0, producing
a violated (F). This proves the alternative without an assumed flow solution.

The exact small dual fixtures only check this algebra. They do not pay
(F) on the huge actual all-N count graph. Phase total balance alone, a
relaxed drift, fitted capacities or a small unrelated LP is not enough.

## 6. Interpret only the result actually obtained next

A constructed uniformly capped feasible sequence gives the explicit IBA
obstruction above. An infeasibility certificate rejects THIS partition and
cap/sequence, not all possible counterfamilies and not universal IBA.
No feasible/infeasible all-N result is claimed now. Any asymptotic sequence
argument must retain the exact microscopic drift, not just b0(Z).

This seed changes the next research question from another generic Poisson
bound to a concrete native count-graph feasibility theorem. It supplies no
physical activity/A7 source, spatial graph, new substrate, experiment or
quantum/gravity/ToE claim.
