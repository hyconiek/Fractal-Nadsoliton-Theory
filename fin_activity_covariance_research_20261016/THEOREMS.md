# Derived finite activity/noise objects and the critical-law handoff lemma

Date:2026-10-16. Status: **ANALYTICAL_MODEL_THEOREMS_WITH_BOUNDED_CHECKS**.
Not formal proof-assistant verification or peer review. No new microscopic
law is adopted; all comparison rates below are diagnostic constructions.

## 1. Exact inverse: predictable occupation covariance determines activity

Fix N>=1, q>=2, a known finite count domain D, and strictly positive
unnormalized stationary weights W(n). The only nonnull jumps are single
replacements n->n-e_i+e_j within D. Assume detailed balance with W and use

    Gamma(f,g)=L(fg)-f Lg-g Lf, p_i=n_i/N.

For i!=j no other single replacement changes both coordinates. Consequently

    -N^2 Gamma(p_i,p_j)(n)=Q(n,n-e_i+e_j)+Q(n,n+e_i-e_j),

with absent edges contributing zero. Fix all counts outside {i,j}. The
remaining states in D split into consecutive path segments. Order a
segment n(0),...,n(M) by increasing n_i and write

    S_k=-N^2 W(n(k)) Gamma(p_i,p_j)(n(k)),
    C_k=W(n(k)) Q(n(k),n(k+1)), C_-1=C_M=0.

Detailed balance gives S_k=C_(k-1)+C_k. Therefore

    C_k=sum_(ell=0)^k (-1)^(k-ell) S_ell, 0<=k<M.

Divide by W at the appropriate endpoint to recover both directed rates.
Every single-replacement edge belongs to exactly one such fibre, proving
uniqueness of the whole off-diagonal generator. Diagonals are negative row
sums. No partition function, irreducibility or global connectedness is needed.

### Compatibility, not just a formal inverse

A proposed exact tensor comes from a reversible single-replacement chain
on the known D if and only if:

1. It is symmetric, has zero row sums and nonpositive off-diagonal entries.
2. Every reconstructed fibre conductance C_k is nonnegative.
3. At each terminal endpoint S_M=C_(M-1), including S_0=0 for a singleton.

Necessity follows from the jump formula. For sufficiency construct rates by
the recursion. Their conductances are symmetric and nonnegative. They give
every off-diagonal tensor entry; zero row sums then give every diagonal
entry. The resulting tensor is automatically positive semidefinite.
If all admissible edges must be present, require strict positivity; known
forbidden edges must reconstruct as zero. Unknown exits or a different
reflection rule invalidate the contract, not just one numerical value.

Nonnegative node sums are insufficient. S=(1,1,1) fails the terminal check.
S=(2,1,1,2) passes that check but reconstructs C=(2,-1,2), which is invalid.
The code retains and rejects both fixtures, never projects them to success.

For FIN use the **count** weights

    W_N(n)=N!/(product n_i!) exp[g n^T A7 n/(2N)].

Omitting multiplicity generally changes the inverse. An exchangeable
labelled realization can use rates Q(n,n-e_i+e_j)/n_i for an occupied
source slot. This does not identify arbitrary hidden or labelled mechanisms.

### Stability and a minimal datum

For two compatible datasets on an M-edge fibre with max|delta S|<=epsilon,
recursion from both endpoints gives

    |delta C_k|<=min(k+1,M-k) epsilon.

Conversion to rates divides by W(n); no uniform rate-error guarantee follows
at tiny weights. A common unknown multiplier of the whole tensor changes
the reconstructed generator by that multiplier; normalization of W cancels.

At a genuine endpoint n_i=0, n_j>0,
Q(n,n+e_i-e_j)=-N^2 Gamma(p_i,p_j)(n). At n=N e0 the ratio
Gamma(p0,p3)/Gamma(p0,p6) is the outgoing count-rate branching ratio.
It changes by h3(g) under the mandatory distance-three thinning. Thus an
independently derived restriction on that noise ratio would be a genuine
activity restriction. It is NOT supplied by computing Gamma from HB.
The fully aligned endpoint is outside the small MT critical domain.

The support assumption is essential. On q2,N3, permit the extra jump0<->3
in the count coordinate k. Conductance changes on edges01,12,23,03 equal
epsilon,-epsilon,epsilon,-epsilon/9 preserve every node's second jump
moment while changing rates. Small epsilon keeps positivity. The path-fibre
inverse must not be applied to such a multi-count-jump model.

This is a weighted path-incidence inversion specialized to FIN, not a claim
of an original general inverse theorem. Targeted duplicate search found
the old forward covariance/Onsager identities, not this explicit recursion.

## 2. New obstruction: even the full pointwise drift leaves activities free

Let q>=3,N>=3 and choose distinct labels a,b,c with n_a,n_b,n_c>=1.
Set the cyclic roots

    r0=e_a-e_b, r1=e_a-e_c, r2=e_b-e_c,
    r3=-r0, r4=-r1, r5=-r2.

They sum to zero and r_(j-1)+r_(j+1)=r_j. The six neighbors n+r_j are
legal; each consecutive pair differs by a single replacement.

Let S have symmetric value+1 on the six centre-to-neighbor spokes and-1
on the six ring edges. Then at the centre sum_y S(n,y)(y-n)=0; at ring
vertex n+r_j the sum is r_j-r_(j-1)-r_(j+1)=0. It vanishes elsewhere.

For any native reversible generator Q0 with positive rates on these edges,
put c0_xy=W(x)Q0(x,y) and

    Qs(x,y)=Q0(x,y)+s S(x,y)/W(x).

For sufficiently small positive s, conductances remain positive and
reversible. For EVERY occupation coordinate,

    (Ls-L0)p_i=0 exactly at EVERY state.

But at the centre

    (Gammas-Gamma0)(p_i,p_j)
       =s/(N^2 W(n)) sum_ell (r_ell)_i(r_ell)_j.

The nonzero three-label matrix is 6I3-2*11^T. In particular its aa entry
is4s/(N^2W(n))>0. Centre exit rate increases by6s/W(n), so even its holding
law changes. Equal pointwise coordinate drift does NOT imply equal finite
time mean trajectories, because the common nonlinear drift is averaged
under different state laws. No stronger equality is claimed.

### Actual FIN support, symmetries, proposal clocks and g0

Use the normalized average S_bar=(1/24)sum_(h in D12) h.S over D12 label
rotations/reflections. Each term still annihilates
affine functions. For q12,N3, spoke endpoints have support sizes3 and2;
rim endpoints both have support size2, so opposite signs cannot cancel.
The averaged stress has60 affected states and132 undirected edges,72 positive.
Exchangeable labelled rates divide count rates by departure multiplicity.
The labels remain analysis labels, not spatial cells.

Native HB has positive null-refresh margin at every occupied source label.
For the averaged S let B_(x,i)=sum_j max(S(x,x-e_i+e_j),0). Define

    d(g)=min{c0_xy/(2|S_xy|) on negative edges,
             W(x)[n_i-sum_j Q0(x,x-e_i+e_j)]/(2B_(x,i)) when B>0},
    s(g)=g^2/(1+g^2) d(g).

Every entry in the finite minimum is positive. Negative changes consume
at most half an original conductance; positive sourcewise changes consume
at most half the null-refresh margin. Thus positivity and per-slot unit
attempt caps survive without refitting time. The invariant minimum and
D12 average preserve symmetry. At g0 the whole generator equals native HB.
This is an explicit countermodel family, not an independently sourced law.

The exact q3,N3 fixture uses W=Mult(n)*2^(sum n_i^2), the corresponding
conditional HB, and s=2/3. It preserves all drifts while changing centre
Gamma(p0,p0) by1/162. The q12 uniform algebra fixture uses s=1/16 and
changes centre total activity by1/192 and Gamma(p0,p0) by1/2592. The latter
nonzero uniform fixture is NOT the g0 member of the g0-preserving family.
A separate nominal A7/g3 fixture checks the construction numerically only.

## 3. What the common drift DOES fix: stationary averaged covariance

For any stationary finite generator, product expansion gives

    E_pi Gamma(f,g)=-E_pi[f Lg+g Lf].

Hence generators with the same stationary law and the same coordinate
drifts have identical equilibrium averaged coordinate covariance tensors,
although their pointwise tensors may differ. Reversibility is not needed
for this identity, but is retained for the licensed comparison class.
For a single observable, matching its exact drift alone suffices.

This explains why the positive centre change in the wheel must be compensated
elsewhere. The verifier pays all nine exact averaged q3 tensor identities.
An averaged covariance is therefore NOT the complete pointwise tensor
required by Theorem1 and does not reconstruct the generator.

## 4. Derived bridge lemma for stationary critical parity

Use the admitted actual MT domain D_N, N in12N, exact
g_N=g_c+lambda/sqrt(N), fixed bounded lambda and stationary pi_D,N.
Consider a diagnostic count generator Ltilde_N with:

- the same pi_D,N stationary, single-replacement jumps and total rate<=CN;
- exact reflected parity drift matching Ltilde_N m=L_HB,D,N m,
  where m=sigma^T p and sigma_i=(-1)^i.

These are class hypotheses, NOT newly derived physical laws. Let
G_N=N Gammatilde_N(m,m). Then0<=G_N<=4C and exactly

    E_piD G_N=E_piD[N Gamma_HB,D,N(m,m)] ->2.

The last limit follows from pi_D concentration at u, native HB proposals
tending to uniform and the paid vanishing departure-layer mass.

If also

    E_piD|G_N-E_piD G_N| ->0,                       (NC)

then E_piD|G_N-2|->0 and stationarity gives the growing-time payment

    E integral_0^S |G_N(p(sqrt(N)s))-2| ds ->0.

This is exactly the bracket scale for Z_N=N^(1/4)m(sqrt(N)s).
No mixing rate, iid sampling, independent environment or new noise is used.
NC is a sufficient analysis condition, not established for every reversible
drift-matched generator. A convenient stronger sufficient condition is

    lim_(delta->0) limsup_N sup_(p in D_N, ||p-u||<=delta)
                         |G_N(p)-G_N(u)|=0.

Boundedness and pi_D=>delta_u imply E|G_N-G_N(u)|->0; the averaged identity
then gives G_N(u)->2 and NC. No global boundary equicontinuity is required.

### A genuine non-HB class member already exists

Multiply only same-parity single-copy edges of native reflected HB by
eta(g)=1/(1+g^2); keep opposite-parity edges unchanged. This preserves
pi_D, D12/slot symmetry, positivity and the unit attempt cap; at g0 it is HB.
Each altered jump has delta m=0. Therefore parity drift AND Gamma(m,m)
are exactly unchanged. Other coordinates/rates change. It is a nontrivial
diagnostic member of the parity-matched class, not a new source for HB.
The exact uniform fixture tests the algebra at eta1/10, not its g0 member.

### Equicontinuity is not automatic, and its failure is not nonuniversality

For large N in12N place the D12-averaged wheel at the balanced type u.
All its finitely many vertices lie in the interior of D_N. Nearby W ratios
tend to1 and native single edge rates divided by N tend to positive values.
Thus s_N=epsilon N W_N(u), with sufficiently small fixed epsilon>0,
preserves positivity and per-slot caps for large N. It preserves even full
coordinate drift and reflected Gibbs reversibility.

For consecutive labels {0,1,2}, sum_ell(sigma^T r_ell)^2=16, also after
D12 averaging. Hence Gtilde_N(u)-G_HB,N(u)=16epsilon. Outside the finite
wheel support it is zero, at points only O(1/N) away. Local equicontinuity
fails. This is an analytical construction with existential thresholds,
not a numerical large-N simulation or a paid finite N0.

It is NOT a different critical-limit counterexample: the finite support's
pi_D mass tends to zero by the continuous joint scaled MT limit, and
bounded G then gives vanishing L1 difference from HB. It actually satisfies
NC. Do not confuse failure of a sufficient regularity condition with failure
of the proposed process theorem.

## 5. Exact handoff boundary

The inverse theorem, drift-preserving stress and stationary bracket lemma
are proved here. The FULL diagnostic-class critical path theorem is **not
claimed executed**: its actual drift remainder, maximal containment,
fourth/square-martingale UI and joint-filtration identification belong to
the next frozen research order. The native UF/MT result alone cannot be
renamed that class theorem. Removing NC is a further separate question.

No independent S07 activity source has been derived: a covariance calculated
from the postulated generator is circular as its source. The original
AS-002/003 physical bridge stays inactive. The present intrinsic mathematics
and the next agent's theorem need no external experiment or new microscopic
law, and do not promise physical QM, gravity, locality or ToE.
