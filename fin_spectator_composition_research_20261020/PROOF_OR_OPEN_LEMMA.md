# SC: composition, hidden steps and autonomy of the selected source

Date: 2026-10-20. Self-audited conditional finite-model mathematics. The
maximum-entropy source and RT closure criterion are inherited, not derived
again as new postulates. No independent peer review or physical source claim.

## 1. Projection with uniform hidden coordinates has an exact source law

Let E have q>=2 labels. On E^N define K_N(x,y)=1 if x=y or x,y differ
in exactly one coordinate, otherwise zero. Each row has degree

    d_N=1+N(q-1).

Retain M>=1 coordinates, hide L=N-M>=1, and supply the static law
pi_N(x,h)=pi_M(x)/q^L, with arbitrary positive pi_M. This product is a
formal decoupled-state test of the generic hypothesis; it is not asserted
to be the interacting native A7 N-body equilibrium.

Every retained-state fibre has q^L elements and the given law is uniform
inside it. From each fine state, a different visible one-slot neighbour
has exactly one continuation. The same visible state has one true self
continuation plus L(q-1) changes of hidden coordinates. RT therefore gives
the exact induced reference and selected visible kernel

    Kvis=K_M+c I, c=L(q-1),
    projected source = P*_(K_M+cI,pi_M).                 (1)

The output is strongly lumpable and Markov. Internal steps were summed,
not rejected, renormalized away or exposed as extra observations. Their
extra diagonal is **derived from the unchanged fine model**, not a marked-
self repair of the primary postulate.

The isolated M-coordinate source is instead P*_(K_M,pi_M). These are two
different model-selection tasks; comparing them tests autonomy. RT already
paid consistency of descriptions of one fixed fine world. It did not say
that deleting coordinates and selecting a new smaller world commute.

## 2. When can the added hidden-step weight be only a clock change?

The following finite lemma is useful. Let K have constant diagonal 1,
constant row sum d, a connected off-diagonal graph containing an odd cycle,
and nonnegative symmetric entries. Fix positive pi. Compare its selected
kernel P with that for K+cI, c>0. Their actual-change destination kernels
are identical **if and only if pi is uniform**.

Proof. Both sources are reversible at the same pi. Write their stationary
flows as f_xy=K_xy z_x z_y and f'_xy=K_xy u_x u_y off the diagonal.
Equality of actual-change destination laws gives f'_xy=alpha_x f_xy,
where alpha_x is the ratio of the outgoing off-diagonal masses at x.
Symmetry of both flows and connected support force alpha_x=alpha_y=alpha
along every off-diagonal edge. Thus r_x r_y=alpha there, r_x=u_x/z_x.
An odd cycle, followed by connectivity, forces every r_x=sqrt(alpha).

The two marginal equations now imply

    pi_x=alpha[pi_x+c z_x^2],
    z_x^2/pi_x=(alpha^(-1)-1)/c=constant.               (2)

Consequently Kz=lambda z for one positive vector z. At its maximum and
minimum coordinates the regular row sum gives lambda<=d and lambda>=d.
Equality at a maximum propagates through all positive edges, so z is
constant and pi is uniform. Conversely uniform pi gives constant scaling,
P=K/d and P'=(K+cI)/(d+c), hence

    P'=I+[d/(d+c)](P-I).                               (3)

They differ only by independent idle steps and have the same change order.
No physical waiting-time law is supplied by (3).

For q>=3, the off-diagonal graph of K_M has triangles (change one coordinate
among three labels) and is connected and regular. Therefore (1) has the
same actual-change law as the isolated model exactly when pi_M is uniform.
For nonuniform pi_M, a scalar clock cannot repair the discrepancy in next
destinations. This is a failure of an **autonomous-component test**, not
empirical falsification or an assumption that the universe has fundamental
independent parts.

The odd-cycle condition is essential. In two states there is only one next
different state, so every positive Markov chain has the same alternating
jump order, even at a nonuniform stationary law. Do not extend the iff
statement to arbitrary bipartite off-diagonal supports.

## 3. Product equilibrium can also create a non-Markov visible record

For two q-state coordinates supply pi(i,a)=p_i nu_a, with p,nu positive.
Retain the first coordinate. Within-fibre weights are w_(i,a)=nu_a, with
fibre masses one. The RT all-positive-marginals criterion says stationary
Markov output at any p is equivalent to closure at the auxiliary base
p_i=1/q. This is an analysis normalization, not a physical preparation.

At that base, visible-label permutation symmetry makes the positive scaling
z_(i,a)=u_a. With S=sum_a u_a, its equation is

    nu_a/q=u_a[S+(q-1)u_a].                            (4)

For a different visible j there is one legal fine move (i,a)->(j,a),
whose probability is

    b_a=q u_a^2/nu_a.                                  (5)

Strong visible closure requires b_a to be independent of a. If so,
u_a=t sqrt(nu_a) for t>0; substituting in (4) gives

    1/q=t^2[(sum sqrt(nu))/sqrt(nu_a)+q-1],

which forces nu to be uniform. Conversely a uniform nu makes the base
source uniform on its regular graph and closure follows; RT then extends
it to every p. Hence, for q>=2,

    visible stationary history is Markov iff hidden nu is uniform.       (6)

For q>=3 there is also a direct check without the RT family argument.
Strong closure of each visible off-diagonal row makes
z_(i,a) z_(j,a)/nu_a independent of a for every i!=j. The visible complete
graph has triangles, so positivity forces z_(i,a)=s_i sqrt(nu_a).
The scaling equation then reads

    p_i=s_i(sum_j s_j-s_i)+s_i^2 (sum_b sqrt(nu_b))/sqrt(nu_a).

Its right side can be independent of a only at uniform nu. Uniform nu
gives the converse by exact block counts. This confirms the q12-relevant
necessity directly; the RT argument also covers q=2.

This concerns the selected source, not a generic theorem that product
equilibria forbid autonomous dynamics. Earlier replication packages supply
compatible families after their laws have been chosen; the present question
is whether this particular new variational postulate selects one.

In particular the statically coherent product family pi_N=nu^tensorN
fails Markov projection already at two->one if nu is nonuniform. Uniform
static marginals are the exceptional coherent case below. A stochastic
count-deletion kernel alone is not used to infer tagged-history Markovity.

## 4. Exact memory and an essential clock-only countercontrol

At p=uniform, put bbar=E_nu b. The compressed one-step visible kernel
has off-diagonal bbar and diagonal 1-(q-1)bbar. The reversible compression
Gram identity gives exactly

    tr(C2-C1^2)=q^2(q-1) Var_nu(b).                    (7)

For a current visible state, off-diagonal residuals are b_a-bbar and its
diagonal residual is -(q-1)(b_a-bbar). Sum their squares divided by the
uniform visible weight 1/q, then average over pi(i,a)=nu_a/q to obtain
(7). By (4)--(5), the variance is positive for nonuniform nu. This is a
finite stationary-history memory witness, not just failure of a worst-case
initial-state test.

Yet in this same case the **actual-change order is uniform**: from any
hidden a all q-1 visible destinations have the same probability b_a in
each step. The constant next-exit probability 1/(q-1) per destination solves
the killed-fibre hitting equations. Since min_a b_a>0 on the finite set,
nonexit tails are bounded by [1-(q-1)min b]^n, paying uniqueness and eventual
exit. The source Markov property then makes successive visible-change
destinations uniform given the full visible past. Thus non-Markov timing
does not by itself prove a clock-free branching discrepancy. No new
physical stochastic clock is adopted by this countercontrol.

For an exact nine-state fixture choose q=3, a=(1,2,3), A=6 and

    nu=(1/8,5/16,9/16), p=(1/3,1/3,1/3),
    z_(i,a_index)=a_(a_index)/sqrt(192).

The unmarked source is rational: hidden changes to b have probability
a_b/[6+2a_a], and visible changes to a different i have probability
a_a/[6+2a_a]. Its stationary law is p times nu, and it satisfies the
positive scaling equations for the original unweighted reference. Then

    b=(1/8,1/5,1/4), E_nu b=7/32,
    Var_nu(b)=9/5120,
    tr(C2-C1^2)=81/2560>0,                             (8)

while the next different visible label has probability exactly 1/2.

## 5. An exact clock-free autonomy failure in the transpose fixture

Swap the roles of the two coordinates: visible p=(1/8,5/16,9/16), hidden
nu=uniform. Equation (1) gives a Markov visible source with reference
K_1+2I=ones+2I and scaling psi=(1,2,3)/8. Its visible transition matrix is

    [[3/8,1/4,3/8],
     [1/10,3/5,3/10],
     [1/12,1/6,3/4]].

At visible state 2, its next different state is 0 with probability 1/3.
The isolated one-coordinate model has complete reference ones and selects
iid p; its corresponding probability is p_0/(1-p_2)=2/7. Their difference
is exactly 1/21. This discrepancy cannot be changed by any order-preserving
clock with fixed origin and no event erasure.
With the same unconditioned stationary visible preparation p, the joint
record (initial label, next different label) therefore has TV at least
p_2/21=3/112. This is not a postselected starting-state guarantee or a
physical finite-time performance bound; eventual exit is paid by the
finite positive kernels.

This uses the hypothesis's **own isolated law**, not an inserted HB target.
The static measure still factorizes: -log pi is a sum of separate-coordinate
functions. No physical conserved-energy interpretation is asserted. What
couples the update choices is the global entropy-per-step rule.
This is a mathematical consequence of that rule, not a contradiction of
the informational-whole ontology. Dynamic autonomy, if required, would
need a separate justification; it is not adopted as a new axiom here.

## 6. Positive case: uniform product families and count deletion

At pi_N=uniform on E^N the reference is regular and the selected source is
P_N=K_N/d_N. Forgetting N-M coordinates gives the Markov source

    P_(N->M)=I+[d_M/d_N](P_M-I).                       (9)

Visible legal changes each have probability 1/d_N; hidden changes add
(N-M)(q-1)/d_N to the visible self probability. This proves (9) directly
and pays full-history closure, not merely marginal agreement.

The corresponding count kernel has diagonal 1/d_N and off-diagonal
n_a/d_N on n->n-e_a+e_b. Let

    D_(N,M)(n,m)=prod_a binom(n_a,m_a)/binom(N,M)

be the actual sampling-without-replacement count map. In this uniform case,

    P_N,count D_(N,M)
      =D_(N,M) {I+[d_M/d_N](P_M,count-I)}.              (10)

To prove (10), retain a uniformly chosen M-slot subset fixed once before
evolution, independently as an analysis selection. It stays a subset of
the same parent history, not a fresh draw at each time. Its count law at
one time is D_(N,M); its future source is (9) by the paid labelled closure.
Conditioning on the initial count proves (10). This argument is not licensed
for a nonclosed tagged source simply because the count sampling formula
exists. A small exact N=3,M=2,q=3 fixture checks the stochastic identity.

Uniform compatibility requires only an idle-step normalization, whose
overall scale remains free. It is not a derived physical clock or per-slot
attempt law. At nonuniform product nu, (6) rejects full labelled Markov
projectivity even though the equilibrium family itself is consistent.
No nonuniform count-sampling identity is asserted without its own proof.

## 7. Scientific disposition and next mathematical input

The new results distinguish: consistency of equilibrium, closure of a
projection, autonomy under selecting a smaller carrier, and equality modulo
clock changes. These four properties cannot be identified. The unchanged
hypothesis has coherent uniform limits but generically permits kinetic
context dependence even in product-equilibrium tests.

These diagnostic product measures are not the actual positive-g A7
all-pairs Gibbs law; primary FIN is unchanged. The previous HB failure is
not repaired by a larger-N search. No physical subsystem, bath, clock,
independence principle or empirical falsification is added.

Equation (1) supplies a concrete next input for a distinct analysis of many
uniform hidden coordinates: the fixed visible reference K_M+cI. Any future
claim about a scaled continuous-time limit must pay the scaling equation,
normalization, preparation and actual path limit; it is **not proved by
this packet** and must not be supplied by fitting HB. Strict F01/F03/F06
and all physical/ontology gates remain open; the broad goal stays active.
