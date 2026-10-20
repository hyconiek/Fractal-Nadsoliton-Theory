# Conditional consequences and adverse test of FIN-HYP-HIST-ME-001

Date: 2026-10-20. Self-contained finite-model proof, self-audited rather than
formalized or independently peer reviewed. The maximum-entropy method is
known. The contribution here is the frozen hypothesis, exact FIN specialization,
clock-free compatibility test and explicit failure/reference ledger.

## 1. A unique law really follows, conditional on the new postulate

Let S be finite, pi_x>0, sum pi=1. Let K be the symmetric 0/1 adjacency
matrix of a connected graph, **K_xx=1**. Histories are stationary S-valued
processes indexed by the integers, supported on K. No Markov property is
initially assumed. Their finite-alphabet entropy rate exists: by stationarity,
H(X_1|X_0,...,X_(-n)) decreases to a nonnegative limit h.

Every such law has stationary one-step flow c_xy=Pr(X_0=x,X_1=y), satisfying

    c>=0, c 1=pi, c^T 1=pi, c_xy=0 if K_xy=0.

Conversely, each such c defines a stationary Markov law P_xy=c_xy/pi_x.
Define

    F(c)=-sum_(K_xy=1)c_xy log(c_xy/pi_x).

Conditional entropy gives h<=F(c). The feasible flow set is compact and
convex. F is strictly concave on its variable entries, so it has a unique
maximizer c*. It is positive on every allowed entry: a feasible strictly
positive flow exists by assigning small equal symmetric masses to the
off-diagonal edges and the positive remainder to each diagonal. Mixing this
flow into a putative boundary maximizer gains an order -t log t term on
zeros, dominating the order-t changes on positive entries.

Since both margins equal pi, transposition preserves F and feasibility.
Uniqueness implies c*=(c*)^T; detailed balance is therefore a consequence
of this postulate and symmetric support, not an extra premise.

At the interior maximum, elementary Lagrange multipliers for the linear
row/column constraints give c*_xy=r_x s_y on allowed edges, with r,s>0.
Symmetry gives r_x/s_x=r_y/s_y along each edge. Connectivity propagates
this constant, so for one positive vector z,

    c*_xy=K_xy z_x z_y,
    pi_x=z_x (Kz)_x,
    P*_xy=K_xy z_x z_y/pi_x.                            (1)

Redundant row/column constraints do not obstruct multipliers: the gradient
is orthogonal to their common nullspace, hence lies in their row span.
The positive solution z exists by the preceding maximizer construction and
is unique since z_x=sqrt(c*_xx). Thus (1) is a well-defined source kernel.

The exact optimality certificate is, for every feasible c,

    F(c*)-F(c)=KL(c||c*)>=0.                            (2)

Indeed log c*_xy=log z_x+log z_y on allowed edges and both margins are
pi, so sum c log c*=2 sum pi log z=sum c* log c*. Positivity and the
equality condition follow directly from -log u>=1-u. For an arbitrary
stationary history law,

    F(c*)-h=KL(c||c*)+[H(X_1|X_0)-h].                  (3)

The stationary P* Markov law achieves equality. Equality for another law
forces c=c* and equality of all finite conditional entropies with the
one-step entropy. Their conditional mutual informations are zero, so every
finite past has the same next-state conditional as X_0 alone. All its finite
cylinders are therefore pi_(x0) product P*_(xk,x(k+1)), the same path law.
This proves uniqueness among histories, not just among Markov kernels.

Stationary finite cylinders are normalized and mutually consistent by row
stochasticity and pi P*=pi, yielding the ordinary finite-alphabet path law.
Connectivity, positivity on allowed edges and self loops give irreducibility
and aperiodicity. No physical clock or arbitrary-start preparation theorem
is inferred. The hypothesis selects a stationary source, not nature's pi.

## 2. The exact two-slot FIN reduction

Now S=(Z/qZ)^2, q>=3, with one self edge and edges changing exactly one
slot. Supply pi_(a,b)=p_(b-a)/q, p_d>0. For native N2,

    p_d=softmax_d[(g/2) A7_(0,d)],

is the inherited all-pairs Gibbs identity; it is not derived by the new
postulate. Simultaneous label translation preserves pi and K. Uniqueness
in (1) implies z_(a,b) depends only on d=b-a. Write
z_(a,b)=psi_d/sqrt(q), S_psi=sum_d psi_d.

For a fixed microscopic state of difference d there is **one** self
continuation and **two** distinct legal continuations to each v!=d,
one for each changed slot. Hence

    p_d=psi_d[psi_d+2 sum_(v!=d)psi_v]
       =psi_d(2S_psi-psi_d),                           (4)
    Pbar*_dd=psi_d^2/p_d,
    Pbar*_dv=2 psi_d psi_v/p_d, v!=d.                  (5)

The row totals to each difference fibre are independent of the absolute
labels. Thus the selected source is **strongly lumpable** to Pbar*, not
merely stationary-averaged. Its actual-change destination probabilities are

    J*_dv=Pbar*_dv/(1-Pbar*_dd)
          =psi_v/(S_psi-psi_d), v!=d.                  (6)

All off-diagonal entries are positive; eventual exit in intrinsic steps is
geometric and almost surely finite. This is a fully specified finite
history rule conditional on the authorized axiom and input data.

The positive solution of (4) also has a scalar representation. Put
R=2S_psi and t=R^2. Since psi_d<S_psi,

    psi_d=(R-sqrt(R^2-4p_d))/2,
    sum_d sqrt(1-4p_d/t)=q-1.                          (7)

The left side is strictly increasing for t>4 max p. At the lower endpoint
it is less than q-1; at t=4 it is greater than q-1 because
sqrt(1-p_d)>1-p_d for 0<p_d<1. Thus exactly one t lies in
(4 max p,4). This is an existence/uniqueness statement, not a numerical
root fit. The general finite source theorem already pays the same root.

## 3. Main adverse result: this is not native HB modulo a clock

For the supplied native N2 HB generator B_dv=2p_v, the first actual-change
law is J_HB,dv=p_v/(1-p_d). Scalar clocks change holding times, not the
odds between possible next destinations at the same prepared state.

For a fixed S_psi, x->x(2S_psi-x) is strictly increasing on (0,S_psi).
Consequently p_j>p_k implies psi_j>psi_k and

    J*_dj/J*_dk=psi_j/psi_k
       =(p_j/p_k)(2S_psi-psi_k)/(2S_psi-psi_j)
       >p_j/p_k=J_HB,dj/J_HB,dk,                       (8)

whenever d is distinct from j,k. For nonuniform p such a choice exists
since q>=3. Therefore the two source laws differ even modulo any
order-preserving onto time reparametrization with fixed origin and no event
erasure. This reuses IC's elementary first-exit invariant; it does not
compare raw spectral gaps after an arbitrary time fit.

For the inherited A7 row, p is positive and nonuniform at nonzero gain.
Thus (8) applies to the actual q12 N2 FIN comparison, conditional on the
admitted static input. No floating matrix experiment is needed for its sign.
At uniform p the branch laws agree; their self probabilities differ, so
this is only a zero-g clock-class control, not a positive-g recovery.

This is a **failure of HB compatibility**, not empirical falsification of
the new hypothesis: no physical data were supplied. We retain the failure
and do not change the frozen self convention to make HB fit.

## 4. A genuine restriction on the old activity freedom

The new source does make a falsifiable conditional choice: J*_dj/J*_dk is
independent of the current d whenever d differs from both destinations.
The mandatory distance-three thinning of native HB has quotient rates

    Bthin_dv=2 h_dist(d,v) p_v,
    h_3=1/(1+g^2), other h=1.

For q=12, use destinations j=0,k=1 and current states d=3 and d=4.
The two thinning odds are h_3 p_0/p_1 and h_3^(-1)p_0/p_1. Their double
ratio is h_3^2, whereas the new source predicts exactly 1. For g>0 these
are different, even after scalar clocks. At g=3 the thinning value is
1/100. This is a supplied-model arithmetic control, not a measured number.

Hence the axiom genuinely constrains a previously free kinetic direction.
It excludes that old alternative **and also disagrees with native HB's
nonuniform destination odds**. Selecting the former exclusion while hiding
the latter would be invalid. No uniqueness of nature's kinetics follows.

## 5. Why coarse-graining the entropy objective changes the answer

If the same bare Shannon-maximization phrase is instead applied directly
to the complete q-state **difference** carrier with one self continuation,
there is no one-slot support restriction left. The unique maximizer is
the iid law P_dv=p_v, since H(D_1|D_0)<=H(p), with equality precisely at
independence. Its change-of-state branch law is the HB one.

That is not the projection of the frozen labelled-carrier hypothesis.
The omitted combinatorial reference is explicit: for the uniform
microscopic lift of a macro row,

    H_micro(row | difference d)
       =H_macro(row | d)+sum_(v!=d)Pbar_dv log 2.      (9)

The selected source (5) is that uniform lift. Conversely every stationary
macro kernel with law p has a uniform two-choice microscopic lift, with
stationary law p_d/q; count the two predecessors per other fibre to verify
it. Its conditional entropy is bounded above by (9), with equality for
uniform microscopic choices. Therefore the full microscopic variational
problem reduces exactly to maximizing

    H_macro + (1-sum_d p_d Pbar_dd) log 2,              (10)

not plain H_macro. This is an exact reference-multiplicity payment, not
physical spatial composition or an optional correction chosen after fit.

For a general fixed symmetric reference K_xy>0 on the same support,
maximizing F_K(c)=F(c)+sum c log K selects
c*_xy=K_xy z_x z_y, with the same exact KL gap. Symmetric multipliers
K'_xy=b_x K_xy b_y do not change the optimizer: their log contribution
is the fixed constant 2 sum pi log b. Other reference changes need not be
gauge changes. In particular replacing one self continuation by two is
not of this gauge form on a graph with any off-diagonal edge: all diagonals
would require b_x=sqrt(2), changing those off-diagonal entries as well.

## 6. Non-adopted stress controls: no silent repair

Counting the two possible labelled **null attempt marks** as two distinct
self histories changes the microscopic reference to K_xx=2, K_xy=1 on
single-slot edges. The quotient reference then has every entry equal to 2.
The weighted problem gives Pbar_dv=p_v and microscopic off-diagonal
probabilities p_v/2: it recovers the N2 HB attempt kernel conditionally.
This is a different marked-history assumption, not a derived natural
reason for the primary frozen unmarked axiom. It is not adopted as a fix.

Nor does that success supply a general HB source for more slots. An exact
three-slot binary control has unnormalized Gibbs weights 8 on all-equal
states and 2 on mixed states. With uniform single-slot HB attempts,
for x=000,y=100 the stationary flows are

    c_xx=8/35, c_yy=1/35, c_xy=c_yx=2/105.

Any maximum-entropy flow with three marked self channels would require
c_xy=z_x z_y and c_xx=3z_x^2, c_yy=3z_y^2, hence

    9 c_xy^2/(c_xx c_yy)=1.

The HB control instead gives exactly 1/2. This eight-state rational model
is a stress test of a **general** marked-recovery claim, not an actual A7
N3 theorem. Do not extrapolate the special N2 recovery to the native family.

## 7. Exact positive fixture and scientific disposition

For q=4 choose a=(1,2,3,2), A=sum a=8 and D=2A^2-sum a^2=110.
Set psi_d=a_d/sqrt(D), p_d=a_d(2A-a_d)/D=(15,28,39,28)/110.
The resulting sixteen-state labelled source and its difference quotient
are exactly rational, positive, stationary and reversible. At current
d=2, destination 0 has probability 1/5 instead of HB's 15/71, difference
-4/355. This is a tiny exact fixture, not the numerical A7 law.

The authorized postulate thus supplies a unique stationary **intrinsic-step
source and branching rule conditional on its premises**. It is neither an
empty provider gate nor a derivation from old FIN. Its dependence on a
specified representation/reference and its adverse HB comparison are paid.
Strict F01/F03/F06, natural A7, carrier/composition, physical units/clock,
calibration/EV and every selector/legacy/QM/GR/SM/L_total/ToE gate stay open.
No physical adoption, model replacement or automatic axiom repair is made.
