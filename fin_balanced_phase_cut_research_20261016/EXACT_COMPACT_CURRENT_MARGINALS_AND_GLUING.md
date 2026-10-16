# A new typed proof object: nonnegative exact current marginals

Date: 2026-10-16. Analytical construction for compact critical levels of
the actual finite count model. It constructs **marginals**, not the missing
edge couplings, a full generator or an IBA counterfamily.

## 1. Current-carry ansatz with exact microscopic source

For phase P=A or B and parity level l define

    X_P,l=P intersect {m=m_l}, W_P,l=sum_(x in X_P,l)W_x,
    J_P,l=sum_(x in P,m<=m_l)s_x,
    I_P,l(x)=J_P,l-1 W_x/W_P,l,
    O_P,l(x)=I_P,l(x)+s_x.

The index on J is the cumulative **unnormalized exact native source**, not
the limiting integral. Therefore

    sum I_P,l=J_P,l-1, sum O_P,l=J_P,l.               (1)

For l and its neighbors in any fixed compact critical interval, the new
phase-prefix/local-slice theorem gives uniformly

    I_P,l(x)/(N W_x) -> J_P(z_l)/(4q_P(z_l)).         (2)

The left ratio is independent of the location x within the fibre.
There is a useful uniform bound on the limiting right side. Since
`3/8<=p_P(z)<=5/8`, evenness and the sign of b0 give

    (3/8)rho0(z)<=J_P(z)<=(5/8)rho0(z).

Thus `3/5<=J_P/q_P<=5/3`, and (2) lies in [3/20,5/12]. In particular,
for all sufficiently large N in the compact window,

    N W_x/8 <= I_P,l(x) <= N W_x/2.                 (3)

This conclusion is for the current-carry ansatz and static native source;
it is not a selected physical mobility or an adopted microscopic law.

## 2. Pointwise positivity away from the true departure layer

At g_c, write p=p_m+w in the actual cylinder. Since g_c A p_m=m sigma,
the unreflected mean-field parity drift on the parity line is tanh(m)-m.
For the fast perturbation, ||g_c A w||_infinity<=12||w||_2<12/10000.
The derivative of sigma.softmax(h+t delta) is a categorical covariance.
Its absolute value is bounded by the range product/4, at most
||delta||_infinity when sigma has range two. Hence

    |sigma.softmax(g_c A p)-m|
       <=epsilon^3/3+12delta <1/100,

with epsilon=1/100, delta=1/10000. The finite LOO correction is uniformly
O(1/N), so the native **actual reflected** drift obeys |d0(x)|<1/100 for
large N at every x outside E_N; there all proposals remain in D_N.
Therefore |s_x|<=N W_x/200 off E_N. From (3),

    O_P,l(x)>=3N W_x/25   for x not in E_N.          (4)

Reflection at E_N can violate this small drift estimate. It must be paid,
not discarded because its mass is small.

## 3. An exact, nonnegative departure correction

Within a compact phase fibre set

    a_x=[-s_x-I_P,l(x)]_+ for x in E_N intersect X_P,l,
    A_P,l=sum a_x, W_bulk=sum_(x in X_P,l\E_N)W_x.

Define I'_P,l by adding a_x at the departure states and subtracting
`A_P,l W_x/W_bulk` at every bulk state. Set O'_P,l=I'_P,l+s_x.
Only the incoming carry is redistributed; the exact source remains fixed.
Since |s_x|<=N W_x everywhere and I>=0,

    A_P,l<=N W(E_N intersect X_P,l).

The slice theorem gives a polynomial lower bound for pi(X_P,l), whereas
pi(E_N) is exponentially small. Uniformly on compact levels,
`W(E_N intersect X_P,l)/W_bulk ->0`. Thus the bulk subtraction is eventually
smaller than the positive margins in (3)--(4). At each departure state,
the very definition of a_x ensures O'>=0 and I'>=0. Consequently

    I'>=0, O'>=0, O'-I'=s_x exactly,
    sum I'=J_P,l-1, sum O'=J_P,l,
    ||I'-I||_1=2A_P,l<=2N W(E_N intersect X_P,l),
    I'_P,l(x)<=2N W_x, O'_P,l(x)<=3N W_x.           (5)

The constants are coarse; the N0 is asymptotic, not numerical. All
functions are G_+-invariant. This is an actual formula from existing
counts, weights and rates, with every departure state retained. The
construction does not turn an approximate equality into an exact one.

## 4. What coupling would finish on a band

On a whole parity fibre X_l use the phasewise I',O' just constructed;
write these as I'_l(x),O'_l(x). Suppose for every adjacent pair of compact
levels one has a nonnegative matrix c_l(x,y) on **all legal upward H edges**,
including the already allowed interphase edges incident to E_N, with

    sum_y c_l(x,y)=O'_l(x),
    sum_x c_l(x,y)=I'_l+1(y).                         (6)

Total masses agree by summing (1) over the phases. These phasewise prescribed
marginals require equal A-to-B and B-to-A interphase mass within each layer
pair; they do not require each mass to vanish. This is a restriction of the
carry ansatz, not a new condition imposed on every original feasible flow.
At every internal band vertex, the difference between
its outgoing and incoming conductance is exactly s_x, and the total incident
conductance is at most 5N W_x by (5). Thus (6) supplies the exact native
parity drift and a uniform total-rate cap **on that band**, by detailed
balance rates q_xy=c_xy/W_x. Top/bottom band carries remain prescribed
analysis boundary data, not reflecting physical controllers or a complete
closed process. They must be joined to the rest of D_N before any global
flow or temporal theorem can follow.

For these fixed marginals the precise next question is weighted Hall:
for every subset T of the **whole** left fibre X_l,

    sum_(x in T)O'_l(x)
       <=sum_(y in N_H(T))I'_l+1(y).                 (7)

Together with equal mass, (7) is necessary and sufficient for (6): attach
a source to left vertices with capacities O', give allowed middle edges
capacity exceeding the total mass, and connect right vertices to the sink
with capacities I'. Finite cut balance gives exactly (7). No edgewise
capacity needs to be added to (6): the already bounded row/column marginals
pay the original vertex cap. Native ellipticity then implies a corresponding
relative bound if desired.

Neither positive marginals nor nine available neighbors prove Hall for
all subsets. This is the missing coupling object, not a new expansion
postulate. Failure of this particular carry ansatz would reject the ansatz,
not arbitrary H-supported flows or universal IBA. Even success on |z|<=1
would still need far-level and actual-domain gluing before SC-002/003 can
be activated. Within-phase couplings would be sufficient, but are not required
by (6). Do not silently remove the true E_N interphase edges or turn a failure
of that stricter coupling into failure of (7) on the full frozen H support.
