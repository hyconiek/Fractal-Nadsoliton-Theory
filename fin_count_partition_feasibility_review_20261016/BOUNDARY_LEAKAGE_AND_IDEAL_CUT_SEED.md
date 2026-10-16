# A different construction gate: rare boundary leakage and order-ideal cuts

Date:2026-10-16. Derived geometry and IF-feasible temporal implications.
**No closed or leaky admissible generator sequence is constructed here.**
The original frozen FP closed-flow problem remains open, not rewritten.

## 1. Phase membership alone cannot create a directional dead end

Retain the exact original phase functions Z_N, R_N and T. Let v=(C,S),
L=|C|+|S| and tau=T_min=2 log2/gamma3>1. From gamma3>4/25 and
|z|exp(-z^2)<=1,

    |T'(z)|<=4/(3gamma3)<=25/3=:K.

For EITHER desired parity increment +/-2/N, the four integer fast changes
w=(+/-1,+/-1) are all realizable as single-copy substitutions: even k3
labels have Fourier values(+/-1,0), odd labels(0,+/-1), and all source
labels are occupied in the actual MT positive cylinder. Exactly

    Delta Z=+/-2N^(-3/4),
    Delta R=4(v.w+1)/N.

Choose w inward with v.w=-L or outward with v.w=L. At either sign of the
parity change the same choices are available.

For phase A, if R<=tau/2, the postjump radius is at most
tau/2+4sqrt(tau/(2N))+4/N<tau for large N, hence remains in A. If R>tau/2,
L>=sqrt(N tau)/2 and the inward choice has Delta R<=-N^(-1/2) eventually,
whereas |Delta T|<=2K N^(-3/4). It therefore remains in A, including ties.
For phase B, R>T>=tau and the outward choice increases R by more than
N^(-1/2), dominating the same threshold change, so remains in B.

The coarse bounds N>=10^8 suffice for these geometric inequalities:
N^(1/4)>=100 and 2K/100<=1/6. This is NOT a numerical MT N0; the actual
MT theorems may require a much larger unknown threshold.

Thus any count state whose full single-copy neighborhood stays in D_N
has a same-phase neighbor in BOTH parity directions for sufficiently large
N. The phase surface by itself does not yield a singleton sign obstruction.
Any lost direction in the original phase graph must involve the actual
MT energy-domain departure layer or a nonlocal cut, not just the smooth
fast-radius cutoff. No zero-label population is silently included.

## 2. Explicitly change the auxiliary graph, not the parent law

Let E_N be the actual native departure layer inside D_N: states with a
legal original single-copy neighbor outside D_N. Its admitted Gibbs mass
satisfies pi_D(E_N)<=poly(N)exp(-bN).

Define a NEW diagnostic support graph H_N:

- retain all legal within-phase count edges;
- additionally allow cross-phase edges if at least one endpoint is in E_N;
- retain the original D_N and forbid every departure from D_N.

The old FP target required ALL cross-phase rates zero. This new route is
not that old target and does not declare it solved. It remains a subgraph
of the existing count graph and a subclass of the original KR/BA kinetic
class. No physical graph, hidden phase state, source law or controller is added.

At a vertex in E_N all native internal edges are now allowed, so any native
positive/negative drift has an up/down neighbor as required. Outside E_N,
the phase lemma supplies both directions. Hence this graph has no elementary
one-vertex drift-sign obstruction. This is NOT proof of a simultaneous
nonnegative flow: edge reversibility couples different vertex demands.

Seek symmetric c_xy>=0 supported on H_N with EXACT original parity drift
and one uniform total CN cap. Neither those rates nor a uniform capacity
bound have been derived by the geometry lemma.

## 3. Why rare leakage would still preserve the intended IBA obstruction

Assume such a family exists and starts in exactly pi_D,N. Reversibility
and the rate cap bound its equilibrium cross-phase jump intensity by

    lambda_cross=2 sum_(cross edges)c_xy <=2CN pi_D(E_N).

Over model time sqrt(N)S the chance of a phase change is at most
2CS N^(3/2)pi_D(E_N)->0. This is a compensator/Markov bound, not iid sampling.

Exact conditional stationarity of a phase is no longer true and must NOT
be reused unchanged. Instead reversibility gives

    E[1_A L m^2]=E[m^2 L1_A],
    |N E[1_A L m^2]|<=N lambda_cross<=2C N^2 pi_D(E_N)->0.

The stationary product identity and the already paid weighted native drift
error therefore give the SAME asymptotic two conditional bracket means
as the old closed-phase seed. At lambda0, with
a=1/2+E exp(-Z^2)/8, b=1-a and M=E[Z^2exp(-Z^2)]>0,

    E[r_N|A] ->-M/(2a), E[r_N|B] ->M/(2b).

Let Y_N=integral_0^S r_N(X_sqrt(N)s)ds. Since |r_N|<=4C, replacing
1_A(X_0) by 1_A(X_sqrt(N)s) inside its expectation has error at most
4CS times the phase-change probability, tending to zero. Exact global
stationarity then yields E[1_A(X_0)Y_N]=S E_pi[1_A r_N]+o(1), and similarly
for B. Jensen conditional on the INITIAL phase gives

    liminf E sup_(t<=S)|integral_0^t r_N ds| >= S M>0.

The initial law is the full unpostselected Gibbs law, not a conditioned
preparation. This implication requires the as-yet unconstructed rates.
As before, changes restricted to a vanishing-mass set cannot realize the
required O(1) phase bracket means under a uniform cap; bulk rate changes
are still needed. Adding rare native edges alone is not a flow construction.

The previously proved N^-2 native-mixture coupling may conditionally extend
a successful leaky family to strictly positive irreducible rates. It is
not used to manufacture missing flows or a uniform lower activity bound.

## 4. A sharper first gate: capacity-free lower ideals

Orient only the parity-changing edges of H_N upward in m. Same-parity
edges may be omitted for nonnegative drift feasibility. Set

    s_x=(N/2)W_x d0(x), (Ac)_x=outflow_x-inflow_x.

Call L a lower ideal if y in L and an allowed upward edge x->y imply x in L.
There is nonnegative c with Ac=s, IGNORING upper capacities, if and only if

    sum_x s_x=0, and sum_(x in L)s_x>=0 for every lower ideal L.

Necessity is the cut balance: lower ideals have no incoming edges. For
sufficiency the capacity-free finite Farkas alternative tests every real
phi with phi_x>=phi_y on upward edges. Any such phi is a constant plus
a positive linear combination of indicators of lower ideals (its upper
level sets). The displayed cut inequalities therefore pay every dual test.
This is a finite support criterion, not an all-N proof for H_N.

A single negative ideal is a capacity-INDEPENDENT dual obstruction, with
phi=1_L and beta=0. Ordinary cuts m<=a are only a small subset of these
ideals; their positive continuum limits do not suffice. The geometry lemma
rules out the simplest singleton sign witnesses, not larger fibre-dependent
ideals. Invariant components and their complements are included in the tests.

Even a capacity-free solution at every N may require caps C_N diverging
with N. A separate uniform CN bound is mandatory before any counterfamily
or temporal conclusion is promoted. A proved obstruction to this enlarged
graph also obstructs its strictly closed subgraph with the SAME N/cap
quantifiers. A merely failed proof attempt does neither. No such obstruction
settles universal IBA, and no construction or impossibility is asserted here.
