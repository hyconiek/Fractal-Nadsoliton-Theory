# FP-002 — exact finite primal/dual and bounded all-N attempt

## Primal on the **actual** D_N (not a fictitious domain)

For each frozen phase P in {A,B}, retain **only** edges `{x,y}` of the actual legal count graph for which both endpoints are in P and `|m(y)-m(x)|=delta=2/N`. Same-parity edges have zero effect on scalar drift and can be set to zero in a nonnegative feasibility test. Orient each surviving edge `x->y` from lower to higher parity. Put `c_xy=W_x q_xy=W_y q_yx>=0` and `s_x=(N/2)W_x d0(x)`.

The exact equations are:

`sum_(e leaving x upward)c_e - sum_(e entering x downward)c_e=s_x` **at every vertex**;

`sum_(e incident x)c_e <= C N W_x` **at every vertex**.

All other rates vanish; the allowed diagnostic reflection is pure departure rejection. The unknown *single* finite C must hold for all sufficiently large N. The sums and rates must use the exact unrounded `g_c`, and all boundary vertices are included.

## Necessary phase cut and cap certificates

For phase P and adjacent occupied parity levels `a_l,a_(l+1)`, write `M_l^P=sum_(x in P,m(x)=a_l)W_x` and `J_l^P=sum_(x in P,m(x)<=a_l)W_x d0(x)`. Summing exact incidence equations gives

`K_l^P=(N/2)J_l^P>=0`;

`K_l^P<=CN min(M_l^P,M_(l+1)^P)`;

`K_(l-1)^P+K_l^P<=CN M_l^P`.

At any vertex, `d0(x)>0` requires a legal same-phase parity-*up* neighbour; `d0(x)<0` requires a legal same-phase parity-*down* neighbour. These local cone checks and the phase cut signs are necessary, not sufficient.

## Exact dual alternative (with all capacity vertices)

Define `A_{x,e}=+1` at the lower endpoint and `-1` at the upper; `B_{x,e}=1` at both endpoints. The primal is `Ac=s, Bc<=CNW, c>=0`. By finite-dimensional Farkas, infeasibility for fixed N and fixed C is equivalent to some arbitrary real vertex potentials `phi_x` and nonnegative `beta_x` satisfying, for **every retained phase edge** lower x to upper y,

`phi_x-phi_y+beta_x+beta_y>=0`, and

`sum_x phi_x s_x + CN sum_x beta_x W_x <0`.

(The witness may involve states away from the parity cut. A witness at one N does not rule out an all-sufficiently-large subsequence.) A proof of asymptotic infeasibility would need admissible witnesses uniformly excluding **every** fixed C on all sufficiently large N, or another exact argument. No such witnesses have been established here.

## Actual attempted resolution and the first remaining lemma

Phase-total source cancellation is exact by D12 rotation. The continuum parity-cut necessary currents from the admitted MT Gaussian are **strictly positive in both phases**, so they provide no macroscopic sign contradiction. A new cap lower bound (next document) is nontrivial, but finite and therefore does not contradict the permitted arbitrary finite C. The actual microscopic `Ac=s` and `c>=0` constraints near the moving fast-radius boundary, **including actual MT energy-boundary vertices**, remain unsolved. Neither the native q12 rate law nor a small artificial LP is substituted as a feasibility proof. 

**Scientific disposition: `OPEN_UNIFORM_FLOW_LEMMA`**. An all-N construction/certificate has NOT been produced. The frozen phase candidate has not been killed; it has also not been realized.
