# Uniform native edge ellipticity on the *actual* local MT count domain

Status: **exact finite-N consequence of the frozen local cylinder and heat-bath rule**. Not a claim about full/global count space or physical kinetics.

Let `N` be an admissible multiple of 12, `D_N` the genuine MT connected sublevel in the cylinder, `x=(n_0,...,n_11)`, `p=n/N`, `sigma_j=(-1)^j`, `m=sigma.p`, `p_m,j=(1+m sigma_j)/12`, and `w=p-p_m`. The inherited cylinder satisfies `|m|<1/100` and `||w||_2<1/10000` at **every** vertex in `D_N`. Therefore, for every label,

`p_min := (1-1/100)/12 -1/10000 = 103/1250 < p_j`

and `p_j < p_max := (1+1/100)/12 +1/10000=158/1875`. These strict bounds use the **actual** MT cylinder, not a fabricated larger state space. The additional MT free-energy sublevel only restricts the domain further.

Let `A=A7` denote the fixed real 12x12 matrix, `g=g_c=12/Lambda6>0` exactly. Define the finite model constant

`Delta_A = max_{b,j,k in {0,...,11}} (A[b,k]-A[j,k]) >=0`,

and `a_* = (103/15000) exp(-g_c Delta_A) >0`.

For any two **distinct count states x,y in D_N** joined by an allowed native single-copy replacement `a -> b`, write `r=n-e_a`, so `r_k>=0`, `sum_k r_k=N-1`, and use the native **leave-one-out** rate with *no allowed-state renormalization*:

`q0(x,y) = n_a exp[(g/N)(A r)_b] / sum_j exp[(g/N)(A r)_j]`.

For every `j` we have `(g/N)[(Ar)_j-(Ar)_b] <= g Delta_A (N-1)/N <= g Delta_A`, because `r_k` are nonnegative. Consequently

`exp[(g/N)(Ar)_b]/sum_j exp[(g/N)(Ar)_j] >= (1/12)exp(-g Delta_A)`.

Because `n_a>=N p_min`, it follows **edgewise, uniformly in N and the location of the edge** that

`a_* N <= q0(x,y) <= p_max N`.

The same holds for the reverse move, since it also starts inside `D_N`. Thus for the exact native conductance `u_xy=W_x q0(x,y)=W_y q0(y,x)`:

`a_* N max(W_x,W_y) <= u_xy <= p_max N min(W_x,W_y)`

and `a_*/p_max <= W_y/W_x <= p_max/a_*`. The estimates apply to retained H edges and deleted F edges **when both endpoints are in D_N**. Edges rejected for leaving `D_N` are not asserted to have positive rates.

Summing original copy-attempt clocks gives `sum_{y in D_N, y != x} q0(x,y) <= N`, with outgoing proposals rejected rather than renormalized. The uniform ellipticity estimate is entirely compatible with `<=N` because the number of distinct count neighbours is at most 132, and `a_*` may be very small.

Because the *admitted versioned* A7 is positive semidefinite (`A7=XX^T`) and its largest eigenvalue is `Lambda6`, `||A7||_op=Lambda6`. Thus `g_c Delta_A <= 2 (12/Lambda6) Lambda6=24`, yielding the **fully explicit** universal lower bound `a_0=(103/15000)exp(-24)>2.59e-13` for every retained q12 MT edge. It is deliberately conservative; it does **not** need a rounded spectral value or a newly fitted constant. The exact `a_*` above is at least `a_0`. Neither constant is a new fitted quantity or an independent activity-provider law: both are inequalities for the **already supplied** native heat-bath generator.

**Scope separation.** The global count simplex allows `n_a=0` and has no corresponding edgewise positive uniform lower bound. Artificial toy graphs can likewise have native reference edges of order `1/N`, invalidating the conclusion. The result below relies fundamentally on the actual MT cylinder's positive occupation floor.
