# RB-001 — native source, leaky graph and exact finite cuts

**Status:** `EXACT_FINITE_GRAPH_CONTRACT_AND_NEW_CUT_DECOMPOSITION`; no all-N flow theorem.

## Frozen objects

For every sufficiently large `N in 12*N`, fix the *actual* connected MT count sublevel `D_N` containing the uniform type, with `p=n/N`, `m(n)=sum_j(-1)^j n_j/N`, and

`F_g(p)=sum_j p_j log(12p_j)-(g/2) p^T A7 p`, `g=g_c=12/Lambda6` **exact**.

The outer cylinder is `|m|<1/100`, `||p-p_m||_2<1/10000`, all `p_j>0`; the sublevel threshold is `F_g(p)<=10^-10`. The boundary/departure layer `E_N` consists exactly of points in `D_N` with at least one legal single-copy proposal outside `D_N`. All edges leaving `D_N` are **rejected**, not renormalized. `E_N` is not a physical controller. Exact membership of `D_N` is a mathematical predicate involving logarithms, not a binary64 approximation.

Let `C(n)=sum_j n_j cos(pi*j/2)` and `S(n)=sum_j n_j sin(pi*j/2)` (both integers), `R_N=2(C(n)^2+S(n)^2)/N`, `Z_N=N^(1/4)m(n)`, `gamma3=1-Lambda3/Lambda6>0`, `p(z)=1/2+exp(-z^2)/8`, `T(z)=-2 log(1-p(z))/gamma3`. Phase `A_N` includes ties `R_N<=T(Z_N)`; phase `B_N=D_N\A_N` includes strict `>`. The 1-label cyclic permutation maps `(m,C,S)->(-m,-S,C)` and preserves `F`, `R`, both phases and `E_N`.

Let `U_N` denote all **native** undirected opposite-parity count edges `{x,y}` in `D_N`. Orient from smaller to larger `m`; all increments are `delta=2/N`. `H_N` retains precisely edges of `U_N` whose endpoints have the same phase **or** at least one endpoint lies in `E_N`. Same-parity edges can be retained in the full diagnostic generator; for the scalar-drift nonnegative feasibility problem their conductances may be zero without loss. No multi-copy jump, proxy domain or new noise is admitted.

The full native reflected count rates use actual `r=n-e_a`, legal `y=n-e_a+e_b`:

`q0(n,y)=n_a exp[(g_c/N)(A7 r)_b] / sum_{j=0}^{11}exp[(g_c/N)(A7 r)_j]`.

Define exact weights `W_n=(N!/prod_j n_j!)exp[(g_c/(2N)) n^T A7 n]` and native drift `d0(n)=sum_{y in D_N}q0(n,y)[m(y)-m(n)]`. The **source** is `s_n=(N/2)W_n d0(n)` at **every** vertex, including the departure layer. The original native edge conductance `c0_xy=W_x q0(x,y)=W_y q0(y,x)>0` is fixed by detailed balance and the constant diagonal of A7. Then `s=A_U c0` exactly, where `A_U` has `+1` at an edge's lower endpoint, `-1` at the upper. The source is *not* an inferred cubic drift.

## Exact order-ideal defect theorem (new bounded result)

A subset `L subset D_N` is a **lower ideal of H_N** when `y in L` and `x->y` in H_N implies `x in L`. Let `F_N=U_N\H_N` denote the forbidden, interior-to-interior cross-phase edges. Define `cut_up(E,L)=sum_{x in L, y notin L, x->y in E}c0_xy`, and `cut_down(E,L)=sum_{x notin L,y in L,x->y in E}c0_xy`.

By summing the exact native continuity equation over `L`, all internal edges cancel. Because a lower ideal has **no incoming H-edge**, we have the exact identity

`sum_{x in L}s_x = cut_up(H_N,L) + cut_up(F_N,L) - cut_down(F_N,L)`.

All three summands use the native, strictly positive, **fixed** conductances. Consequently `sum_L s>=0` is equivalent to the *native defect test*

`cut_down(F_N,L) <= cut_up(H_N,L)+cut_up(F_N,L)`.

No conductances for the sought diagnostic generator are assumed in this identity. It isolates precisely how deleting interior phase-crossing edges can generate a negative ideal source. In particular, checking the sign of the supplied full microscopic `s` can be reduced to this exact cut inequality, rather than to an approximate expansion of `d0`.

**Automatic global parity-cut positivity:** For `L={x:m(x)<=ell}`, no native upward edge enters `L`, hence `cut_down(F_N,L)=0`. Therefore `sum_L s=cut_up(U_N,L)>=0` **for any support deletion**. These ubiquitous cuts cannot test leaky flow feasibility. Every remaining potential violation must use a more irregular ideal which selects a strict subset of some parity level (or several levels).

**Component payment:** Every undirected component `K` of `H_N` is simultaneously a lower ideal and an upper ideal, hence existence of an exact flow requires `sum_{x in K}s_x=0`. Connectivity of the *native* `D_N` is not connectivity of `H_N` and does not prove this. Rotating a lower ideal by one label and taking its complement gives another lower ideal with the same signed source; the D12 symmetry does not force all ideals nonnegative.

## Scope caveats

The user-supplied FP review states that opposite-parity moves in both directions remain available off `E_N`, with a coarse geometric bound `N>=10^8` and a further unknown MT `N0`. I did **not** obtain `BOUNDARY_LEAKAGE_AND_IDEAL_CUT_SEED.md` as original bytes and therefore do not independently promote that directional-availability estimate. Even if paid, local positive neighbours do not imply global feasible routing.

The usual MT proof shows `E_N` lies in an energy layer `F_g(p)>=h/2`, so its stationary mass is exponentially small after polynomial factors; this is inherited. Such rarity does not by itself build nonnegative within-phase conductances or prove a uniform cap.
