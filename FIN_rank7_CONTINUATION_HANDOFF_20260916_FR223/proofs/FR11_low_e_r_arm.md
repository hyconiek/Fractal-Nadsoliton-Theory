# FR11 — low-odd-mass extension in the `r` direction

FR10 permits a large simultaneous box in `(u,v,e)` but its symmetric `r` radius
remains the original `1/8192`.  That limitation is caused by asking the same
box to tolerate `e` as large as `1/3072`.

On the physically important low-odd-mass slice, reuse the identical
second-order rational interval-AD / R7P-044 Schur-cone checker with

- `|r-r_*| <= 1/7600`,
- `u=1-exp(-3J4/2) <= 1/4800`,
- `v=1-exp(-J5/2) <= 1/4800`,
- `e=1-q_even <= 1/81,920`.

All sign prerequisites, `c2>0`, and both Schur endpoint tests pass with the
same margin `1/10000`.  Hence `lambda2(M4)<=sigma_*` on the whole box.

The same conservative checker rejects `e<=1/75000` with the other three
bounds unchanged; it also rejects `|r-r_*|<=1/7500` at the certified e bound.  This is retained only as a checker negative control, not as
a counterexample.
