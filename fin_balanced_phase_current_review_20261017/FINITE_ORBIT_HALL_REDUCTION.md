# Reusable finite tool for the remaining PC gate: orbit Hall deficits

Date: 2026-10-17. Specialization of the already admitted orbit-flow method,
not a new all-N Hall theorem or a new physical bridge.

The actual row r=O'_l, column c=I'_l+1 and H adjacency are G_+-invariant.
G_+ fixes each parity level and phase, and preserves the exact E correction.
Odd rotations generally exchange levels and are not automatically symmetries
of this fixed ordered fibre pair.

For a left subset T set delta(T)=r(T)-c(N_H(T)). This finite deficit is
supermodular: N(T union S)=N(T) union N(S), while N(T intersect S) is a
subset of N(T) intersect N(S). Thus

    delta(T union S)+delta(T intersect S)>=delta(T)+delta(S).

All translates of a maximizing T have the same maximum deficit. Unions
and intersections of maximizers remain maximizers by this inequality.
The union of all G_+ translates is therefore an invariant maximizing set.
Consequently **every Hall inequality holds iff every orbit-union left-set
Hall inequality holds**. A positive maximum has an invariant witness.

For an actual quotient coupling f_OP, use positive native reference u_xy
and U_OP=sum_(x in O,y in P)u_xy to lift

    c_xy=f_OP u_xy/U_OP.

Transitivity gives row U_OP/|O| and column U_OP/|P|, hence exactly the
original vertex marginals even with stabilizers and unequal orbit sizes.
Zero aggregate edges are absent. The previously derived marginal bounds
would pay the vertex cap; no edgewise cap is secretly substituted.

This identifies where the next proof should work: **dense orbit unions,
including departure tails**, not another local degree/overlap refinement.
Orbit sizes are at most twelve; the number of orbits is not bounded by
the 36 replacement pairs. No enumeration or uniform estimate for those
remaining dense sets has been supplied here.
