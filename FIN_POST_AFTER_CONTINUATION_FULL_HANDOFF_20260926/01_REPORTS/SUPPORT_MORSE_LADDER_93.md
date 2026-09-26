# SUPPORT-MORSE-LADDER-93
## Two independent branch networks realize the same 2-label / 3-label / 4-label Morse ladder

Date: 2026-09-26

The large-g support theorem predicts

    index = |S|-1

for affinely independent strict supports.

Reports 90-92 now show that this is not merely an endpoint classification:
actual finite-g symmetry events connect supports of successive complexity.

### d=1 network

Pair:
    {0,1}
    support size 2
    index 1
    stabilizer Z2.

Reflection-breaking daughter:
    {0,1,9}
    support size 3
    index 2
    trivial stabilizer.

Fold-return sheet:
    {0,1,4,9}
    support size 4
    index 3
    reflection stabilizer.

### d=6 / k4×k6 network

Pair:
    {0,6}
    support size 2
    index 1
    stabilizer order 4.

Transverse daughter:
    {0,3,6}
    support size 3
    index 2
    stabilizer order 2.

Pure-k4 parent:
    {0,3,6,9}
    support size 4
    index 3
    larger k4 symmetry.

So in two unrelated symmetry sectors, increasing asymptotic support complexity
by one corresponds to gaining one negative covariance/Hessian direction.

This does NOT establish a universal law that every local bifurcation must add
exactly one support label.  The d=2 family is already a counter-warning: it
enters a near-degenerate k5 angular cluster rather than displaying the same
simple genealogy.

The robust statement is:

    support cardinality fixes asymptotic Morse index;

the new empirical structural observation is:

    several concrete finite-g symmetry events connect those support strata in
    the expected index order.

That observation is now strong enough to motivate a support-incidence graph:
nodes are strict zero-temperature supports modulo D12, and edges are finite-g
fold/symmetry-breaking components.
