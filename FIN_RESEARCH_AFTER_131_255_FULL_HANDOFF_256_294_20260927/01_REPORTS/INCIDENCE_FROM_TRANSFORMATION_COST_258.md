# INCIDENCE-FROM-TRANSFORMATION-COST-258
## Adjacency is operationally reconstructible once primitive pair interventions are given, but equal pair costs cannot select a unique graph

Date: 2026-09-27

Status:
- exact positive reconstruction theorem after task 256;
- exact symmetry no-go for sourcing incidence from current one-unit data alone.

## 1. Positive part: reconstruct adjacency without coordinates

After task 256 reconstructs factors

    Omega_1,...,Omega_m,

take any primitive pair intervention T.

Its operational support is the set of reconstructed factors whose quotient labels can change under T.

This definition uses only:
- the intervention action;
- the reconstructed partitions.

It is invariant under arbitrary global state recoding.

Declare an operational edge {i,j} when a primitive elementary intervention has support exactly {i,j}.

### Scrambled C4 replay

Use four hidden Z3 factors, hence 81 global states.

Globally recode all 81 states by

    i -> 7i+11 mod81.

From recoded local interventions task 256 first reconstructs four 3-state factors.

Then four recoded pair SWAP operations are supplied.

Their recovered supports are exactly

    {0,1}
    {1,2}
    {2,3}
    {0,3}.

So the C4 adjacency is recovered exactly after a monolithic non-product recoding.

No coordinates are used.

## 2. Negative part: current FIN data do not select which pairs are primitive

Suppose:
- m effective Z3 slots are indistinguishable;
- every pair admits the same SWAP gate;
- the gate cost is derived only from one-unit data.

Then slot-permutation symmetry S_m makes all unordered pairs cost-equivalent.

There are only two natural outcomes.

### Include every equally elementary pair

The incidence graph is

    K_m,

the complete graph.

This returns the old mean-field locality problem.

### Minimize total cost subject to connectedness

Every spanning tree has the same cost.

The number of labeled minimizers is Cayley's number

    boxed:
    m^(m-2).

Examples:

    m=4:   16 minimizers
    m=6:   1296 minimizers
    m=8:   262144 minimizers.

No unique geometry is selected.

## 3. Consequence for the planned transformation-cost programme

A minimum-cost incidence principle can work only if FIN supplies pair-dependent relational data, for example:
- a pair state;
- a barrier;
- a response cost;
- a conserved relational charge;
- or a restricted primitive-pair algebra.

The accepted one-unit Z3 generator and unique SWAP gate are insufficient.

## Verdict

Task 258 is a controlled NO-GO under the currently sourced data.

The problem is no longer "how to read a graph from operations".

That part is solved.

The blocker is:

    boxed:
    what in FIN makes one pair relation elementary/cheap and another non-elementary/expensive?

Until that is sourced, physical dimension-flow task 259 must not be promoted.
