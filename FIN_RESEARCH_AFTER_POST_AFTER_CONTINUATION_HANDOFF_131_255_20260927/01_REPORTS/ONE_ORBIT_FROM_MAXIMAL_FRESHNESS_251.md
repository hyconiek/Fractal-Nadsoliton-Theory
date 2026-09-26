# ONE-ORBIT-FROM-MAXIMAL-FRESHNESS-251
## For a fixed finite environment capacity, one cycle uniquely maximizes the exact heat-bath reset horizon

Date: 2026-09-26

Status:
exact permutation theorem.

Report 250 used one cyclic slot permutation.

Now remove that assumption.

Let P be ANY permutation of L simultaneous environment/system slots.

Observe one distinguished slot s.

Let

    ell(s)

be the length of the permutation cycle containing s.

Under repeated application of P, the observed slot visits exactly

    ell(s)

distinct initial slot contents before repeating.

If all other visited slots were prepared as independent uniform ancillas, then
the exact fresh-reset horizon is

    boxed:
    ell(s)-1.

## 1. Maximum possible horizon

Since

    ell(s)<=L,

the longest possible exact reset window for fixed slot capacity L is

    L-1.

Equality holds iff the observed slot lies in a cycle of length L.

But a length-L cycle contains every slot.

Therefore:

    boxed:
    maximal use of finite environment capacity
      iff
    P is one transitive L-cycle.

This supplies an operational reason for the "one orbit" condition of reports
183-184.

It is no longer merely connectedness chosen for geometric convenience.

It maximizes how long a finite reversible environment can impersonate an
irreversible MaxEnt bath before recurrence.

## 2. Exact counting

The number of L-cycles on L labelled slots is

    (L-1)!.

The fraction among all L! permutations is

    1/L.

Exhaustive enumeration for L=2,...,7 reproduces this exactly.

So one-cycle structure is not statistically generic.

It is selected by the freshness-horizon objective.

## 3. Epistemic boundary

FIN has not yet derived the principle

    "for fixed environment capacity, maximize exact fresh-reset duration."

But unlike an arbitrary graph axiom, this extremal criterion is directly tied
to a previously established physical/mathematical task:

    reproduce the accepted one-unit heat bath with globally reversible
    information bookkeeping.

That makes it a concrete source candidate for transitivity.
