# FINITE-DILATION-OPTIMALITY-282
## Minimum capacity and maximal freshness make the cyclic register Pareto-optimal, but not unique without persistent record identity

Date: 2026-09-27

Status:
- exact entropy lower bound;
- exact coordinate-routing theorem;
- exhaustive q=2,m=2 counterexample to unrestricted uniqueness.

Goal:
test whether the finite cyclic register can be promoted from a natural construction
to a unique finite reversible approximation of the fresh heat bath.

## 1. Information lower bound

To generate m exact independent q-state fresh outputs from a fixed subsystem state,
the environment must carry entropy at least

    m log q.

Hence

    boxed:
    |E| >= q^m.

A register of m q-state environment records saturates this lower bound.

So it is capacity-optimal.

## 2. Freshness upper bound at minimum capacity

At minimum capacity the environment contains exactly enough information for m
independent q-state outputs.

No exact reversible construction can provide more than m independent fresh q-state
outputs without reusing information.

Thus freshness horizon m is maximal.

The cyclic register attains it.

## 3. Zero-content-rewrite objective

If elementary evolution only transports record contents and never rewrites their
values, the rewrite cost is nonnegative and the cyclic coordinate shift has cost zero.

Thus the cyclic register simultaneously achieves:

    minimum environment capacity;
    maximum exact fresh-reset horizon;
    zero content rewrite.

It is a Pareto-optimal construction.

## 4. Unrestricted uniqueness FAILS

Exhaustive test:

    q=2,
    m=2,

so the global system has

    2^(m+1)=8

states and there are

    8! = 40320

time-homogeneous bijections.

Requirement:
for each fixed initial system bit, the two observed outputs after steps 1 and 2
must run bijectively over all 4 possible fresh bit pairs as the minimal environment
varies.

Number of admissible bijections:

    boxed:
    9216.

So:

    minimum capacity + maximal freshness

is nowhere near unique.

## 5. Even multiset preservation is not enough

Among those 9216 maps, require the global multiset of three record values to be
preserved at every update.

Survivors:

    boxed:
    4.

But only:

    boxed:
    2

are fixed coordinate 3-cycles.

The other two are state-dependent content-preserving permutations.

Therefore the condition

    "do not rewrite values"

does not by itself give a fixed geometric scheduler.

## 6. Persistent record-identity theorem

Now strengthen the transport premise:

    elementary evolution permutes PERSISTENT RECORD IDENTITIES
    by one fixed coordinate permutation, independent of their contents.

Then the fresh-reset horizon equals

    cycle_length(observed record)-1.

With m+1 total records, achieving the maximal horizon m requires the observed record
to lie in a cycle of length m+1.

Therefore the coordinate routing is a single (m+1)-cycle.

All such cycles are conjugate by relabeling record identities.

Hence:

    boxed:
    within the fixed-record-routing class,
    the maximal-freshness topology is uniquely C_(m+1) up to relabeling/orientation.

## Verdict

The finite cyclic register is Pareto-optimal but NOT globally unique.

To promote it to a sourced carrier FIN needs the stronger notion:

    persistent record identity + content-independent routing.

This is a sharper source requirement than merely:
- minimum entropy;
- maximal freshness;
- zero symbol rewrite.
