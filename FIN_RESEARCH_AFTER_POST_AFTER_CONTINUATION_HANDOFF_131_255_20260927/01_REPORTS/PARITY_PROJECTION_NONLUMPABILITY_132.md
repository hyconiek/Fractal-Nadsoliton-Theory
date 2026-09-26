# PARITY-PROJECTION-NONLUMPABILITY-132
## Exact leave-one-out Gibbs does not close Markovianly on parity count

Date: 2026-09-26

Status:
EXACT FINITE-STATE COUNTEREXAMPLE to strong lumpability.

Microscopic process:
exact leave-one-out Gibbs heat-bath.

Resolved observable:

    Y(n) = number of copies on even labels.

A necessary condition for strong lumpability is that, for all microstates n
with the same Y, the total transition rate into each target Y' block is the
same.

This condition fails.

## 1. Explicit N=4, Y=2 counterexample

At

    g = 5.145228719489142,

different count vectors with the same

    N=4,
    Y=2

have different total rates Y=2 -> Y=3.

Across the complete Y=2 fiber:

    min rate ≈ 1.104539312219
    max rate ≈ 1.533899342869

so

    boxed:
    max/min ≈ 1.388723177.

The same spread occurs for the reverse Y=2 -> Y=1 rate by parity symmetry.

Therefore Y(t) is not an autonomous continuous-time Markov chain.

## 2. The defect does not disappear trivially with small-N growth

For the central parity fiber:

    N=3:
      one direction has ratio 1.0000,
      the other 1.1845;

    N=4:
      ratio ≈ 1.3887;

    N=5:
      ratios ≈ 1.8974 and 1.3598;

    N=6:
      ratio ≈ 1.7275.

There is no evidence in N=3..6 that the hidden rate dependence simply vanishes.

## 3. Why the drift can still look closed

Two microstates can have the same parity drift while having different:
- total activity;
- waiting-time law;
- second-step transition statistics.

Therefore agreement of a first projected drift or a relation of the form

    P L J = L_coarse

does not imply equality of projected processes.

The defect appears at second order and in multi-time statistics.

## 4. Consequence for task 132

Exact autonomous closure must use a stronger condition such as exact
intertwining/lumpability.

If that fails, the correct reduced object is either:
- an approximate Markov process with a proved error bound; or
- a non-Markov process with an explicit memory kernel.

Report 133 tests the second option directly.
