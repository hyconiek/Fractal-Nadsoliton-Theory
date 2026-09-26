# Z2-SUPPORT-SELECTION-105
## In two independent FIN branches, transverse Z2 breaking selects one label from a symmetry-related pair

Date: 2026-09-26

Status:
- exact support/stabilizer statement for the two branch families already
  continued in reports 92 and 102;
- this is not claimed as a universal theorem for every FIN Z2 event.

## 1. First example: d=6 pair -> three-label daughter

Parent large-g support:

    S2={0,6}.

Its stabilizer is the order-4 subgroup

    H2={1,R,T6,T6 R}.

The candidate labels

    {3,9}

form a two-element orbit under H2.

The transverse critical mode at

    g≈5.538559619700

is odd under the symmetry exchanging 3 and 9.

The two broken-symmetry daughters therefore choose one of them:

    {0,3,6}
or
    {0,6,9}.

These two supports are D12-equivalent.

For the representative

    S3={0,3,6},

the stabilizer is

    {1,T6 R},

of order 2.

Hence:

    support size:      2 -> 3
    stabilizer order:  4 -> 2
    D12 orbit size:    6 -> 12
    Morse index:       1 -> 2.

## 2. Second example: four-label parent -> five-label daughter

Parent support:

    S4={2,4,8,10}.

Its stabilizer is again

    H4={1,R,T6,T6 R}.

The labels

    {0,6}

form a two-element orbit under this parent symmetry.

At

    g≈5.764771557084

the transverse mode is even under R and odd under T6.

The two daughters therefore select one member of {0,6}.

A representative is

    S5={2,4,6,8,10}.

Its stabilizer is

    {1,R},

again of order 2.

So once more:

    support size:      4 -> 5
    stabilizer order:  4 -> 2
    D12 orbit size:    6 -> 12
    Morse index:       3 -> 4.

## 3. Common combinatorial mechanism

The two independent examples realize the same support-level rule:

    parent support S
      + residual Z2 exchanging labels a<->b
      + odd transverse mode

        ->

    daughter support S union {a}
or
    daughter support S union {b}.

Breaking the residual Z2:
- selects one member of a two-label symmetry orbit;
- halves the stabilizer;
- doubles the full D12 orbit;
- increases support cardinality by one;
- increases asymptotic Morse index by one.

## 4. Why this matters for the emergence picture

This gives a concrete meaning to “symmetry breaking creates a new state” inside
the zero-temperature support language.

The new daughter is not just a rotated continuous perturbation.

At large g it resolves a previously symmetry-equivalent pair of candidate
labels and retains one of them in the asymptotic support.

So in these two examples the finite-g odd mode has a direct combinatorial
large-g interpretation:

    odd order parameter
      <-> binary label selection.

## 5. Boundary

This mechanism is exact for the two branch families analyzed here.

It should not yet be generalized to:
- D3 cubic branching;
- D6 two-dimensional branching;
- k5 degree-12 angular selection;
- folds without symmetry breaking.

Those events have different invariant grammars.
