# K6-TRUNK-SUPPORT-LADDER-103
## One connected symmetry trunk organizes asymptotic support sizes 2,3,4,5,6

Date: 2026-09-26

This is a synthesis of reports 92, 100-102.

The pure-k6 stationary trunk has large-g support

    size 6:
      {0,2,4,6,8,10},
      index 5,
      stabilizer order 12.

At the exact k4 crossing on this trunk, the exact k4×k6 factorization produces
two nonzero global sheets.

## Positive-large k4 sheet

Large-g support:

    size 2:
      {0,6},
      index 1,
      stabilizer order 4.

This is the d=6 pair branch.

It later undergoes the transverse event at

    g≈5.538559619700

and emits the daughter

    size 3:
      {0,3,6},
      index 2,
      stabilizer order 2.

## Negative k4 sheet

Large-g support:

    size 4:
      {2,4,8,10},
      index 3,
      stabilizer order 4.

It later undergoes the transverse event at

    g≈5.764771557084

and emits the daughter

    size 5:
      {2,4,6,8,10},
      index 4,
      stabilizer order 2.

## Result

Within one symmetry-organized neighborhood of the pure-k6 trunk, FIN therefore
realizes every asymptotic support size

    boxed:
    2,3,4,5,6

with the exact Morse sequence

    boxed:
    1,2,3,4,5.

This is the clearest finite-g realization so far of the large-g theorem

    index = support size - 1.

The structure is not a single linear path 2->3->4->5->6.

It is a branching tree:

                         S6 / idx5
                            |
                     exact k4 crossing
                      /             \
                     /               \
              S2 / idx1          S4 / idx3
                  |                  |
        transverse breaking   transverse breaking
                  |                  |
              S3 / idx2          S5 / idx4

So one parent symmetry sector contains two complementary support-complexity
ladders:
- 2 -> 3;
- 4 -> 5;

while the trunk itself supplies 6.

This provides strong evidence that support cardinality is not merely an
asymptotic bookkeeping device: finite-g representation-theoretic events
actively organize transitions between support strata.
