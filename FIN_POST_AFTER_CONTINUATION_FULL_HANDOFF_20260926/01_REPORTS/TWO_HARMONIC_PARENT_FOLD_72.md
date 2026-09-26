# TWO-HARMONIC-PARENT-FOLD-72
## The small MP7-039 parent and the large R7P-036 counterexample are two sheets of one two-harmonic fold

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact two-harmonic stationary reduction inherited from R7P-034;
- new fold location is numerical;
- full-H7 index change is numerical with large spectral margins;
- no interval certificate yet.

## 1. Exact two-harmonic family

Use

    h_j
      = J cos(pi j/2)
        +K (-1)^j,

with exact stationary equations

    J = g lambda3 x/6,
    K = g lambda6 y/12,

where

    x = sinh(J)/(cosh(J)+exp(-2K)),

    y = [cosh(J)-exp(-2K)]/
        [cosh(J)+exp(-2K)].

This is the exact invariant family already proved in R7P-034.

## 2. New fold

Solve the two stationary equations together with singularity of their
two-by-two Jacobian.

The numerical root is

    boxed:
    J_fold ≈ 0.718505246005630
    K_fold ≈ 0.490908854614046
    g_fold ≈ 4.621196599489125.

In full X7 coordinates this corresponds approximately to

    theta3c ≈ 1.256671377194,
    theta6  ≈ 1.111171682901,

with all other coordinates zero.

## 3. Full Hessian at the fold

The full H7 eigenvalues are approximately

    -0.06024604,
    -0.06024604,
     0,
     0.11696230,
     0.11818645,
     0.11818645,
     0.14190333.

Thus:
- the fold direction is one-dimensional;
- two transverse negative directions already exist;
- all remaining directions are well separated and positive.

Hence the two sheets joined by the fold have full Morse indices

    2 and 3.

## 4. Identification of the two sheets at g=5

The large-amplitude sheet is the previously certified R7P-036 stationary
counterexample, represented in the g=5 full atlas by orbit `i=3`:

    index = 2,
    stabilizer size = 6,
    orbit size = 4.

The small-amplitude sheet is atlas orbit `i=10`:

    index = 3,
    stabilizer size = 6,
    orbit size = 4.

The small sheet is precisely the D3-symmetric parent that later reaches the
MP7-039 transverse crossing near

    g_D3 ≈ 5.171841831943.

Therefore the structural chain is

    large two-harmonic index-2 parent
        |
    fold g≈4.62119659949
        |
    small two-harmonic index-3 parent
        |
    D3 crossing g≈5.17184183194.

## 5. Consequence for the old counterexample

R7P-036 was previously a certified isolated witness showing that unrestricted
stationary index<=1 is false.

The new continuation result gives that witness a branch role:

> it is one sheet of the same two-harmonic component whose other sheet later
> undergoes the D3 symmetry-resolved crossing.

So the counterexample is not an accidental remote stationary point. It is part
of the symmetry-breaking branch network.

## 6. Updated branch graph fragment

The branch graph now contains

    R7P-036 large parent
        index 2
          |
          | two-harmonic fold
          | g≈4.62119660
          v
    small D3 parent
        index 3
          |
          | MP7-039
          | g≈5.17184183
          v
    index-1 D3 parent above crossing
        + two triplets of index-2 daughters.

One daughter triplet then connects, through report 69, to the main
localization saddle/minimum component.

## 7. Proof boundary

The fold location and branch join are currently numerical.

The clean next certification target is only three variables `(J,K,g)`, so an
interval Krawczyk proof should be inexpensive compared with a full seven-
dimensional branch certificate.

### NEXT: TWO-HARMONIC-FOLD-CERTIFICATE-73

Certify:
1. a unique fold root near the values above;
2. nonzero fold transversality;
3. the two negative transverse H45 directions;
4. local index-2/index-3 branch exchange.

That would turn this new graph edge from numerical to interval-assisted.
