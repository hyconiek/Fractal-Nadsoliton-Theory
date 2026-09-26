# METASTABLE-C3xC4-SKELETON-205
## The observed low-barrier localized-state graph aligns exactly with the CRT product coordinates

Date: 2026-09-26

Status:
exact graph identity using the previously identified d=3 and d=4 localized
transition families.

The 12 localized minima are labelled by j in Z12.

Under

    j -> (j mod3, j mod4),

the two important transition distances act as:

### d=3

    (a,b)
      ->
    (a,b±1)

because

    3 mod3=0,
    3 mod4=-1.

So d=3 moves only in the Z4 coordinate.

### d=4

    (a,b)
      ->
    (a±1,b)

because

    4 mod3=1,
    4 mod4=0.

So d=4 moves only in the Z3 coordinate.

## 1. Product graph

Keep all ±d3 and ±d4 transitions.

Every localized minimum then has:
- two Z4 neighbors;
- two Z3 neighbors.

After CRT relabelling the adjacency is EXACTLY

    boxed:
    C3 square C4

(the Cartesian graph product).

The matrix residual of the graph identification is zero.

Its degree is four.

Its unweighted Laplacian spectrum is

    {0,
     2,2,
     3,3,
     4,
     5,5,5,5,
     7,7}.

## 2. Relation to the successful 3-state coarse variable

The previously selected metastable sectors are exactly

    C_a
      =
      {j : j mod3=a}.

So the three-state coarse variable is precisely the projection

    Z3 × Z4
      ->
    Z3.

Its eliminated internal coordinate is Z4.

This is not an arbitrary clustering discovered after the fact.

It coincides with an exact quotient of the original internal cyclic carrier.

## 3. Boundary

The 12 vertices are alternative localized states of one statistical FIN unit.

Therefore C3×C4 is currently a STATE-SPACE transition geometry.

It is not twelve simultaneous spatial sites.
