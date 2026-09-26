# SECOND-D3-DAUGHTER-FOLD-70
## The opposite D3 daughter triplet forms a second signed-C4 fold and two known g=5 atlas orbits

Date: 2026-09-26

Status:
- numerical pseudo-arclength continuation;
- numerical fold root;
- exact D12 orbit matching of numerical roots to the existing atlas to
  machine precision.

## 1. Opposite daughter set

For

    g<g_D3,

the nonzero cubic D3 normal form selects

    cos(3phi)=-1.

The three local daughters are rotated by pi/3 relative to the g>g_D3
triplet.

A representative is continued away from the crossing.

## 2. Identification at g=5 before the new fold

On the first sheet, at g=5 the representative is

    theta ≈
    (0.271673875059, 0,
     0.083546704081, 0.144707136273,
     0.121216308472,-0.209952804979,
     0.396436330675),

with

    Phi ≈ 0.000697351440512
    index = 2.

A D12 translation maps it to

    (-0.271673875059,0,
     -0.167093408162,0,
      0.242432616944,0,
      0.396436330675),

which is atlas orbit `i=9` to about `7e-14`.

So report 66's local index-2 daughter is directly visible in the existing
g=5 stationary census.

## 3. Second fold

Continuing toward lower g reaches a turning point.

After a D12 shift into a reflection-fixed C4 chart, solving stationarity plus

    det(H_C4)=0

gives

    boxed:
    g_second ≈ 4.395526393543094

with signed C4 coordinates

    (s3,s4,s5,s6)
      ≈
      (-1.330174026970,
       -0.651860029116,
        0.720774033303,
        1.069731466438).

At this point the full Hessian has:
- one negative eigenvalue;
- one zero fold eigenvalue;
- five positive eigenvalues.

Across the fold the branch changes

    index 2 -> index 1.

## 4. Return branch at g=5

After the fold, the index-1 sheet returns to larger g.

At g=5 it gives a D12-equivalent representative of atlas orbit `i=2`:

    Phi ≈ -0.105513112212
    index = 1.

The D12 matching error is below `1e-13`.

Thus one connected component contains both existing atlas orbits:

    i=9 (index 2)
      -> fold at g≈4.39552639
      -> i=2 (index 1).

## 5. High-g continuation

The index-1 return sheet was numerically continued to at least

    g=10

without another detected index change.

At g≈10 its Hessian still has exactly one negative eigenvalue.

So within the explored range this second component remains a saddle family
rather than becoming a stable minimum.

This is a numerical nonconclusion beyond the explored interval; it is not an
asymptotic theorem.

## 6. Combined D3 picture

The D3 crossing therefore attaches two distinct local daughter triplets to two
different global branch components:

### + triplet
connects into:
    upper fold
      -> main index-1 localization saddle
      -> simple fold
      -> stable localized minimum.

### - triplet
connects into:
    signed-C4 index-2 branch
      -> fold at g≈4.39552639
      -> index-1 saddle branch.

This explains several previously separate entries of the g=5 stationary atlas
as parts of one symmetry-organized branch graph.
