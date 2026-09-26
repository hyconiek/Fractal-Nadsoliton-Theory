# DEGENERATE-LARGE-G-SUPPORTS-87
## Three non-isolated zero-temperature support classes and entropy selection

Date: 2026-09-26

Status:
- exhaustive singular-support classification for the supplied rank-seven mediator;
- entropy-selection mechanism is analytic;
- finite-g branch replays at g=100 are numerical checks.

The isolated atlas of report 84 excludes singular equal-field systems. Exhaustive
enumeration finds exactly three D12 classes that admit an interior positive
equal-field solution despite singularity.

## 1. Eight-label class

Canonical support:

    S8 = {0,1,3,4,6,7,9,10}.

D12 orbit size:

    3.

The equal-field manifold has dimension one.

Its affine feature rank is 6.

The leading-energy manifold contains the uniform-on-support point

    p_i*=1/8 on S8.

Because Shannon relative entropy is strictly convex, the O(1) entropy term
selects this equal-weight point uniquely along the flat leading-energy
manifold.

The outside field gap at the selected point is approximately

    0.274946106167.

At g=100 the full softmax stationary solve differs from this support law only
at about 1.1e-12 in L1, because the first algebraic 1/g correction vanishes
and the leading correction is exponentially small.

The full H7 Morse index is 6.

## 2. Nine-label class

Canonical support:

    S9 = {0,1,2,4,5,6,8,9,10}.

D12 orbit size:

    4.

The equal-field manifold has dimension two and affine feature rank 6.

Entropy selects the unique interior point

    p* approximately
    (0.1079280076,
     0.1174773181,
     0.1079280076,
     0.1079280076,
     0.1174773181,
     0.1079280076,
     0.1079280076,
     0.1174773181,
     0.1079280076).

The strict outside gap is approximately

    0.230420817792.

A full stationary solve at g=100 approaches this support law and has H7
Morse index 6.

Unlike S8, the entropy-selected weights are nonuniform, so there is a visible
algebraic 1/g correction.

## 3. Full twelve-label class

The full support has a four-dimensional equal-field degeneracy because A7 has
rank seven.

Entropy selects the exact uniform distribution

    p_i=1/12.

This is not merely asymptotic: the uniform state is an exact stationary branch
for every g.

At sufficiently large g all seven retained Hessian directions are negative, so

    Morse index = 7.

## 4. Why these are special

For an isolated support of size m, affine independence gives covariance rank
m-1.

For these singular supports, support size exceeds affine rank+1.

The leading zero-temperature interaction therefore leaves probability
directions that cost no O(g) energy.

The entropy term becomes the next selector:

    leading interaction:
        selects an equal-field affine manifold;

    entropy:
        selects one point inside that manifold;

    finite-g terms:
        perturb that selected point.

This is a distinct asymptotic mechanism from the 76 isolated branches.

## 5. Complete large-g support taxonomy

The current zero-temperature classification is therefore:

    isolated strict classes:
        76 D12 classes, support sizes 1..7;

    degenerate entropy-selected classes:
        one size-8 class,
        one size-9 class,
        one size-12 class.

No additional strict positive singular class was found in the exhaustive
4095-subset enumeration.
