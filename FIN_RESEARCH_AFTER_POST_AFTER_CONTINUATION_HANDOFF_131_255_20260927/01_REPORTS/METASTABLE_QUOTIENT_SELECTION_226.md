# METASTABLE-QUOTIENT-SELECTION-226
## Barrier connectivity selects Z3 as the metastable unit even though Z4 has a slightly smaller global relaxation gap

Date: 2026-09-26

Reports 220-224 resolve the apparent quotient ambiguity.

## Algebraic level

Both

    Z12 -> Z3
    Z12 -> Z4

are natural and exactly lumpable after the 12-state Markov reduction.

So algebra alone does not choose one.

## Spectral level

For N=3..8:

    gap_Z4 < gap_Z3

by a few percent.

Thus Z4 carries the slower global quotient relaxation mode.

## Metastable barrier level

The saddle filtration says something different and more directly relevant to
persistent basin identity:

At B3:
- the system splits into exactly three mod3 components.

At B4:
- those components merge.

No mod4 partition appears as a connected component at any intermediate
barrier threshold.

Indeed a mod4 class needs d4 edges to mix internally, but d3 edges of LOWER
barrier already carry it out of that class.

So mod4 is not a metastable basin partition in the mapped landscape.

## Capacity level

Consistently:
- raw microscopic cap/pi is smaller for mod3 than mod4 at N=6 and N=8;
- effective exit/conductance is smaller for mod3 for every N=3..8.

## Verdict

For the physical role

    "metastable unit that preserves identity before escape",

the current evidence selects

    boxed:
    Z3 sectors.

For the role

    "slowest global quotient eigenmode",

Z4 remains slightly slower in the tested finite-N range.

There is no contradiction.

The two criteria measure different properties.

This restores Z3 as the preferred metastable coarse unit while preserving the
important correction that it is not the unique algebraic or spectral quotient.
