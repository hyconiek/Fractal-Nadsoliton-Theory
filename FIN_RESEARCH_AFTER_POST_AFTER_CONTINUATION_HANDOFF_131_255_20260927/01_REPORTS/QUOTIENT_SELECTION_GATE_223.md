# QUOTIENT-SELECTION-GATE-223
## The first effective level is not unique until the physical role of the coarse variable is specified

Date: 2026-09-26

Reports 220-222 correct the previous one-quotient narrative.

The exact internal algebra supplies two natural quotient coordinates:

    Z12 -> Z3
    Z12 -> Z4.

Both:
- are symmetry-defined;
- are exact lumpings of the 12-state effective chain;
- admit accurate direct microscopic Mori-Zwanzig reductions.

They differ in what they optimize.

## Z3

Advantages:
- lower basin escape/conductance;
- groups the lowest-barrier d3 moves internally;
- longer sector residence time;
- matches the earlier three-sector metastability construction.

## Z4

Advantages:
- slightly smaller global quotient spectral gap over N=3..8;
- therefore carries the slowest quotient relaxation eigenmode among these two
  factors.

## Scientific gate

A coarse variable should be selected by a preregistered physical role:

### metastable unit / identity carrier

Use:
- capacity;
- committor;
- residence time;
- internal-vs-external mixing.

This presently favors Z3.

### slow collective coordinate

Use:
- spectral gap / eigenmode separation.

This presently favors Z4 slightly.

### geometry coordinate

Requires an additional operational incidence/distance bridge.
Neither quotient wins automatically.

## Potential asymptotic crossover

In the idealized d3+d4-only product limit:

    gap_Z4 ~ 2 q3,
    gap_Z3 ~ 3 q4.

Then Z3 becomes slower only when

    q3/q4 > 3/2.

Using the bare saddle difference with unit prefactor would place this around

    N≈23.3.

A descriptive N=3..8 rate-ratio extrapolation instead gives roughly

    N≈47.6.

Neither is a theorem.

Therefore no large-N hierarchy selection should be claimed until larger-N or
a controlled capacity/Eyring-Kramers result is available.
