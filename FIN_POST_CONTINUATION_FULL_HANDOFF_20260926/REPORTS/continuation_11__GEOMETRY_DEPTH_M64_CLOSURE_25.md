# GEOMETRY-DEPTH-M64-CLOSURE-25 — complete depth decomposition at the 8-phase M=64 peak

Status: **EXACT VARIANCE RECURSION + VALIDATED FFT PRODUCER; FINITE-SIZE NUMERICAL RESULT**.

This closes the explicit open item left in GEOMETRY-DEPTH-SUSCEPTIBILITY-22: the internal 13.38% of the M=64 susceptibility is now resolved by hierarchy depth.

At `alpha=0.8245`, very near the independently localized M=64 maximum,

    C_H = 300.60981480883402.

The exact nonnegative depth decomposition is

| depth | C_H contribution | fraction |
|---:|---:|---:|
| 0 | 260.38871755973042 | 0.86620165 |
| 1 | 25.451273250807844 | 0.08466548 |
| 2 | 9.535376490061706 | 0.03172011 |
| 3 | 3.755732468354507 | 0.01249371 |
| 4 | 1.478715039879552 | 0.00491905 |

The depth sum closes to the total `C_H` at numerical precision.

For comparison, the root fractions at the finite-size peaks are

    M=16: 0.65814111
    M=32: 0.77974095
    M=64: 0.86620165.

The deficit from the root therefore decreases

    0.34186 -> 0.22026 -> 0.13380.

The corresponding effective two-size decay exponents of the deficit are only numerical scouts,

    p_16->32 ~= 0.634,
    p_32->64 ~= 0.719,

and are not promoted to an asymptotic exponent.

Inside the M=64 non-root tail, the successive ratios are

    depth2/depth1 ~= 0.37465,
    depth3/depth2 ~= 0.39387,
    depth4/depth3 ~= 0.39372.

This near-geometric tail is striking but currently only a single-size finite-depth observation.

The result strengthens the interpretation that the rapidly growing M=64 susceptibility is increasingly controlled by reorganization of the coarsest split, with a smaller multiscale tail rather than a single isolated two-state switch.

No thermodynamic-limit theorem, first-order transition theorem, physical geometry, `D_H=2`, SI scale or laboratory claim follows.
