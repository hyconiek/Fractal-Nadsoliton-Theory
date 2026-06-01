# P2378 S1328: unit-normalized transport coupling insufficiency theorem

Status: `OPEN_PROGRESS_UNIT_NORMALIZED_TRANSPORT_INSUFFICIENT_SOURCE_STRENGTH_STILL_OPEN`

## Result

P2378/S1328 proves that the exact P2377 transport primitive still lacks enough source strength if its total transport mass is unit-normalized.
For `K_strict(d)+M*C(d)`, the d5 chamber requires `M>(3*K1-K5)/(C5-3*C1)`.

## Certificate

- Numerator `3*K1-K5`: `1.3858257945714303`.
- Denominator min: `0.7517322065151169` at `{'eta': 1.8, 'beta_tors': 0.1}`.
- Denominator max: `1.1786549963416464` at `{'eta': 2.0, 'beta_tors': 0.0}`.
- Threshold range: `{'minimum_tau_gt': 1.1757688202848233, 'maximum_tau_gt': 1.8435099395246706}`.
- Unit mass insufficient on rectangle: `True`.
- Unit grid scans select: `[{'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}, {'3,3': 24}]`.
- Just-above-threshold grid scans select: `[{'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}, {'0,4': 12}]`.

## Hard limits

- This is an insufficiency theorem for normalized transport strength, not a source theorem fixing a super-unit coupling.
- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.
