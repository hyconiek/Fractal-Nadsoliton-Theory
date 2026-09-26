# D6-FACTORIZED-GENEALOGY-92
## Exact k4×k6 factorization plus a transverse daughter gives a second 2→3→4 support network

Date: 2026-09-26

Status:
- k4/k6 factorization and the k6 pitchfork location are exact;
- the transverse crossing location is high-precision numerical;
- the large-g daughter support is solved exactly.

For a field containing only k=4 cosine and k=6 parity,

    h_j=J4 cos(2πj/3)+J6(-1)^j,

the uniform twelve-label average factorizes because the mod-3 and parity
coordinates are independent.

Hence the stationary equations decouple exactly:

    J4 = g(lambda4/6) R4(J4),

    J6 = g(lambda6/12) tanh(J6),

where, with `t=3J4/2`,

    R4=(exp(t)-1)/(exp(t)+2).

The k6 factor therefore has an exact pitchfork at

    boxed:
    g6=12/lambda6
      =5.123427551398616.

On the large positive k4 sheet, the k4 dual coordinate at this event is

    s4=1.775824710981854.

The pure-k4 parent has two negative directions plus the k6 zero mode at the
crossing.  For g>g6 the nonzero ±k6 daughters exist and have index 2 locally.

These mixed k4+k6 daughters are precisely the finite-g ancestors of the
distance-6 pair support

    S2={0,6},

whose large-g data are

    p*=(1/2,1/2),
    mu=0.561776644984393,
    gap=0.390363673524383,
    stabilizer order 4,
    orbit size 6.

Farther up the same pair branch a second, transverse zero occurs at

    boxed:
    g_T≈5.538559619699902.

The crossing eigenvector lies in the k3-sine/k5-sine sector and breaks the
order-4 pair stabilizer down to an order-2 reflection stabilizer.

A nonzero daughter exists on the high-g side and has index 2.

Its large-g destination is the strict three-label support

    boxed:
    S3={0,3,6}

with exact weights

    (0.369069585031288,
     0.261860829937423,
     0.369069585031289),

common field

    mu=0.459555689457877,

and outside gap

    Delta=0.371843524942702.

Its stabilizer order is 2, so its orbit
size is 12.

The pure-k4 parent itself tends at large g to

    S4={0,3,6,9},
    p*=1/4 on each label,

with common field

    mu=0.366594808222202

and asymptotic index 3.

Thus the same symmetry block contains three asymptotic complexities:

    S2 pair       -> index 1,
    S3 daughter   -> index 2,
    S4 pure-k4    -> index 3.

This is a second independent realization of the support-size/Morse-index
hierarchy.
