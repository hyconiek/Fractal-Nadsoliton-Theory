# MP7-015 — exact g=37/10 aligned stationary problem

Scientific state: **DOMAIN PROVED; THREE CANDIDATES IDENTIFIED; TWO NONZERO ROOTS
LOCAL-INTERVAL-CERTIFIED; GLOBAL EXHAUSTION OPEN**.

Fix the exact gain

```
g = 37/10.
```

By MP7-011 any global minimizer of the full dual is D12-equivalent to an aligned
nonnegative C4 representative.  By MP7-013/014, in the entire Target-P window
`0<g<=250/67`, every nonzero aligned stationary point is interior.  Thus at
`g=37/10` the stationary equation is

```
s = g E_s[C4],       s_i>=0,
```

with either `s=0` or all four coordinates strictly positive.

## Compact root domain

Each normalized feature coordinate is bounded above by its coordinate maximum
`m_i`, hence every nonnegative fixed point lies in the exact order box

```
0 <= s_i <= g m_i.
```

This is an invariant box for the isotone fixed-point map `T(s)=g E_s[C4]`.
The lower fixed point is exactly zero; monotone iteration from the upper corner
converges to the greatest fixed point.

## Candidate navigation

A deterministic numerical search using 2048 Sobol starts plus structured and
boundary starts found exactly three clusters in this box:

1. `s=0`, `Phi=0`, local minimum;
2. an interior localized local minimum near
   `(1.79742293,1.89006397,1.88931419,1.34918293)`, with
   `Phi≈0.00823514169`;
3. an interior index-one saddle near
   `(0.96053196,1.02089340,0.98183683,0.70059551)`, with
   `Phi≈0.04879007510`.

The numerical search is navigation only.

## Local interval certificates

Using the accepted outward strict spectral intervals and exact `g=37/10`,
separate 4D Krawczyk tests strictly isolate both nonzero roots in radius-`1e-8`
boxes.  The localized root has positive C4 Hessian at its center and a positive
sine block; the saddle has one negative C4 eigenvalue at its center and a
positive sine block.  The strict Krawczyk inclusions are stored in
`results/MP7-015_g37_local_interval_roots.json`.

These local certificates do **not** prove that no additional interior roots
exist.  MP7-016 therefore remains a genuine global complement/objective task.
