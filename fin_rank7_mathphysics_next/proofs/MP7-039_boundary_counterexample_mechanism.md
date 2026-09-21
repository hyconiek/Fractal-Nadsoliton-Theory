# MP7-039 — symmetry-resolved transverse mechanism on the two-harmonic family

Scientific state: **PROVED_INTERVAL_ASSISTED LOCAL TRANSVERSE CROSSING**.

Use the exact invariant period-four, reflection-even family

`h_j = J cos(pi j/2) + K (-1)^j`

with

`x=sinh(J)/(cosh(J)+exp(-2K))`,
`y=(cosh(J)-exp(-2K))/(cosh(J)+exp(-2K))`,

and stationary equations

`J=g lambda3 x/6`, `K=g lambda6 y/12`.

The full X7 Hessian contains two symmetry-related `(4,5)` blocks.  Their common
determinant is

`Delta45=(1/g-lambda4/12)(1/g-lambda5/12)-lambda4 lambda5 x^2/144`.

A three-variable parametric interval Krawczyk calculation for the two stationary
equations together with `Delta45=0`, uniform over the accepted outward spectral
intervals, proves a unique root in

```
J in [0.038434553978898216, 0.03843495397889821],
K in [0.18527186140188867,  0.18527226140188868],
g in [5.171841631942819,     5.171842031942818].
```

The Krawczyk image is strictly interior and refines the gain to approximately

`g*=5.17184183194`.

Implicit differentiation of the stationary branch, performed with interval
arithmetic on the same box, gives

`d Delta45/dg in [0.0016957351021672582, 0.0016960387698923191]`.

Thus the crossing is simple along the two-harmonic stationary branch.

At `Delta45=0`, one eigenvalue vanishes in each of the real `(4cos,5cos)` and
`(4sin,5sin)` blocks.  The critical space is therefore two-dimensional.  The
base period-four reflection-even field is fixed by translation `T^4` and a
reflection, so its isotropy subgroup is `D3` of order six.  On the critical
2D space, `T^4` acts as rotation by `2*pi/3`, while reflection acts as complex
conjugation.  Hence, for a complex critical coordinate `z`, the cubic invariant

`Re(z^3)`

is symmetry-allowed.

Consequently this event must not be described as a generic one-dimensional
pitchfork.  The certified result is the symmetry-resolved transverse crossing;
a nonlinear branch classification would require the actual cubic/quartic
normal-form coefficients and is not asserted here.

Replay: `scripts/mp7_039_transverse_crossing.py`.
Machine record: `results/MP7-039_transverse_crossing.json`.
