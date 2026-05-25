# P2027 S977 p1950-next Rank-3 / GB-Null Adaptive Quadrature Witness

Status: `OPEN_OBSTRUCTION_WITH_TRACE`

Local verdict: `PASS_RANK3_QUOTIENT_IDENTIFIABLE_ON_SCALAR_B1_WITH_GB_NULL_TRACE`

This packet consumes the direct P1848 B1 operator profiles and keeps the strict-only
kernel lane.  The four scalar profiles obey the explicit relation

`R2 - 4*Ric2 + Riem2 - GB = 0`.

Therefore the full B1 scalar Gram matrix has rank `3` and nullity `1`.
The exact null vector in `(R2,Ric2,Riem2,GB)` coordinates is
`(1,-4,1,-1)`.

## Rank-3 Quotient

Independent representative channels: `R2, Ric2, Riem2`.

Rank-3 quotient coefficients:
- a_R2=9.9999999999999922e-01
- a_Ric2=1.1656308464946203e-15
- a_Riem2=6.0663244882429037e-17
- a_GB_representative=0.0000000000000000e+00

Full-row residual for this GB-zero representative:
`1.110e-16` in L-infinity norm.

The equivalent four-channel family is

`a_full = (1,0,0,0) + tau*(1,-4,1,-1)`.

The minimum-norm family point has `tau=-5.2631578947368418e-02`.

## Adaptive Quadrature

P1950 used a Simpson refinement replay and exposed a large discretization
bound.  P2027 replaces that gate with adaptive `scipy.integrate.quad`.
Any matrix/RHS entry involving `Riem2` or `GB` uses the left-endpoint transform

`x = epsilon + (1-epsilon)*t**power`

to resolve the `Kdd**2` endpoint profile.

Primary power: `5`.
Check power: `7`.
Relative primary/check gap: `1.233e-16`.
Gate tolerance: `1.000e-08`.
Adaptive gate pass: `True`.

## Honest Interpretation

The scalar B1 quotient is locally identifiable after removing the GB null
direction.  This does not identify an independent `a_GB`, does not prove
four-channel counterterm uniqueness, and does not replace the need for a
tensor-resolved projection if the theorem target requires the full curvature
operator basis.
