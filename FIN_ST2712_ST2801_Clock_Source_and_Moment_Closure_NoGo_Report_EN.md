# FIN ST2712–ST2801 — Clock Source and Moment-Closure No-Go

## Scale theorem

For every configuration generator \(Q\) and \(c>0\), the rescaled generator

\[
Q_c=cQ
\]

has the same stationary state, detailed-balance ratios and embedded jump chain.
Every timestamp-free sequence of visited states has the same law. Only holding
times change:

\[
T_x\sim\operatorname{Exp}(c\lambda_x).
\]

Therefore the positive scale of a generator is operationally invisible without
a calibrated temporal record. Generators with the same jump chain form an
\(\mathbb R_+\) clock torsor.

## Moment-closure obstruction

For the ternary energy, flipping spin \(i\) gives

\[
\Delta E_i
=2\theta x_i\sum_{f\ni i}\tau_f\prod_{j\in f\setminus i}x_j.
\]

Hence

\[
Qx_i=-2x_ir_i(x)

\]

contains higher Walsh characters. First moments couple to pair, triangle and
higher moments and do not close exactly.

At \(\theta=0\), heat-bath rates give

\[
Qx_i=-x_i
\]

for every degree-one mode. This identity decay cannot reproduce the nontrivial
strict spectrum of \(A\). Configuration dynamics therefore does not project
canonically to FIN vertex heat dynamics.

## Spectral-clock ambiguity

Several equally natural internal conventions exist:

- set the clock by the strict spectral gap,
- set it by \(\lambda_{\max}(A)\),
- use trace density,
- normalize the maximum transition rate.

They disagree. Each provides a dimensionless convention, not seconds. Unitary,
heat and wave channels also require different refinement scalings.

## Time-arrow obstruction

An oriented configuration circulation is a valid carrier of a kinetic arrow.
Rotation averaging can preserve it, but reflection reverses it. Full
\(D_{12}\) averaging cancels every reflection-odd current. The pair \(J,-J\)
has the same stationary state and entropy-production magnitude.

Thus neither equilibrium information nor entropy production selects an arrow
orientation. A boundary drive, chiral preparation or nonequilibrium resource
would supply it conditionally.

## Minimal clock record

For normalized exposure \(u_r=\lambda_{x_r}T_r\), one has

\[
u_r\sim\operatorname{Exp}(c).

The likelihood and estimator are

\[
\log L=n\log c-c\sum_ru_r,
\qquad
\widehat c=\frac{n}{\sum_ru_r}.

Fisher information is

\[
I_n(c)=\frac{n}{c^2}.

Consequently the clock rate is statistically identifiable after calibrated
timestamps are supplied, but not before. One hundred waiting times give a
conditional relative one-sigma lower scale of about ten percent.

## Gate

| Requirement | Result |
|---|---|
| rate-scale torsor theorem | PASS |
| moment nonclosure theorem | PASS |
| clock Fisher theorem | PASS |
| strict spectral normalization | FAIL |
| strict timestamp calibration | FAIL |
| strict circulation polarity | FAIL |
| exact \(A\to Q\) moment intertwiner | FAIL |
| physical seconds | FAIL |
| OA custody | FAIL |
| laboratory evidence | FAIL |

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

## Final interpretation

FIN contains dimensionless temporal parameters, mode frequencies and valid
time-arrow carriers. It does not contain a unique kinetic normalization,
physical clock calibration or selected circulation polarity. A soliton
frequency can parameterize internal phase, but it does not become seconds
without a calibration map.

The only results capable of changing this verdict are an independently sourced
scale-charged clock object or an operational calibration event with timestamps.

No selector, SI time, physical arrow, apparatus evidence, Standard Model,
gravity, \(L_{\rm total}\), or Theory of Everything closure follows.
