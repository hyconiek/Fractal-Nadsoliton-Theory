# Release 10.98 — Frozen OA Dual-Dynamics Discrimination and Composite Robustness

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 10.98
- **Publication date:** 2026-08-31
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0

## Abstract

Release 10.98 reports six research rounds, ST1272–ST1361, freezing one
conditional Operational/State package for distinguishing the two FIN dynamics

\[
P_t=e^{-tA},\qquad U_t=e^{-itA}.
\]

Preparing vertex zero and measuring return gives

\[
C(t)=\frac1{12}\sum_k e^{-\lambda_k t},\qquad
Q(t)=\left|\frac1{12}\sum_k e^{-i\lambda_k t}\right|^2.
\]

A million-point scan of \(0\le t\le10\), supported by analytic derivatives,
finds the largest first separation near

\[
t_*\approx0.59453079,qquad |Q(t_*)-C(t_*)|\approx0.4112409998.
\]

The executable design freezes \(t=0.6\), where

\[
p_C=0.4123573972,qquad p_Q=0.8235756226.
\]

For ideal iid return counts, 29 clicks and the rule “choose quantum for at
least 19 returns” give both model-selection errors below one percent. Full
twelve-vertex readout increases the Chernoff information from approximately
0.100286 to 0.136480.

The principal falsification result is that a single time is not identifying
against an open-quantum alternative: energy-basis dephasing can exactly match
the classical return probability. A four-time design
\((0.3,0.6,1.2,2.0)\) remains numerically separated when clock scale is
constrained to \([0.9,1.1]\), but a widely free clock scale permits near
mimicry. Clock calibration and a frozen nuisance library are necessary.

The package includes a hashed JSON protocol, fail-closed validator, explicit
attempt/click/no-click accounting, synthetic fixtures, a one-shot decision
rule and a visualization. These are executable mathematical specifications,
not laboratory observations.

## Scientific boundary

The ideal 29-shot guarantee applies only to two frozen simple hypotheses with
perfect iid preparation and measurement. The four-time composite result is
strong numerical evidence, not an interval-certified minimax theorem. No
apparatus, raw natural events, independent custody, physical state-category
selection or Theory-of-Everything closure is claimed.

## Included files

- English PDF and LaTeX source
- six research scripts and shared OA functions
- frozen protocol and synthetic fixtures
- fail-closed protocol validator
- six result JSON files and summaries
- 90 evidence packets
- figure, combined tests, source bundle and SHA-256 manifest
