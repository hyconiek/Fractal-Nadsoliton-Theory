# Release 10.64 — FIN ST118–ST129: Operational Equivalence, Certified Design, and Typed Tests

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint and reproducibility package
- **Version:** 1.0.0
- **Publication date:** 2026-08-11
- **Language:** English
- **License:** CC BY 4.0
- **Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Abstract

Release 10.64 executes FIN research programs ST118–ST129. It identifies the exact additional operational equivalence needed to obtain the fine-blind observable algebra, classifies the complete circulant section of the hidden lift spectrahedron, independently replays the transcendental-kernel fold certificate, classifies every maximal code recoverable after coarse/fine pinching, and tests selector states, active feedback, robust continuous design, mixed refinements, correlated detector noise, typed thermal dilations, and minimax fine observation.

The principal structural theorem is conditional but sharp. If arbitrary fine-sector basis changes are declared operational gauge, the fixed-point algebra is uniquely

\[
M_{12}\oplus\mathbb C,
\]

and Haar twirling produces its trace-preserving expectation. The theorem isolates the missing premise rather than deriving it: strict (A), swap symmetry, and ordinary locality do not supply the fine gauge.

The real symmetric circulant section of the complete 78-dimensional lift fiber is exactly

\[
[s,\infty)\times\prod_{d=1}^{6}[-w_d,w_d].
\]

It therefore has (64) vertices and one recession extreme ray. This is a complete exact classification of the seven-parameter cyclic section, not of all noncirculant faces.

ST120 independently reproduces the ST108 Krawczyk inclusion using only NumPy and mpmath as third-party dependencies. It imports neither SciPy nor FIN research modules. The accepted radius-(10^{-8}) box again has margin (9.999987060638205\times10^{-9}). The replay corrects an earlier precision description: the 70-digit internal transcendental evaluation is converted to outward binary64 intervals of approximately one ulp, (10^{-16})–(10^{-18}), before the matrix proof.

ST124 proves a unique stationary Chernoff design inside a radius-(10^{-6}) box centered at

\[
(\beta_*,s_*)=(2.186189925005268,0.539833703306343),
\]

but does not prove global uniqueness on (eta\in[0.05,5]). The attempted interval derivative cover leaves 274 boxes unresolved; a 2,001-point floating audit finds one sign change and is reported only as strong evidence.

## Program outcomes

- **ST118 — Proven conditionally:** arbitrary fine-sector unitary gauge equivalence selects (M_{12}\oplus\mathbb C); the equivalence itself remains an axiom.
- **ST119 — Proven:** the circulant lift section is a six-cube times one ray, with 64 vertices.
- **ST120 — Proven in source-code interval scope:** independent two-dependency replay of the transcendental fold certificate.
- **ST121 — Proven:** complete normal form of all maximal twelve-dimensional pinching-recoverable codes.
- **ST122 — Proven:** strict functional-calculus states remain uniform for covariant branch instruments and cannot close QW-2191.
- **ST123 — Proven for a constructed model:** a scalar active branch requires an explicit pump crossing the loss threshold and a saturation term.
- **ST124 — Proven locally; global status strong evidence:** strict stationary-root inclusion, incomplete global cover.
- **ST125 — Proven:** odd refinement degrees preserve (mathbb Z_2) holonomy; the first even degree erases antiperiodicity.
- **ST126 — Proven conditionally:** a gauge-fixed local net is consistent, but the gauge—not locality—does the selecting work.
- **ST127 — Proven:** arbitrary temporal correlation can eliminate every positive finite-count error exponent while preserving all one-event marginals.
- **ST128 — Proven for a restricted class:** a reset-to-Gibbs channel lies outside the declared (H_B=0) bath subclass; standard thermal-operation separation remains open.
- **ST129 — Proven:** (P_f/12) is the unique diagonal minimax detector against the normalized adversarial lift simplex, with value (1/12).

## Central conclusion

FIN now has a mathematically precise conditional definition of its projection layer: it is the fixed algebra of an explicit fine-sector gauge equivalence. This is a substantive clarification, but not strict closure. The theory still lacks a theorem generating that gauge equivalence, a non-premise selector, dimensional units, calibrated preparation and apparatus, and independent experimental records.

No legacy-to-strict completion or physical-role transfer, Standard Model, gravity, (L_{\rm total}), or Theory-of-Everything closure is claimed.

## Included artifacts

- `FIN_ST118_ST129_Operational_Equivalence_Certified_Design_and_Typed_Tests_Report_EN.pdf`
- `FIN_ST118_ST129_Operational_Equivalence_Certified_Design_and_Typed_Tests_Report_EN.tex`
- `fin_st118_st129_research.py`
- `fin_st120_minimal_replay.py`
- `test_fin_st118_st129.py`
- `FIN_ST118_ST129_Results.json`
- `FIN_ST118_ST129_Summary.csv`
- `FIN_ST119_Circulant_Lift_Polyhedron.json`
- `FIN_ST120_Minimal_Transcendental_Replay_Certificate.json`
- `FIN_ST124_Continuous_Beta_Interval_Certificate.json`
- `FIN_ST129_Minimax_Fine_Observable.json`
- `FIN_ST118_ST129_Figures/` — eight figures
- `FIN_ST118_ST129_SHA256SUMS.txt`

## Validation

The complete suite passes **28/28 tests**. Tests include live reconstruction of positive semidefinite circulant vertices, Knill–Laflamme graph-code identities, strict-state uniformity, pump threshold and cost balance, local Krawczyk inclusion, refinement parity, correlated-error constancy, and the minimax detector.

## Recommended continuation

ST130–ST141 are specified in the report. Highest priority is ST130: seek a strict operational source for the fine (U(12)) gauge equivalence without naming the fine sector in advance, or prove a no-go in the strict functional-calculus class. ST132 should remove the imported numerical center from the minimal fold replay. ST139 should replace the arbitrary-correlation impossibility by positive finite-count bounds under a declared Markov or (alpha)-mixing condition.

## Keywords

FIN; operational equivalence; gauge-fixed algebra; conditional expectation; spectrahedron; circulant polyhedron; Krawczyk operator; quantum error correction; active bifurcation; Chernoff information; temporal correlation; thermal operations; minimax observation.
