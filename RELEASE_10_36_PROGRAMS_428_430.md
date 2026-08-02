# Release 10.36 — FIN Research Programs P428–P430

## Type-correct cosine reduction, exact-rational Krawczyk isolation, and global dual feasibility

**Creator:** Żuchowski, Krzysztof  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Resource type:** Publication — Preprint  
**Version:** 1.0.0  
**Publication date:** 2026-08-01  
**Language:** English  
**Publisher:** Zenodo  
**Access:** Open access  
**License:** CC BY 4.0  
**Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Abstract

Release 10.36 reports a bounded local research batch, FIN Programs P428–P430. It closes the principal exact-mathematics frontier left by Release 10.35: the simultaneous seven-contact signed-moment candidate O144.

P428 repairs the formal typing of the cosine trust boundary. The previous abstract type `Rat -> Rat` cannot denote the ordinary real cosine. A Lean-compiled rational-cut interface now separates exact FIN-specific rational arithmetic from one standard real-analysis provider. All twelve phase-domain and alternating-sum conditions compile in Lean 4.28.0; the maximum rational Taylor width is approximately (1.913\times10^{-37}). No local Mathlib or equivalent real-trigonometric library is available, so the standard cosine provider remains an explicit local formalization boundary.

P429 establishes a strict exact-rational Krawczyk inclusion for the full 25-variable moment/contact system. A relative rational box radius of (10^{-20}) succeeds for all variables. The exact weighted infinity contraction bound is numerically at most (4.0673\times10^{-11}), the six interior nodes are strictly ordered, and the seven signed weights have pattern ((- + + + + - +)). The objective is enclosed in a rational interval of width (4\times10^{-20}) around (0.7073534677231137).

P430 proves global dual feasibility. Interval curvature certifies six stationary contact neighborhoods, interval monotonicity certifies the endpoint contact, and exact affine Bernstein ranges certify every complementary interval without recursive subdivision. Under the inherited signed-moment duality theorem, the certified primal–dual value is therefore globally optimal for the declared dimensionless moment problem. The theorem does not claim global uniqueness outside the P429 isolating contact box.

The release also integrates the research-backed reconstructed historical class

\[
K^*_{\mathrm{legacy}}(d)
=A\frac{\cos(\pi d/4+\pi/6)}{1+\beta d},\qquad A,\beta>0.
\]

The local coupling diagram contains four mechanism labels and two parameter cross-modulations, but no typed state self-map, feedback equation, or fixed-point theorem. Its displayed path exponent and integer-node derivations are algebraically refuted. Legacy* remains an intermediate historical kernel class; it is not silently identified with the strict kernel and transfers no physical roles.

## Main results

- **[Proven]** Type-correct rational-cut reduction of the twelve FIN cosine enclosures.
- **[Refuted]** `Rat -> Rat` as an adequate type for ordinary real cosine.
- **[Blocked locally]** Proof-kernel compilation of standard real cosine analysis because no suitable local library is present.
- **[Computer-assisted proof]** Exact-rational existence and local uniqueness of the seven-contact 25-variable KKT zero.
- **[Computer-assisted proof]** Global feasibility (-1\le p(x)\le0) on ([0,1]).
- **[Computer-assisted proof, inherited duality assumptions]** Global optimality of the certified signed-moment value.
- **[Proven audit]** Legacy* is a reconstructed fixed-phase subclass of the intermediate legacy kernel; the diagram does not define typed self-coupled dynamics.
- **[Refuted]** The diagram's (d^{1.6}d^{-0.6}\sim d^{-1}) claim and alleged zero sequence (2,5,8,11).

## New theoretical objects

- **O155 — Rational-Cut Cosine Provider Interface.** A type-correct boundary between rational formal arithmetic and standard real cosine.
- **O156 — Exact-Rational Seven-Contact KKT Isolating Box.** A strict Krawczyk inclusion and contraction certificate in (mathbb R^{25}).
- **O157 — Contact-Aware Global Dual Certificate.** A curvature/monotonicity/Bernstein proof of the full dual inequality.

## Scientific boundary

This release contains only local analytical, formal, and computational research. It contains no laboratory data, external audit, web research, external empirical record, or independent unblinding. It does not establish:

- a complete legacy-to-strict bridge;
- legacy physical-role transfer;
- QW-2191 discharge or a non-premise selector;
- a canonical orientation or phase origin;
- internally generated physical units;
- a physical state, preparation, clock, apparatus, instrument, environment, observer, or record;
- Standard Model, gravity, total-action, or Theory-of-Everything closure.

## Main files

- `FIN_Local_Research_Checkpoint_P428_P430_EN.pdf`
- `FIN_Local_Research_Checkpoint_P428_P430_EN.md`
- `FIN_Local_Research_Checkpoint_P428_P430_EN.tex`
- `FIN_Current_Local_Research_Report_EN.pdf`
- `FIN_Programs_428_430_Results.json`
- `FIN_Programs_428_430_Summary.csv`
- `FIN_Programs_428_430_P428_Cosine.csv`
- `FIN_Programs_428_430_P429_Krawczyk.csv`
- `FIN_Programs_428_430_P430_Global_Dual.csv`
- `FIN_Programs_428_430_LegacyStar_Coupling_Audit.csv`
- `FIN_Programs_428_430_P428_Cosine_Reduction.lean`
- `FIN_Programs_428_430_Figures/p428_p430_exact_closure_and_legacy_context.png`
- `fin_programs_428_430.py`
- `test_fin_programs_428_430.py`
- `FIN_PROGRAMS_428_430_RELEASE_MANIFEST.sha256`

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-mpl-cache python3 fin_programs_428_430.py
.elan/toolchains/leanprover--lean4---v4.28.0/bin/lean FIN_Programs_428_430_P428_Cosine_Reduction.lean
MPLCONFIGDIR=/tmp/fin-mpl-cache python3 -m unittest -v test_fin_programs_428_430.py
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P428_P430_EN.tex
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P428_P430_EN.tex
sha256sum -c FIN_PROGRAMS_428_430_RELEASE_MANIFEST.sha256
```

## Recommended next bounded batch

The selected next batch is P435, P436, and P440:

1. construct a genuine noisy-comb SDP primal/dual certificate, or document the exact local solver obstruction;
2. rigorously certify at least one positive erasure-aware O149 gain cell;
3. construct a detector- and finite-sample-aware optimal JSR estimator or prove a minimax obstruction.

Programs requiring real devices, traceable clocks, independent custody, laboratory records, or unblinding remain external gates and must not be replaced by synthetic evidence.

## Keywords

FIN; exact rational interval arithmetic; Krawczyk method; signed moment problem; Bernstein polynomial; Lean; cosine enclosure; spectral operator; legacy kernel; self-coupling audit; mathematical physics.
