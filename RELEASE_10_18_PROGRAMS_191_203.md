# Release 10.18 — Reference States, Operational Quotients, and Conditional Prediction

## Zenodo release description

**Title:** FIN Programs 191–203: Reference States, Operational Quotients, and Conditional Prediction

**Creator:** Żuchowski, Krzysztof  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Resource type:** Publication — Preprint  
**Version:** 1.0.0  
**Publication date:** 2026-07-27  
**Language:** English  
**Publisher:** Zenodo  
**Access:** Open access  
**License:** Creative Commons Attribution 4.0 International (CC BY 4.0)

## Abstract

Release 10.18 reports thirteen executed mathematical and operational research
programs, numbered 191 through 203. The studies investigate the remaining
state, scale, environment, apparatus, record, and prediction problems around
the dimensionless FIN spectral core.

The release proves that inner-unitary naturality selects normalized trace on
each simple matrix block but leaves an independent probability measure on the
centre of a direct-sum algebra. For the FIN dimension vector
`(1,2,2,2,2)`, uniform central weighting yields `9/5`, while normalized trace
on the represented space yields `17/9`; neither is uniquely selected by the
current axioms.

The positive compression parameter is quotiented under length rescaling. The
complete reduced kernel profile is

```text
k_hat(x) = cos(nu*x + phi)/(1 + x^eta)
x = beta^(1/eta) d
nu = omega/beta^(1/eta)
```

This establishes `beta=1` as a coordinate representative, not a strict
source. Gauge-invariant legacy and strict shape data remain different, so the
quotient does not complete the legacy-to-strict bridge.

The release’s main analytic result is a closed-form boundary for
deterministic one-qubit reflection-covariant conversion. The
Alberti–Uhlmann all-parameter trace-norm criterion reduces to three candidate
points:

```text
r_target,max^2 =
  min_{x in {0,1,r_source^2/(1-z_source^2)}}
  [max(x,r_source^2+x*z_source^2)-x*z_target^2].
```

It exactly recovers the previous counterexample: source `(r,z)=(0.6,0)`
cannot reach target `(0.6,0.8)` because the maximum reachable transverse
coordinate is `0.36`.

A second theorem proves that single, plus, and echo visibility contrasts form
a minimal rank-three tomography frame for a two-step Gaussian dephasing
process with unknown detector blur. The frame identifies blur, increment
variance, and temporal covariance; removing echo lowers the rank to two.

Further contributions include:

- a nonasymptotic DKW theorem for an observed regenerative beta-mixing
  acquisition class;
- a preregistered order-sensitive open-set detector;
- an operational environment equivalence quotient;
- an explicit pair of microscopic environments with the same one-time
  channel and different two-time process;
- a natural analytic phase-source no-go theorem;
- a catalogue of projective scale-free spectral observables;
- an event-level synthetic operational bundle passing ten of eleven intake
  fields;
- a held-out `W0 + CA + OP` conditional prediction dry run.

Arb/python-flint was not available for the requested common
directed-rounding certificate. A general finite Dirichlet Lean source is
included, but no local Lean toolchain was installed, so machine compilation
is not claimed. No independent external dataset passed the intake gate, and
the final prediction is therefore method validation rather than empirical
physics.

## Main scientific results

1. **Central-state classification — Proven**
   - block state: normalized trace;
   - remaining datum: a central probability measure;
   - exact competing expectations: `9/5` and `17/9`.

2. **Beta-gauge quotient — Proven**
   - maximum orbit-collapse residual: `1.11e-16`;
   - gauge-invariant legacy/strict RMS profile gap: `0.3816158693`;
   - bridge closure: not obtained.

3. **Common Arb certificate — Open**
   - one of five dependency nodes closed in the inherited chain;
   - inherited worst formal enclosure: `0.6129831478`;
   - target: `0.03`;
   - no false formal certificate is asserted.

4. **Finite Dirichlet mathematics — Proven on paper and by exact checks**
   - 150 exact rational cycle cases;
   - positivity, constant kernel, unitary exponential, positive stochastic
     heat semigroup;
   - Lean source included but not machine compiled.

5. **Reflection-conversion geometry — Proven**
   - numerical maximization replaced by a finite analytic formula;
   - 500 random cross-checks;
   - maximum discrepancy from the former dense Choi grid: `6.68e-6`.

6. **Regenerative dependent inference — Proven in the declared class**
   - nominal-record coverage: `0.85625`;
   - valid regenerative coverage: `1.0`;
   - mean retained fraction: `0.05029`.

7. **Order-sensitive open-set challenge — Strong synthetic evidence**
   - iid Gaussian rejection: `0.01429`;
   - iid stable rejection: `0.00714`;
   - sorted, AR(1), repeated-block and drift rejection: `1.0`.

8. **Minimal process tomography — Proven in the declared family**
   - full design rank: `3`;
   - rank without echo: `2`;
   - condition number: `7.2874`.

9. **Environment nonidentifiability — Proven**
   - common one-time dephasing factor: `0.8`;
   - two-time characteristics: `0.28` versus `1.0`.

10. **Analytic phase provenance — Proven no-go in declared classes**
    - continuous circle endomorphisms map roots of unity to roots of unity;
    - target-producing analytic candidates are target-coded or branch-coded.

11. **Scale-free observable algebra — Proven**
    - raw gap changes by `10^12`;
    - normalized observable residuals remain below `8e-14`.

12. **Operational intake — Proven audit**
    - synthetic method bundle: `10/11`;
    - failed field: independent analysis boundary;
    - external bundle admitted: no.

13. **Conditional prediction — Method validation**
    - held-out predicted visibility: `0.56121`;
    - held-out synthetic visibility: `0.59768`;
    - residual: `0.03646`;
    - external physical test: not executed.

## Scientific interpretation

The same self-adjoint generator continues to support wave, heat, and Green
dynamics. The new results show why this common engine is not by itself a
physical experiment. A non-simple observable algebra requires a central
state; positive compression requires a scale section; reduced dynamics
identify only an environment quotient; memory requires interventions;
measurement requires a POVM and detector calibration; physical time requires
a clock; and validation requires an independent ordered record.

The strongest interpretation surviving falsification is:

> FIN is a dimensionless spectral-information core whose route to physics is
> an operational quotient-and-section problem. The core determines
> nontrivial dynamics and conversion geometry, while state, scale,
> environment, apparatus, and empirical record remain independent typed
> objects until fixed by theorem or experiment.

## Scope guardrails

This release does **not** claim:

- a strict central state or temperature source;
- a target-independent strict `beta` source;
- an Arb-certified sub-three-percent fractional bound;
- machine compilation of the Lean library;
- a non-premise selector or discharge of `QW-2191`;
- a canonical physical unit;
- a completed legacy-to-strict bridge;
- transfer of legacy physical roles;
- a role-bearing `L_total`;
- Standard Model or general-relativistic generation;
- external empirical validation;
- Theory-of-Everything closure.

## Included files

- `FIN_Programs_191_203_Reference_Process_Prediction_Monograph.pdf`
  — 57-page English preprint.
- `FIN_Programs_191_203_Reference_Process_Prediction_Monograph.md`
  — editable source.
- `FIN_Programs_191_203_Reference_Process_Prediction_Monograph.tex`
  — archival LaTeX.
- `FIN_Programs_191_203_Reference_Process_Prediction_Results.json`
  — machine-readable results and claim boundaries.
- `fin_programs_191_203_reference_process_prediction.py`
  — executable research suite.
- `test_fin_programs_191_203_reference_process_prediction.py`
  — 65 regression and scientific-boundary tests.
- `FIN_Programs_191_203_Order_Sensitive_Preregistration.json`
  — frozen Program 197 protocol.
- `FIN_Programs_191_203_Synthetic_Operational_Bundle.csv`
  — ordered event-level method-validation record.
- `FIN_Programs_191_203_Synthetic_Operational_Bundle.json`
  — bundle metadata and immutable data hash.
- `FIN_Programs_191_203_Dirichlet_Library.lean`
  — uncompiled formalization source.
- `FIN_Programs_191_203_Reference_Process_Prediction_Figures/`
  — thirteen figures.
- `FIN_Programs_191_203_Release_10_18_Manifest.json`
  — archival hashes.

## Reproduction

```bash
python3 fin_programs_191_203_reference_process_prediction.py
python3 -m unittest -v \
  test_fin_programs_191_203_reference_process_prediction.py
python3 fin_programs_191_203_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error \
  FIN_Programs_191_203_Reference_Process_Prediction_Monograph.tex
```

The result suite uses seed `20260727`. It does not use Firecrawl, web search,
or external datasets.

## Next research programs

Programs 204–216 are ranked in the monograph. Their leading targets are:

1. Morita-natural central-measure classification;
2. a target-independent `eta=9/5` renormalization cocycle;
3. a reproducible Arb certificate environment;
4. a machine-compiled Dirichlet proof library;
5. tensor conversion and catalysis;
6. hidden-refresh confidence sequences;
7. sequential temporal conformal calibration;
8. finite-shot process-tomography confidence regions;
9. trigonometric moment bounds for environment equivalence classes;
10. formal analytic phase-source no-go;
11. a scale-free falsification protocol;
12. an independent 11/11 experimental bundle;
13. the first held-out conditional prediction on admitted external data.

## Keywords

FIN; spectral operator; operator algebras; central states; normalized trace;
spectral graph theory; Dirichlet forms; reflection covariance;
Alberti–Uhlmann theorem; quantum resource theory; beta mixing; empirical
processes; open-set detection; process tensors; quantum channels;
Stinespring dilation; environment identifiability; scale invariance;
operational physics; conditional prediction; information dynamics.

## Suggested citation

Żuchowski, K. (2026). *FIN Programs 191–203: Reference States, Operational
Quotients, and Conditional Prediction* (FIN Research Monograph, Release
10.18; Version 1.0.0) [Preprint]. Zenodo.

