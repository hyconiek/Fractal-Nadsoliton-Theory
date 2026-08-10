# Release 10.51 — Programs P507–P516

## Nonlinear-Source Obstruction, Certified Localization, and Exact Frequency Arithmetic

This release records ten new local analytical and computational studies of the
FIN strict-kernel programme.  The batch deliberately separates what follows
from the supplied finite operator from what still requires an additional
dynamical, informational, operational, or dimensional premise.

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint and research software
- **Version:** 1.0.0
- **Publication date:** 2026-08-10
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0
- **Repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory

## Abstract

Programs P507–P516 test the nonlinear and operational frontier opened by the
finite strict operator.  P507 proves an obstruction: the quadratic strict core
fixes only the Hessian of an action and cannot determine the sign or coefficient
of a cubic DNLS evolution.  P508 then treats the focusing DNLS as an explicit
additional premise and certifies, with 401 interval Krawczyk charts and 400
interface certificates, one unique fold-free real localized branch over the
complete coupling interval.  P512 closes the earlier finite-denominator clock
test by proving that the declared first two nonzero strict frequencies have a
transcendental, hence irrational, ratio.  The remaining programs delimit
orbital-stability evidence, the failure of a localized equal-power
Peierls–Nabarro comparator, survival under a stationary memory loading, a
conservative context-limit error envelope, a global interval PSD phase
diagram, a synthetic robustness budget, and conditional identifiability of
finite-degree parameter laws.  None of these results supplies dimensional
units, selects an operational clock, completes the legacy-to-strict bridge,
transfers historical physical roles, or constitutes laboratory evidence.

## Main results

1. **P507 / O211 — quartic-jet source obstruction [Proven].**  Two actions can
   have the same value, gradient, and Hessian at the vacuum while generating
   opposite cubic nonlinearities.  The strict quadratic core therefore cannot
   select focusing versus defocusing DNLS or its nonlinear coefficient.
2. **P508 / O212 — parametric Krawczyk continuation tube [Proven, conditional
   on the supplied DNLS].**  All 401 parameter charts pass strict inclusion;
   all 400 shared-parameter root boxes are nested in both adjacent charts.  The
   largest Krawczyk defect is 0.006276 and the smallest strict margin is
   (1.33\times10^{-7}).
3. **P509 / O213 — finite-DNLS orbital-stability ledger [Strong numerical
   evidence].**  At the declared frequency the finite GSS/VK sign ledger is
   favourable, but the power slope changes sign near 0.72214.  Stability is
   therefore not licensed over the whole frequency family.
4. **P510 / O214 — equal-power PN audit [Refuted/undefined in the declared
   comparison].**  The localized bond branch does not reach the site-centred
   target power before merging into the uniform branch, so the conventional
   localized equal-power barrier is not defined.
5. **P511 / O215 — memory-loaded localized branch [Strong numerical evidence,
   conditional].**  Localization survives the declared stationary
   Stieltjes–Bernstein memory loading; this is not a temporal non-Markovian
   stability theorem.
6. **P512 / O216 — exact frequency-ratio theorem [Proven].**  A
   Lindemann–Weierstrass/Laurent-polynomial argument proves that
   λ₂/λ₁ is transcendental.  This gives nonperiodic relative phase dynamics
   but not an absolute time unit or an internally selected clock.
7. **P513 / O217 — context error envelope [Partial certificate].**  A rigorous
   conservative (O(n^{-4/5})) bound is obtained, but it becomes informative
   only at substantially larger context sizes than the observed numerical
   convergence.
8. **P514 / O218 — global interval PSD diagram [Proven for the declared
   homotopy].**  Adaptive interval subdivision leaves only a
   (3.50\times10^{-10})-wide transition bracket and excludes hidden later PSD
   losses.
9. **P515 / O219 — one-context robustness budget [Conditional].**  The chosen
   synthetic context retains positive pairwise separation under stated
   preparation, measurement, and timing budgets.  These budgets are not lab
   calibrations.
10. **P516 / O220 — polynomial parameter-law design [Conditional].**  Constant,
    linear, and quadratic parameter velocities become locally identifiable
    from finite observations only after the finite-degree law class is
    supplied; the unrestricted bridge-law no-go remains.

## Included files

- `FIN_Programs_507_516_Nonlinear_Source_and_Certified_Localization_Report_EN.pdf`
- `FIN_Programs_507_516_Nonlinear_Source_and_Certified_Localization_Report_EN.tex`
- `fin_programs_507_516_research.py`
- `test_fin_programs_507_516.py`
- `FIN_Programs_507_516_Results.json`
- `FIN_Programs_507_516_Summary.csv`
- `FIN_Programs_507_516_Figures/`
- `FIN_PROGRAMS_507_516_RELEASE_MANIFEST.sha256`

## Reproduction

```bash
python3 fin_programs_507_516_research.py
python3 -m unittest -v test_fin_programs_507_516.py
pdflatex -interaction=nonstopmode -halt-on-error \
  FIN_Programs_507_516_Nonlinear_Source_and_Certified_Localization_Report_EN.tex
sha256sum -c FIN_PROGRAMS_507_516_RELEASE_MANIFEST.sha256
```

The principal script finished locally in approximately 57.1 seconds.  The
release test suite contains eleven tests.

## Recommended next batch

The monograph defines P517–P526.  The leading priorities are P517 (test whether
an information functional can supply the missing fourth-order jet), P519
(interval certification of the finite orbital-stability ledger), and P521
(temporal stability for the memory-loaded model).  These are bounded local
tasks and do not require laboratory access or external audit.

## Scope and non-claims

This is a finite, dimensionless mathematical and computational release.  It
does not derive the nonlinear law from FIN, discharge QW-2191, select a
canonical orientation, provide SI units, complete the historical
legacy-to-strict bridge, transfer legacy physical claims, furnish an apparatus
or independent record, derive a Standard Model or gravitational sector, or
establish a Theory of Everything.

