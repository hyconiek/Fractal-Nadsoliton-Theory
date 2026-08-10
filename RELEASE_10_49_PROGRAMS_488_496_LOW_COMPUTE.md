# Release 10.49 — FIN Programs P488–P496

## Low-Compute Analytical Research

### Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint and reproducible computational package
- **Version:** 1.0.0
- **Publication date:** 2026-08-10
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0
- **Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Abstract

Release 10.49 executes the nine studies from the FIN Theory Compendium that
do not require large exact elimination, global cone optimization, laboratory
access, or external validation. The P485/P487 Gröbner campaign and the global
exact complex phase-face problem are deliberately excluded.

The release constructs a saturating strict spectral flow

\[
\dot\psi=(r-g|\psi|^2)\psi-iA\psi.
\]

Every constant-modulus strict Fourier eigenmode is an exact persistent
rotating state. Five excited modes have strictly negative transverse real
exponents in the numerical linearization; mode three has an additional
neutral direction. All modes have inverse participation ratio `1/12`, so the
result constructs persistent coherent cycles but not localized solitons.

Further results include:

- a relational spectral clock with dimensionless period
  `8.3317982424005` and an exact affine-scale obstruction;
- a positive three-pole Stieltjes learning flow recovering the declared
  Schur trace response to `1.31e-15`, with condition number `3.52e7`;
- exclusion of an exact dynamic fixed point in the declared normalized
  `C12/C24/C48/C96` attenuation-envelope family, whose distinct visible pole
  counts are `3/6/12/24`;
- a theorem that maximum coupling resonance of `W` is exactly minimum
  Dirichlet cost of `A=sI-W`, while maximum temporal frequency is a different,
  oppositely ordered criterion;
- an exact exponent-lift kernel automaton containing Legacy* and strict as
  supplied parameter points, without deriving the parameter transition;
- a topological obstruction proving that scalar `C^12` states do not carry
  protected winding sectors without a nonvanishing phase field and refinement
  rule;
- a synthetic operational atlas in which coherence-sensitive Fourier records
  raise four-model classification accuracy from `0.5775/0.745` to `1.0`;
- a dependency-minimal exact rational replay of the Dirichlet, resonance,
  Green–Schur, and Legacy* recurrence identities.

The complete numerical batch ran locally in approximately 32 seconds. Ten
regression tests passed in 0.18 seconds.

## Main document

- `FIN_Programs_488_496_Low_Compute_Analytical_Report_EN.pdf`
- `FIN_Programs_488_496_Low_Compute_Analytical_Report_EN.tex`

## Reproducibility package

- `fin_programs_488_496_low_compute.py`
- `fin_programs_488_496_exact_checker.py`
- `test_fin_programs_488_496.py`
- `FIN_Programs_488_496_Results.json`
- `FIN_Programs_488_496_Exact_Check.json`
- `FIN_Programs_488_496_Summary.csv`
- `FIN_Programs_488_496_Figures/`
- `FIN_PROGRAMS_488_496_RELEASE_MANIFEST.sha256`

## Scientific boundary

This release supplies conditional finite dynamical models, exact algebra, and
synthetic identifiability results. It does not derive a localized nadsoliton,
an absolute time unit, the adaptive target/loss, a unique legacy-to-strict
parameter law, particle topology, physical units, a selector, laboratory
evidence, Standard Model, gravity, total action, or a Theory of Everything.

## Keywords

spectral limit cycle; relational time; Stieltjes function; adaptive memory;
Schur complement; fractal fixed point; resonance; graph Dirichlet energy;
topological winding; operational identifiability; FIN.
