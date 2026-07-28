# Release 10.21 — FIN Laboratory Transfer Package

## Title

**FIN Laboratory Transfer Package: P240 Optimal Tomography, P241 Validator,
and P242 One-Shot Holdout Pipeline**

## Metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Technical report and executable
  specification
- **Version:** 1.0.0 (FIN release series 10.21)
- **Publication date:** 2026-07-27
- **Language:** English
- **License:** CC BY 4.0 for the report; repository code remains subject to the
  repository license
- **Repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory

## Abstract

Release 10.21 is a single scientific transfer package for theoretical and
experimental physicists.  It translates the finite FIN strict-generator model
into a falsifiable laboratory protocol without claiming that a realization
already exists.

The package contains:

1. **P240:** a frozen twelve-state heat-kernel tomography design, including a
   dimensionless optimal acquisition time and nonasymptotic shot guidance;
2. **P241:** a fail-closed external-bundle validator for preparations, clock,
   vertex-resolving measurement, finite event records, and independent
   custody;
3. **P242:** a one-shot semigroup and projective-spectrum holdout pipeline that
   cannot execute before the external gate passes;
4. apparatus hypotheses for a photonic graph walk, qubit interferometry with
   dephasing, and an event-level double-slit control;
5. a blind-custody protocol separating provider, registrar, and analyst roles;
6. machine-readable locks, schemas, fixtures, figures, and regression tests.

The principal candidate is a photonic continuous-time quantum walk on a
twelve-vertex graph, subject to an experimental proof that the implemented
preparations, generator, clock, and vertex POVM match the declared model.
The double-slit platform is retained as an event-level coherence and custody
control; by itself it cannot validate the twelve-state heat-semigroup claim.

## Scientific boundary

The report proves conditional reconstruction and software validation
properties.  It does not prove that any laboratory apparatus realizes FIN,
does not create independent custody by code, and does not supply a physical
clock, selector, dimensional action, completed kernel bridge, or
experimentally validated physical theory.

## Primary files

- `FIN_Laboratory_Transfer_Package_P240_242.pdf`
- `FIN_Laboratory_Transfer_Package_P240_242.md`
- `fin_lab_p240_optimal_tomography.py`
- `fin_lab_p241_validator.py`
- `fin_lab_p242_pipeline.py`
- `test_fin_lab_p240_242.py`
- `FIN_Lab_P240_242_Transfer_Manifest.json`

## Keywords

heat-kernel tomography; continuous-time quantum walk; photonic graph;
Markov semigroup; projective spectrum; blind custody; preregistration;
event-level double slit; falsification; FIN.
