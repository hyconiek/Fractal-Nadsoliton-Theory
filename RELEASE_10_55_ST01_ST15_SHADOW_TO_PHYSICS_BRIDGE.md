# Release 10.55 — ST01–ST15 Shadow-to-Physics Bridge Research

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint and accompanying research software
- **Version:** 1.0.0
- **Publication date:** 2026-08-10
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0
- **Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Abstract

Release 10.55 executes the fifteen local analytical and computational studies
ST01–ST15 proposed by the preceding counterfactual FIN-as-ToE shadow atlas. The
studies ask which parts of the frozen twelve-state strict FIN operator can be
carried toward physical interpretation and which additional structures remain
logically independent.

The release contains a typed physical-shadow certificate, common-generator
round trips across five spectral channels, five classes of null operators,
algebra-completion calculations, explicit incompatible continuum
continuations, propagation-tail theorems, a conditional
information-to-thermodynamics bridge, a gauge-covariant U(1) receiver, a new
globally coercive saturating-mediator completion, a Hamiltonian-memory/Krein
audit, a Schur-reduction audit, a legacy-role transfer gate, frozen synthetic
cross-channel tests, a hashed prediction ledger, and a six-part formal-core
replay.

The central conclusion is negative but constructive. No single missing theorem
turns the current finite spectral object into physics. The remaining bridge
separates into independently necessary selector/algebra, dimensional,
refinement, operational, and external-record obligations. The report therefore
does not establish FIN as a Theory of Everything.

## Principal results

1. **Typed physical-shadow certificate — proven schema.** A proposed physical
   reduction must type its sector, units, refinement/coarse graining,
   observables, preparations, environment, instruments/apparatus, composition,
   and record. The composition identity for approximate intertwiners is proved.

2. **Common spectral generator — closed-loop validation.** Unitary, heat,
   Green, Gibbs, and wave transforms recover the frozen generator with relative
   matrix errors below `4.3e-15` on certified branches. Altered-generator and
   aliased-unitary controls are rejected. These are generated-transform round
   trips, not physical tomography.

3. **Isospectral obstruction.** In a 6,000-operator, five-ensemble null atlas,
   orthogonal isospectral rotations preserve every eigenvalue to `4.9e-15`
   while generically destroying FIN vertex geometry and the Markov/M-matrix
   interpretation. The spectrum alone cannot identify the substrate.

4. **Conditional algebra completion.** `C*(A)` has dimension 7. Adding one
   diagonal observable with twelve distinct vertex labels generates
   `M_12(C)`, dimension 144. This supplies rather than derives labels and
   orientation, so it does not solve the selector obstruction.

5. **Continuum nonidentifiability.** Two continuations agree exactly at
   `N=12` but yield fitted gap exponents `-1.96485` and `-0.08635`. The current
   endpoint therefore cannot choose a unique continuum family. A literal
   integer-distance continuation first develops a negative weight at distance
   8 and loses the Markov property.

6. **No exact cyclic propagation cone.** Because the frozen strict kernel has
   nonzero all-distance couplings, remote unitary and heat amplitudes appear at
   first order and wave amplitudes at second order. Approximate locality in a
   future scaled refinement remains open.

7. **Conditional physics receivers.** The release proves the standard
   relative-entropy/free-energy identity after supplying thermodynamic
   semantics, and constructs a gauge-covariant U(1) link receiver after
   supplying link phases. Neither construction sources its physical objects
   from strict FIN.

8. **New coercive nonlinear completion.** A saturating auxiliary-field coupling
   preserves the attractive fourth-order jet of Program P527 and produces a
   finite-dimensional energy bounded below and coercive. The saturation map,
   scale, and positive coercive term are new axioms, not strict-core outputs.

9. **Memory stability remains open.** The implemented Hamiltonian-memory
   Jacobian has exact block form `J=J0 L`. A scan locates a numerical transition
   in mediator speed within `[1.24, 1.28]`, but no exact Krein or collision
   certificate is claimed.

10. **Prediction boundary.** A preregistered kernel-law fingerprint has zero
    false accepts in a separate 5,000-instance synthetic holdout, but it merely
    tests the defining strict formula. It is useful for code integrity and is
    not an experimentally meaningful prediction.

## Evidence inventory

- 15 executed research programs;
- 6,000 null-ensemble matrices in five declared classes;
- 5,000 frozen synthetic holdout matrices after a 2,000-sample training phase;
- 6 generated figures;
- 1 formal-core certificate with six accepted checks;
- 1 hashed preregistration record;
- 17 regression/acceptance tests, all passing.

## Files included

- `FIN_ST01_ST15_Shadow_to_Physics_Bridge_Research_Report_EN.pdf` — complete
  English research report;
- `FIN_ST01_ST15_Shadow_to_Physics_Bridge_Research_Report_EN.tex` — report
  source;
- `fin_st01_st15_research.py` — deterministic analytical and numerical study;
- `test_fin_st01_st15.py` — live acceptance tests;
- `FIN_ST01_ST15_Results.json` — machine-readable results;
- `FIN_ST01_ST15_Summary.csv` — compact program ledger;
- `FIN_ST01_ST15_Formal_Core_Certificate.json` — ST15 replay certificate;
- `FIN_ST14_Prediction_Preregistration.json` — frozen ST14 record;
- `FIN_ST01_ST15_Figures/` — six figures;
- `FIN_ST01_ST15_INPUTS.sha256` — hashes of principal source inputs;
- `FIN_ST01_ST15_RELEASE_MANIFEST.sha256` — release-integrity manifest.

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-st-mpl python3 fin_st01_st15_research.py
MPLCONFIGDIR=/tmp/fin-st-test-mpl python3 -m unittest -v test_fin_st01_st15.py
pdflatex -interaction=nonstopmode -halt-on-error FIN_ST01_ST15_Shadow_to_Physics_Bridge_Research_Report_EN.tex
sha256sum -c FIN_ST01_ST15_RELEASE_MANIFEST.sha256
```

The generated transforms, ensembles, and holdout are deterministic under the
recorded seed. Numerical evidence remains dependent on the documented binary64
and library environment. The symbolic/rational replay is not proof-assistant
machine checking.

## Recommended next programs

ST16–ST27 are ranked in the report. The leading four are:

1. an exact Krein/zero-mode certificate for the Hamiltonian memory family;
2. classification of selector-free additional observable algebras;
3. search for a FIN-internal refinement functor;
4. a genuinely derived, preregistered cross-sector observable replacing the
   definitional kernel fingerprint.

## Scope and non-claims

Release 10.55 supplies no canonical selector or orientation, physical clock,
absolute units, FIN-internal continuum, legacy-to-strict completion,
legacy-role transfer theorem, laboratory apparatus or record, independent
custody, Standard Model sector, gravitational sector, total physical
Lagrangian, or Theory-of-Everything closure.

All PDF reports generated after this release are to be written in English
unless a later request explicitly specifies another language.
