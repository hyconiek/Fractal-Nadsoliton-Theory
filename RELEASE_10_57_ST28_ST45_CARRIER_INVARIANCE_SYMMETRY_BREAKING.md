# FIN Release 10.57 — ST28–ST45: Carrier-Invariant Information and Symmetry Breaking

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 1.0.0
- **Publication date:** 2026-08-10
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0
- **Repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory

## Abstract

Release 10.57 executes the twelve recommended programs ST28–ST39 and adds six programs, ST40–ST45, addressing two new foundational questions: whether information can remain invariant under a change of physical carrier, and whether the strong symmetry of the FIN strict operator can enable selection through a separate positive-feedback instability.

ST28 certifies the unique stability-collision root for the explicitly frozen decimal coefficient model in

```text
1.278014751804222 ≤ s* ≤ 1.2780147518366733.
```

ST29 proves that a state-dependent density can enlarge the strict invariant algebra from dimension 7 to dimension 74 for a reflection-anchored state and to the full dimension 144 for a generic state. The symmetry orbit, however, does not select one canonical member. ST30 proves that plain dyadic-refinement associativity leaves a new free parameter at every level. ST31 reports a genuine negative result: at the declared 2,500-count budget, the finite-count mixed-channel test has only 4.75% power against its specified alternative. ST32 proves a polynomial long-range influence bound without assuming Markov positivity. ST33 classifies the reflection-compatible cycle fluxes as 0 and π and proves that a link obtained from one state phase has trivial holonomy. ST34 proves a nonzero minimizer set for a constructed saturating DNLS model, but the numerical candidate is uniform rather than localized; no soliton evidence is claimed.

ST36 proves the remaining dimensional scale orbits. ST37 provides a frozen synthetic noise-and-clock-nuisance curve. ST38 proves a one-vector obstruction: no state in the declared construction class jointly supplies canonical selection, the full algebra, and nontrivial holonomy. ST39 organizes the current physical-completion requirements into nine relatively independent axiom groups.

ST40–ST42 give the main informational result. Simultaneous faithful transport of the generator, state, preparations, effects, composition, and records preserves every finite record probability to numerical residual `8.89e-16`. Transporting only an isospectral generator while leaving instruments fixed changes the record by total variation `0.1498688`. Therefore the candidate carrier-invariant object is the operational isomorphism class of the complete experiment, not the spectrum alone. Information can be independent of a particular faithful carrier, but operational information is not defined after every realization structure has been removed.

ST43–ST44 separate symmetry from feedback. Symmetry supplies twelve equivalent possibilities; positive gain is a distinct dynamical premise that amplifies a fluctuation, while saturation bounds growth. The strict operator canonically supplies a two-dimensional degenerate eigenspace, but no preferred axis in it. A conditionally selected axis generates only the 74-dimensional reflection-anchored algebra. ST45 proves that unitary, heat, and wave channels may share the same generator while remaining dynamically inequivalent.

No result in this release supplies laboratory data, SI calibration, a continuum, a strict-derived gain or nonlinear law, a canonical realized history, a legacy-to-strict role-transfer theorem, or Theory-of-Everything closure.

## Principal scientific conclusions

1. **Carrier independence has a precise mathematical meaning.** It is invariance under faithful transport of the complete operational process, not propagation without any realization.
2. **The spectrum is not a complete information invariant.** States, projectors, instruments, composition, and record semantics are necessary.
3. **Symmetry is not positive feedback.** Symmetry supplies degenerate alternatives; a gain law determines whether deviations grow.
4. **State dependence can escape the invariant-algebra no-go.** It does not by itself provide a canonical orbit member or nontrivial holonomy.
5. **The previous ideal finite-count receiver is too weak at the declared budget.** This is retained as a failed test, not reinterpreted as confirmation.
6. **The constructed saturating nonlinear completion did not produce localization.** Its best candidate is a uniform ordered phase.
7. **A shared spectral generator organizes several channels but does not make their dynamics identical.**

## Included files

- `FIN_ST28_ST45_Carrier_Invariance_and_Symmetry_Breaking_Report_EN.pdf` — 16-page English research report.
- `FIN_ST28_ST45_Carrier_Invariance_and_Symmetry_Breaking_Report_EN.tex` — report source.
- `fin_st28_st45_research.py` — executable analytical and numerical study.
- `FIN_ST28_ST45_Results.json` — complete machine-readable result ledger.
- `FIN_ST28_ST45_Summary.csv` — compact outcome table.
- `FIN_ST28_Jordan_Root_Certificate.json` — ST28 certificate data.
- `FIN_ST31_Finite_Count_Preregistration.json` — hashed ST31 synthetic protocol.
- `FIN_ST37_Nuisance_Preregistration.json` — hashed ST37 synthetic protocol.
- `test_fin_st28_st45.py` — twenty acceptance tests.
- `FIN_ST28_ST45_Figures/` — seven report figures.
- `FIN_ST28_ST45_INPUTS.sha256` — frozen source-input inventory.
- `FIN_RELEASE_10_57_MANIFEST.sha256` — release-file integrity manifest.

## Reproducibility

```bash
MPLCONFIGDIR=/tmp/fin-st28-mpl python3 fin_st28_st45_research.py
MPLCONFIGDIR=/tmp/fin-st28-mpl python3 -m unittest -v test_fin_st28_st45.py
pdflatex -interaction=nonstopmode -halt-on-error FIN_ST28_ST45_Carrier_Invariance_and_Symmetry_Breaking_Report_EN.tex
pdflatex -interaction=nonstopmode -halt-on-error FIN_ST28_ST45_Carrier_Invariance_and_Symmetry_Breaking_Report_EN.tex
```

The final local run passed all **20/20** tests. No network access, laboratory data, or external service is required.

## Recommended next programs

ST46–ST57 are ranked in the report. The highest proof-completion priority is full upstream interval regeneration for ST28. The highest conceptual priority is deriving or refuting the positive-feedback gain from the strict adaptive law and testing multi-component/projective objects for joint selector, algebra, and π-holonomy closure. Carrier invariance should continue through classification of operational intertwiners, a two-carrier synthetic transfer protocol, and entropy/data-processing invariants.

## Keywords

FIN; spectral operator; carrier invariance; operational information; representation independence; symmetry breaking; positive feedback; selector obstruction; state-dependent algebra; cycle holonomy; refinement; finite-count falsification; long-range dynamics; mathematical physics.
