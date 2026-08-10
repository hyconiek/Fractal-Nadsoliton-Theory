# Release 10.52 — Programs P517–P526

## Fourth-Jet Obstruction, Stability Boundary, and the Six-Frequency Torus

This release executes the ten bounded local programs recommended by Release
10.51. It tests whether FIN's information language supplies the missing
focusing nonlinearity, upgrades the finite-DNLS stability analysis to outward
interval certificates, distinguishes stationary memory persistence from
temporal stability, and strengthens the exact arithmetic of the strict
spectrum.

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

P517 proves a family-level obstruction: entropy, Fisher information and
compression labels do not determine the focusing fourth-order jet required by
the localized DNLS state. A reference state, extremum direction and quartic
coupling normalization remain additional typed data. P518 provides an exact
rational replay of all 801 exported Krawczyk acceptance inequalities. P519
constructs 208 outward state charts and 207 nested interface certificates over
`0.68 <= omega <= 1.20`, proves the finite `L-/L+` inertia ledger, and isolates
the unique VK sign transition in
`[0.722143668857225, 0.722143708857225]`.

P520 identifies the bond-to-uniform transition with the mode-one uniform
bifurcation at `omega = 0.3770605771035399`; no scanned phase kick remains both
localized and coherently translating. P521 supplies an explicit causal
three-pole hidden-relaxation realization of the stationary Stieltjes memory and
finds a positive growth exponent already at loading `0.01`. It therefore
falsifies the inference from stationary persistence to temporal stability for
that realization. P522 proves that all six distinct positive strict
frequencies are linearly independent over the algebraic numbers; the exact
distance-to-mode determinant is `-3456*sqrt(3)`. The same program proves an
absolute-scale no-go under `(A,t) -> (cA,t/c)`.

P523 transfers a fractional `n^-0.8` context certificate to `C384`, with
maximum normalized bound `0.0670791` under a declared FFT rounding model. P524
exports a dependency-minimal exact replay of the complete 573-box PSD ledger.
P525 constructs a seven-coordinate synthetic detector-error polytope. P526
shows that frozen hold-out ranking changes with regularization and proves a
finite-record no-free-lunch theorem for unrestricted parameter laws.

## Main epistemic results

1. **P517 / O221 — information-functional fourth-jet obstruction [Proven].**
   The finite strict quadratic core and standard information labels do not
   source an attractive quartic response.
2. **P518 / O222 — exact O212 replay [Proven].** All 401 parameter charts and
   400 interface boxes pass exact dyadic-rational acceptance replay.
3. **P519 / O223 — interval stability tube [Proven, conditional on DNLS].** A
   connected certified branch covers the full declared frequency interval;
   the favourable VK sector contains `omega=1`.
4. **P520 / O224 — bond and phase-kick audit [Strong evidence / refuted
   candidate].** The uniform bifurcation explains the bond transition; the
   scanned kicks do not yield a localized translating object.
5. **P521 / O225 — temporal memory realization [Conditional / refuted
   candidate].** The first explicit hidden-relaxation completion is linearly
   unstable even though its stationary localized branch persists.
6. **P522 / O226 — maximal frequency independence [Proven].** Six positive
   frequencies are algebraically linearly independent, but cannot select an
   absolute second or energy unit.
7. **P523 / O227 — C384 fractional certificate [Conditional].** The bound is
   useful but retains an explicit finite-FFT rounding premise.
8. **P524 / O228 — portable PSD replay [Proven].** Exact standard-library
   checks verify the partition and every exported sign condition.
9. **P525 / O229 — detector-error polytope [Conditional].** The nonempty
   region is a synthetic operational budget, not apparatus calibration.
10. **P526 / O230 — frozen-holdout law audit [Proven/conditional].** Finite
    candidate ranking is prior-sensitive, and unrestricted finite-record law
    recovery is impossible.

## Included files

- `FIN_Programs_517_526_Fourth_Jet_Stability_and_Clock_Report_EN.pdf`
- `FIN_Programs_517_526_Fourth_Jet_Stability_and_Clock_Report_EN.tex`
- `fin_programs_517_526_research.py`
- `fin_p518_p524_stdlib_checker.py`
- `test_fin_programs_517_526.py`
- `FIN_Programs_517_526_Results.json`
- `FIN_Programs_517_526_Summary.csv`
- `FIN_P518_Krawczyk_Replay_Certificate.json`
- `FIN_P524_PSD_Replay_Certificate.json`
- `FIN_Programs_517_526_Figures/`
- `FIN_PROGRAMS_517_526_RELEASE_MANIFEST.sha256`

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-mpl-1052 python3 fin_programs_517_526_research.py
python3 fin_p518_p524_stdlib_checker.py
MPLCONFIGDIR=/tmp/fin-mpl-1052 python3 -m unittest -v test_fin_programs_517_526.py
pdflatex -interaction=nonstopmode -halt-on-error \
  FIN_Programs_517_526_Fourth_Jet_Stability_and_Clock_Report_EN.tex
sha256sum -c FIN_PROGRAMS_517_526_RELEASE_MANIFEST.sha256
```

The complete research script ran locally in approximately 24.2 seconds.
Twelve regression and certificate tests passed.

## Recommended next batch

The monograph specifies P527–P536. The leading priorities are P527 (an
axiomatic representation theorem for the fourth-order information jet), P529
(validated nonlinear dynamics on both sides of the VK boundary), and P530 (a
classification of passive, Hamiltonian and metriplectic temporal memory
realizations).

## Scope and non-claims

This is finite, dimensionless mathematics and computation. It does not derive
the focusing law from FIN, discharge QW-2191, select canonical orientation,
provide physical units, complete the historical legacy-to-strict bridge,
transfer legacy physical claims, furnish an apparatus or independent record,
derive a Standard Model or gravitational sector, or establish a Theory of
Everything.

