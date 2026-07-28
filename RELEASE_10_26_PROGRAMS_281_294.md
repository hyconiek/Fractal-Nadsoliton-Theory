# Release 10.26 — FIN Research Programs P281–P294

## Regularized spectral inference, resource-explicit physics, and the pointed-torsor boundary

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 1.0.0
- **Publication date:** 2026-07-29
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0
- **Repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory

## Abstract

Release 10.26 executes fourteen bounded research programs, P281–P294, at the
current FIN proof frontier. It studies regularized Stieltjes inversion, a
machine-checked Schur witness, operational current measurement, coarse
completion maps, multitime information accounting, one legacy-sourced fractal
compression candidate, empirical admission, detector-aware experimental
design, null-ensemble false positives, mechanism discrimination, intervention
design, reservoir replication, resource-explicit information erasure, and a
finite torsor obstruction to selecting an operational reference.

The central constructive result is the **Pointed Physicalization Resource**

\[
\mathfrak P=(T,t_0).
\]

For the tested finite regular torsor, equivariance alone provides no section,
whereas pointing yields a unique pointed equivariant map. The result identifies
a precise operational boundary: a prepared reference can remove a symmetry
ambiguity, but the symmetric spectral structure does not derive that reference.

The release does not claim physical confirmation. It preserves the strict /
legacy kernel split, leaves QW-2191 open, transfers no legacy physical role,
and derives no dimensional constant, Standard Model sector, general
relativity, complete Lagrangian, or Theory of Everything.

## Principal results

1. **P281:** compact pole and residue constraints guarantee existence of a
   constrained Stieltjes fit but do not remove inverse ill-conditioning.
2. **P282:** Lean 4.28 verifies a dependency-free exact rational
   \(3\times3\) nested/direct Schur witness with value \(13/5\).
3. **P283:** the minimal current POVM has an explicit Naimark isometry; a
   frozen loss/crosstalk/dark-count channel preserves response rank five.
4. **P284:** an isotropic complement restores coarse lifted rank and the null
   mode, but fails the strict fingerprint and is not an RG fixed point.
5. **P285:** a two-step quantum process tree closes an exact relative
   information ledger under supplied instruments and channels.
6. **P286:** the single canonical legacy fractal-compression atom
   \(d\mapsto d^{4\ln2}\) is insufficient because its relevant spectra remain
   negative.
7. **P287:** no production P241 laboratory bundle is present; P242 remains
   unauthorized.
8. **P288:** an event-level detector model gives an interior synthetic
   two-flux optimum \(h_\star=0.148111\).
9. **P289:** under the declared
   \(\operatorname{Dirichlet}(1,\ldots,1)\) null, importance sampling places
   the strict-fingerprint false-positive probability near \(2.9\times10^{-8}\)
   at tolerance \(0.02\).
10. **P290:** three times plus one intervention distinguish four frozen
    mechanisms with 96.9% accuracy at reference noise 0.20.
11. **P291:** a balanced long two-level pulse maximizes the weakest local
    information direction in the frozen adaptive model.
12. **P292:** the NARMA10 advantage only partially replicates—13 of 24
    independent realizations beat every matched control—and is not universal.
13. **P293:** a supplied collision-model erasure protocol closes its
    work/heat/information ledger while exposing every dimensional resource.
14. **P294:** a pointed finite torsor yields a unique pointed map; the point
    remains additional structure.

## Methodological correction

The first P292 implementation inherited a NARMA driving convention that
overflowed in some replications. Those non-finite trials were rejected. The
program was rerun with an independent stable NARMA input range
\([0,0.40]\). This correction is recorded because invalid numerical trials
cannot be used as evidence.

## Included files

- `FIN_Programs_281_294_Research_Report_EN.pdf` — archival English monograph;
- `FIN_Programs_281_294_Research_Report_EN.md` — editable report source;
- `FIN_Programs_281_294_Research_Report_EN.tex` — generated LaTeX source;
- `fin_programs_281_294.py` — deterministic research runner;
- `test_fin_programs_281_294.py` — 19-test verification suite;
- `FIN_Programs_281_294_Schur_Core.lean` — Lean 4.28 certificate;
- `FIN_Programs_281_294_Results.json` — machine-readable results;
- `FIN_Programs_281_294_Summary.csv` — compact program matrix;
- five detailed CSV audit tables;
- five figures in `FIN_Programs_281_294_Figures/`;
- `FIN_Programs_281_294_SHA256SUMS.txt` — integrity manifest.

## Verification

Run:

```bash
MPLCONFIGDIR=/tmp/matplotlib-p281-294 \
python3 -m unittest -v test_fin_programs_281_294.py
```

Expected result: **19 passed, 0 failed**. The suite recompiles the Lean
certificate.

## Recommended continuation

The report specifies fourteen programs P295–P308. The highest-priority
mathematical sequence is:

\[
\mathrm{P295}\to\mathrm{P296}\to\mathrm{P303}
\to\mathrm{P304}\to\mathrm{P298},
\]

covering operator-valued Stieltjes stability, a general machine-checked Schur
theorem, analytic rare-event asymptotics, adversarial mechanism
identifiability, and classification of non-target-coded complement RG maps.

The external empirical lane remains blocked until a genuinely independent,
hash-committed P241 bundle is supplied.

## Keywords

FIN; spectral operator; Stieltjes transform; inverse problem; Schur
complement; Lean; POVM; Naimark dilation; detector tomography; renormalization;
quantum process; relative entropy; intervention design; reservoir computing;
Landauer erasure; torsor; operational physics; mathematical physics.

## Suggested citation

Żuchowski, K. (2026). *FIN Research Programs P281–P294: Regularized Spectral
Inference, Resource-Explicit Physics, and the Pointed-Torsor Boundary*
(Release 10.26; Version 1.0.0) [Preprint]. Zenodo.
