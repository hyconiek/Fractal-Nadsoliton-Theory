# Release 10.27 — FIN Research Programs P295–P308

## Legacy physics revisited: positive path mixtures, strict spectral completion, and physical-role transfer

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 1.0.0
- **Publication date:** 2026-07-30
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0
- **Repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory

## Abstract

Release 10.27 executes fourteen research programs, P295–P308, that compare
the historical legacy kernel and its proposed physical roles with the current
strict finite spectral operator.

The principal analytic result is a sharp separation of attenuation
mechanisms. The legacy envelope

\[
\frac{1}{1+\beta d}
\]

is exactly a positive Laplace mixture of exponential distance attenuations.
The strict envelope

\[
\frac{1}{1+d^{1.8}}
\]

is not completely monotone near the origin and cannot belong to the same
positive mixture class. Together with four phase-sign mismatches, forced
flat-band multiplicities, and 32 spectral mode-order inversions, this excludes
positive fixed-phase completion and monotone Bernstein functional calculus
as full legacy-to-strict bridges.

The release also audits historical maps to the electroweak angle,
electromagnetic coupling, gravitational hierarchy, and dimensional scales.
Those maps are conditional parameter formulas, not proved functionals of the
legacy or strict operator. The electromagnetic and gravity formulas fail a
distance-coordinate gauge test; none has a complete strict role-transfer
certificate.

Constructive results include multi-probe resolvent recovery, a
machine-checked staged-minimization theorem, exact complex-unitary
compilation of the strict current POVM, a multitime memory witness, a
full-rank synthetic nuisance-detector design, local strict-fingerprint coarea
asymptotics, an adversarial finite-time clock no-go, intervention-assisted
synthetic coordinate recovery, a scale-orbit quotient, and a pointed
role-transfer theorem.

The release constructs the **Typed Legacy–Strict Completion Span** and the
**Role-Transfer Certificate**. It does not claim full bridge closure, physical
confirmation, a selector, a dimensional constant, Standard Model or
gravitational dynamics, \(L_{\rm total}\), or a Theory of Everything.

## Principal results

1. **P295:** four shared Stieltjes probes increase local inverse rank from
   \(6\) to \(15\) and improve frozen pole recovery, without providing
   uniform stability.
2. **P296:** Lean 4.28 checks general staged-minimization logic; 400 exact
   rational SPD Schur witnesses pass.
3. **P297:** the six-effect strict current POVM compiles into an exact
   \(72\times72\) complex-unitary Givens mesh with unitarity residual
   \(2.12\times10^{-15}\).
4. **P298:** scalar complement closure is obstructed by forced flat bands,
   and monotone Bernstein completion is obstructed by 32 mode-order
   inversions.
5. **P299:** a supplied hidden-rate process has positive conditional mutual
   information and differs from its best one-step Markov fit by TV
   \(0.02157\).
6. **P300:** the historical legacy envelope has an exact positive path-scale
   measure; strict lies outside that cone, and unchanged legacy phase fails
   at four nodes.
7. **P301:** no independently admissible P241 event bundle was accepted; P242
   was not run on fabricated evidence.
8. **P302:** the frozen synthetic five-nuisance detector design is full rank
   and gives gradient RMSE \(0.00906\).
9. **P303:** the strict fingerprint chart has rank five and local small-ball
   exponent five; importance sampling gives \(5.217\).
10. **P304:** a positive nonlinear clock warp exactly matches all three
    observed propagators while changing off-grid predictions.
11. **P305:** interventions raise a supplied bridge-coordinate library from
    rank four to six and recover its synthetic law on hold-out.
12. **P306:** no independently certified physical reservoir record is
    present.
13. **P307:** the strict heat/Gibbs family has a scale orbit; historical
    electromagnetic and gravity roles fail the coordinate-gauge filter.
14. **P308:** a free torsor has no invariant point; a supplied point makes an
    equivariant map unique but does not derive the point or transfer a
    physical role.

## Historical corrections

- \(\cos(\pi d/4+\pi/6)=0\) at \(d=4/3+4k\), not at the integer nodes
  stated in old diagrams.
- \(d^{1.6}d^{-0.6}=d^{+1}\), not \(d^{-1}\); the stated path-counting
  argument would require a \(d^{-2.6}\) path weight.
- The old electroweak-angle manipulation contains a cancellation yielding
  \(\alpha_{\rm geo}\), not \(\alpha_{\rm geo}/12\).
- Historical Wilson-loop plots are internal simulations, not independent
  laboratory measurements.
- The old physical formulas do not factor through a proved kernel/operator
  observable and cannot be silently transferred to strict.

## New theoretical objects

- **O54 — Spectral Observation Signature**
- **O55 — Bernstein Spectral Completion** (refuted in the tested class)
- **O56 — Operational Prediction Pseudometric**
- **O57 — Positive Path-Scale Measure**
- **O58 — Scale-Orbit Prediction Quotient**
- **O59 — Role-Transfer Certificate**
- **O60 — Pointed Operational Kernel Package**
- **O61 — Typed Legacy–Strict Completion Span**

## Included files

- `FIN_Programs_295_308_Legacy_Strict_Physics_Report_EN.pdf` — archival
  English monograph;
- `FIN_Programs_295_308_Legacy_Strict_Physics_Report_EN.md` — editable report;
- `FIN_Programs_295_308_Legacy_Strict_Physics_Report_EN.tex` — generated
  LaTeX source;
- `fin_programs_295_308.py` — deterministic research runner;
- `test_fin_programs_295_308.py` — 17-test verification suite;
- `FIN_Programs_295_308_Formal_Core.lean` — Lean 4.28 certificate;
- `FIN_Programs_295_308_Results.json` — machine-readable results;
- `FIN_Programs_295_308_Summary.csv` — compact result matrix;
- eight detailed CSV audit tables;
- five figures in `FIN_Programs_295_308_Figures/`;
- `FIN_Programs_295_308_SHA256SUMS.txt` — integrity manifest.

## Verification

```bash
MPLCONFIGDIR=/tmp/mpl-fin \
python3 -m unittest -v test_fin_programs_295_308.py
```

Expected result: **17 passed, 0 failed**. The suite recompiles the Lean
certificate.

## Recommended continuation

The report ranks fourteen programs P309–P322. The recommended immediate
sequence is

\[
\mathrm{P311}
\rightarrow
\mathrm{P316}
\rightarrow
\mathrm{P312}
\rightarrow
\mathrm{P315}
\rightarrow
\mathrm{P320},
\]

covering a loss-complete unitary measurement model, the minimum signed/complex
path resource, a nontrivial common-parent search, operational
legacy-versus-strict discrimination, and a clock-sensitive experiment.

External P241 and hardware-reservoir lanes remain blocked until genuine
providers, registrars, calibrated apparatus, event records, and frozen
hold-outs exist.

## Keywords

FIN; legacy kernel; strict kernel; spectral operator; positive Laplace
mixture; complete monotonicity; Bernstein function; Naimark dilation; complex
unitary mesh; Stieltjes inverse; Schur complement; process memory; detector
tomography; coarea formula; scale gauge; torsor; physical-role transfer;
operational physics; mathematical physics.

## Suggested citation

Żuchowski, K. (2026). *FIN Research Programs P295–P308: Legacy Physics
Revisited—Positive Path Mixtures, Strict Spectral Completion, and
Physical-Role Transfer* (Release 10.27; Version 1.0.0) [Preprint]. Zenodo.

