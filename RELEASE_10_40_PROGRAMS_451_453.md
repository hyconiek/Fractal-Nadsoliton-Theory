# Release 10.40 — FIN Local Research Programs P451–P453

## Coherence Advantage, Global Coarse-Erasure Symmetry, and Canonical Jordan Representation

**Creator:** Żuchowski, Krzysztof  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Resource type:** Publication — Preprint  
**Version:** 1.0.0  
**Publication date:** 2026-08-01  
**Language:** English  
**License:** CC BY 4.0

## Abstract

Release 10.40 records three local-only analytical and computer-assisted FIN
programs. P451 solves the complete diagonal face of a declared three-slot
reduced causal tester problem and supplies an exact-rational coherent tester
whose certified half-distance is strictly larger. P452 proves the missing
coarse-erasure reversal inequality and promotes the earlier palindromic-line
interval result to a full four-sector simplex certificate. P453 combines the
P429 Krawczyk solution and P430 global dual with strict complementarity and
Vandermonde injectivity to prove global uniqueness of the minimum-negative-
mass signed measure.

## Principal results

- Global diagonal three-slot half-distance:
  `0.5227633693019748...`.
- Certified coherent witness interval:
  `[0.52332810026048937, 0.52332810026049204]`.
- Certified coherence advantage over the diagonal upper bound:
  at least `0.00056473095851463685`.
- Exact coarse-erasure symmetrization condition:
  `4 q^2 cos^2(theta) >= 1`; certified margin at the declared P446
  parameters: `1.1364871761393385...`.
- Full-simplex coarse-objective bracket:
  `[0.46327828319203235, 0.4642707484333634]`.
- Unique global minimum-negative-mass / fixed-mass minimum-TV signed measure,
  with certified Vandermonde lower bound
  `1.4954513576571674e-9`.

## New typed objects

- O167 — Causal-Coherence Advantage Witness.
- O168 — Coarse-Erasure Symmetrization Certificate.
- O169 — Strict-Complementarity Jordan Gauge Fixer.

## Scope and nonclaims

All work is local and dimensionless. No laboratory record, external audit,
internet result, remote computation, or physical evidence is used. The full
21-dimensional P451 optimum remains open. The minimum-TV rule is a
mathematical variational rule, not a physical preparation theorem. The
legacy and strict kernels remain separated; no legacy physical role is
transferred. QW-2191, canonical orientation, dimensional units, apparatus,
complete legacy-to-strict bridge, Standard Model, gravity, L_total, and ToE
closure remain open or unclaimed.

## Main document

`FIN_Local_Research_Checkpoint_P451_P453_EN.pdf`

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-mpl-1038 python3 fin_programs_451_453.py
python3 fin_programs_451_453_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P451_P453_EN.tex
python3 -m unittest -v test_fin_programs_451_453.py test_fin_checkpoint_p451_p453.py
sha256sum -c FIN_PROGRAMS_451_453_RELEASE_MANIFEST.sha256
```

## Next selected batch

P454 — full-cone primal/dual certificate; P455 — exact coherent symmetry-face
optimization; P457 — sharpen the licensed P452 one-dimensional certificate.
