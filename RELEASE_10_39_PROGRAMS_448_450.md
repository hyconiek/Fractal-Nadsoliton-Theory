# Release 10.39 -- Global Majorization, a Causal Echo Counterexample, and a Representation-Gauge No-Go

## FIN Programs P448, P449, and P450

### Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher -- Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication -- Preprint
- **Version:** 1.0.0
- **Release:** 10.39
- **Publication date:** 2026-08-01
- **Language:** English
- **License:** CC BY 4.0
- **Access:** Open access

## Abstract

Release 10.39 reports three local analytical and computer-assisted FIN research programs. It uses no laboratory data, external audit, network resource, remote computation, or physical validation.

Program P448 constructs a concave fine-grained-erasure majorant for the three-use heralded discrimination objective. Trace-distance contractivity, a nuclear-norm variational formula, reversal symmetry, and outward-rounded interval branch-and-bound give the first rigorous upper bound on the complete Hamming-sector simplex:

$$
0.46327828319203235
\le D_{\mathrm{global}}
\le0.4666305033804779.
$$

The remaining gap is 0.003352220188445554. Exact globality remains open because the majorant discloses the computational values of lost bits.

Program P449 proves the recursive causal-support cone

$$
\mathcal C_n
=\{B_0\oplus B_1:B_0,B_1\succeq0,
B_0+B_1\in\mathcal C_{n-1}\}
$$

and the affine-dimension recurrence $d_1=1$, $d_n=d_{n-1}+4^{n-1}$. It then constructs a full positive three-slot tester extending the rational echo-history law

$$
t=\frac1{40}(9,1,9,1,1,9,1,9).
$$

Rational transcendental enclosures and exact spectral root isolation certify

$$
D_{\mathrm{echo}}
\ge0.49063899018433244
>0.47302632064577882
=D_{\mathrm{GHZ},3}.
$$

The advantage is at least 0.017612669538553616. Thus three-slot GHZ optimality is refuted in the declared reduced phase/dephasing comb. The global three-slot optimum remains open.

Program P450 proves that the detector-aware sampler O163 is not determined by the twelve target moments. The exact finite-difference measure

$$
\nu_{12}=\sum_{j=0}^{12}(-1)^j\binom{12}{j}\delta_{j/12}
$$

annihilates every polynomial of degree at most eleven. Adding a multiple of this null cycle preserves all target moments but makes the declared minimax risk unbounded. A preparation-level representation rule or quotient-invariant estimator is therefore necessary.

## Principal scientific results

1. **O164 -- Concave Fine-Grained-Erasure Majorant.** A rigorous full-simplex upper certificate, with the coarse objective's exact globality explicitly open.
2. **O165 -- Recursive Causal-Support Cone.** Exact recursion and dimension theorem; an explicit full-space echo-history tester strictly refutes three-slot GHZ optimality.
3. **O166 -- Moment-Null Representation Gauge.** Exact representation-dependence and risk-unboundedness theorem for O163.

## Scientific boundary

The programs are downstream of supplied reduced channels or moment representations. They do not derive either FIN kernel, type the historical coupling diagram as a self-evolution, complete the legacy-to-strict bridge, or transfer physical roles. The echo law is a mathematical tester, not a laboratory instrument. No selector, canonical orientation, dimensional unit, detector calibration, physical record, Standard Model, gravity theory, total action, or Theory-of-Everything closure is claimed.

## Included files

- `FIN_Local_Research_Checkpoint_P448_P450_EN.pdf`
- `FIN_Local_Research_Checkpoint_P448_P450_EN.md`
- `FIN_Local_Research_Checkpoint_P448_P450_EN.tex`
- `fin_programs_448_450.py`
- `fin_programs_448_450_to_latex.py`
- `test_fin_programs_448_450.py`
- `test_fin_checkpoint_p448_p450.py`
- `FIN_Programs_448_450_Results.json`
- `FIN_Programs_448_450_Summary.csv`
- `FIN_Programs_448_450_P448_Full_Simplex_Majorant.csv`
- `FIN_Programs_448_450_P449_Three_Slot_Recursion.csv`
- `FIN_Programs_448_450_P449_Three_Slot_Witness.npz`
- `FIN_Programs_448_450_P450_Null_Cycle.csv`
- `FIN_Programs_448_450_Figures/p448_p450_global_majorant_and_gauge_obstruction.png`
- `FIN_Local_Research_Checkpoint_P448_P450_State.json`
- `FIN_Checkpoint_P448_P450_AGENTS_Guardrail.txt`
- `FIN_PROGRAMS_448_450_RELEASE_MANIFEST.sha256`

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-mpl-1038 python3 fin_programs_448_450.py
python3 fin_programs_448_450_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P448_P450_EN.tex
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P448_P450_EN.tex
MPLCONFIGDIR=/tmp/fin-mpl-1038 python3 -m unittest -v test_fin_programs_448_450.py test_fin_checkpoint_p448_p450.py
sha256sum -c FIN_PROGRAMS_448_450_RELEASE_MANIFEST.sha256
```

## Recommended continuation

The preferred next bounded batch is:

1. **P451:** globally optimize the 21-dimensional three-slot support cone or certify a sharper echo witness;
2. **P452:** close or narrow the P448 coarse-graining gap by a direct full-simplex certificate;
3. **P453:** prove uniqueness of a canonical reduced/minimum-total-variation atom representation or prove a no-go.

## Keywords

FIN; quantum comb; causal tester; channel discrimination; echo history; GHZ; heralded erasure; trace norm; concavity; interval arithmetic; signed moment problem; representation gauge; mathematical physics.
