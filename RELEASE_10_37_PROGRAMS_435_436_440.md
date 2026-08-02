# Release 10.37 -- Exact Operational Certificates for FIN Programs P435, P436, and P440

## Process-tester SDP admission, a rigorously certified heralded-code gain, and detector-envelope minimax Jordan sampling

### Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher -- Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication -- Preprint
- **Version:** 1.0.0
- **Release:** 10.37
- **Publication date:** 2026-08-01
- **Language:** English
- **License:** CC BY 4.0
- **Access:** Open access

## Abstract

Release 10.37 advances three local operational research lanes in the FIN repository without using laboratory data, external audits, network resources, or remote computation.

Program P435 formulates finite phase/dephasing discrimination as a process-tester semidefinite program. At coherence (q=4/5) and phase (	heta=\pi/8), it supplies an exact matching one-slot primal and dual certificate with success probability

$$
\frac{1+(4/5)\sin(\pi/8)}2.
$$

It also exports the genuine two-slot memoryless comb and its causal primal/dual constraints. The multiround optimum remains open because no suitable SDP solver or independent checker is installed locally; this is reported as a reproducibility boundary, not a no-go theorem.

Program P436 converts one previously numerical heralded-erasure code advantage into a computer-assisted theorem. For (n=3), coherence and survival (q=\eta=4/5), phase (2\pi/15), and the exact symmetric sector law

$$
\frac1{250000}(56923,68077,68077,56923),
$$

the certified trace-distance gain over the better of product and GHZ baselines is at least

$$
0.022572776021405654.
$$

The proof combines Machin and alternating-series rational enclosures, exact heralded partial traces, exact characteristic-polynomial root isolation for rational Gram matrices, and an explicit nuclear-norm perturbation bound.

Program P440 derives the detector-envelope minimax extension of the variance-optimal Jordan sampler. For supplied efficiency lower bounds and dark-count upper bounds, the exact conditional allocation is

$$
q_i^\star\propto
|w_i|\|f_i\|_2
\sqrt{\frac1{\underline\epsilon_i}
+\frac{\overline d_i(1-\overline d_i)}{\underline\epsilon_i^2}}.
$$

On the certified P429 seven-atom midpoint, the declared worst mean-square-error coefficient is (7.4591\%) below absolute-weight sampling and (0.3401\%) below the detector-blind P422 law. A conservative simultaneous twelve-coordinate Hoeffding ledger is included.

## Principal scientific results

1. **O158 -- Process-Tester SDP Admission Certificate.** Exact for one slot; a complete two-slot instance is exported, while its matching local dual remains open.
2. **O159 -- Rational Heralded-Code Gain Certificate.** A strict, interval-certified gain over product and GHZ baselines in one nonideal three-use cell.
3. **O160 -- Detector-Envelope Minimax Jordan Sampler.** An exact conditional theorem for heterogeneous efficiency and independently subtracted dark counts.

## Scientific boundary

The release preserves the distinction between the intermediate legacy kernel and the later strict kernel. The reconstructed research-backed

$$
K^*_{\mathrm{legacy}}(d)
=A\frac{\cos(\pi d/4+\pi/6)}{1+\beta d}
$$

remains a fixed-phase legacy subclass, not the strict kernel. No physical-role transfer, complete legacy-to-strict bridge, typed self-coupling completion, non-premise selector, canonical orientation, dimensional unit, laboratory realization, Standard Model, gravity theory, total action, or Theory-of-Everything closure is claimed.

The detector envelope and all generated records are mathematical design inputs. They are not apparatus calibration or empirical evidence.

## Included files

- `FIN_Local_Research_Checkpoint_P435_P440_EN.pdf`
- `FIN_Local_Research_Checkpoint_P435_P440_EN.md`
- `FIN_Local_Research_Checkpoint_P435_P440_EN.tex`
- `fin_programs_435_436_440.py`
- `fin_programs_435_436_440_to_latex.py`
- `test_fin_programs_435_436_440.py`
- `test_fin_checkpoint_p435_p440.py`
- `FIN_Programs_435_436_440_Results.json`
- `FIN_Programs_435_436_440_Summary.csv`
- `FIN_Programs_435_436_440_P435_Comb_Certificate.csv`
- `FIN_Programs_435_436_440_P435_Two_Slot_Comb_Instance.npz`
- `FIN_Programs_435_436_440_P436_Erasure_Intervals.csv`
- `FIN_Programs_435_436_440_P440_Detector_Sampler.csv`
- `FIN_Programs_435_436_440_Figures/p435_p436_p440_certificates.png`
- `FIN_Local_Research_Checkpoint_P435_P440_State.json`
- `FIN_Checkpoint_P435_P440_AGENTS_Guardrail.txt`
- `FIN_PROGRAMS_435_436_440_RELEASE_MANIFEST.sha256`

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-mpl-1037 python3 fin_programs_435_436_440.py
python3 fin_programs_435_436_440_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P435_P440_EN.tex
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P435_P440_EN.tex
python3 -m unittest -v test_fin_programs_435_436_440.py test_fin_checkpoint_p435_p440.py
sha256sum -c FIN_PROGRAMS_435_436_440_RELEASE_MANIFEST.sha256
```

## Recommended continuation

The preferred next bounded batch is:

1. **P445:** dependency-free symmetry reduction or analytic dual for the two-slot qubit comb;
2. **P446:** interval branch-and-bound for global optimality in the three-use symmetric heralded-code simplex;
3. **P447:** propagation of the complete P429 atom boxes through the detector-envelope minimax allocation.

## Keywords

FIN; quantum channel discrimination; quantum comb; process tester; semidefinite programming; heralded erasure; dephasing; trace norm; interval arithmetic; exact root isolation; Jordan sampling; importance sampling; detector efficiency; dark counts; finite-sample confidence; mathematical physics.
