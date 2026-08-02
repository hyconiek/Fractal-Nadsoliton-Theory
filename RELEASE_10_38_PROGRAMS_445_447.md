# Release 10.38 -- Structural Closure and Certified Propagation for FIN Programs P445--P447

## Exact two-slot comb optimality, a scoped global heralded-code certificate, and full atom-box detector allocation

### Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher -- Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication -- Preprint
- **Version:** 1.0.0
- **Release:** 10.38
- **Publication date:** 2026-08-01
- **Language:** English
- **License:** CC BY 4.0
- **Access:** Open access

## Abstract

Release 10.38 reports three local analytical and computer-assisted FIN research programs. It uses no laboratory record, external audit, network resource, remote computation, or physical validation.

Program P445 replaces the previously open two-slot semidefinite-program cell by an exact structural theorem. The four-dimensional process support of the declared reduced qubit phase/dephasing channel compresses every causal tester normalizer to four history weights and one constrained block coherence. The complete trace-norm objective is

$$
D^2=a^2(t_{00}+t_{11})(t_{01}+t_{10})+b^2t_{00}t_{11}
-4a^2u^2+2ab(t_{11}-t_{00})u,
$$

where $a=2q\sin\theta$, $b=2q^2\sin(2\theta)$, and $u$ is the real part of the causal coherence. Completing a square proves that, for $q=4/5$ and $\theta=\pi/8$, the GHZ history law is globally optimal over every admissible two-slot tester, including arbitrary intermediate memory. The exact half-distance and equal-prior success are

$$
D_2^\star=\frac{8\sqrt2}{25},\qquad
p_{\mathrm{succ},2}^\star=\frac12+\frac{4\sqrt2}{25}.
$$

Program P446 performs an $S_3$ irreducible-block reduction of the three-use heralded-erasure code problem. Outward-rounded interval branch-and-bound proves a global palindromic-line upper/lower gap of $0.0009924652413310642$. A 24-start full-simplex audit and 16 sampled interior Hessians provide strong evidence for the same reversal-symmetric optimum, but do not prove full-simplex concavity or globality.

Program P447 propagates the complete exact P429 rational atom boxes through the conditional detector-envelope minimax allocation. The maximum probability-interval width is $3.2547095453761505\times10^{-20}$, and the strict worst-MSE ordering

$$
\mathrm{MSE}_{\mathrm{P447}}<\mathrm{MSE}_{\mathrm{P422}}<\mathrm{MSE}_{|w|}
$$

is certified throughout the box. The detector envelope remains a supplied synthetic model, not calibration or empirical evidence.

## Principal scientific results

1. **O161 -- Diagonal-Support Comb Reduction.** An exact characterization and optimization theorem for every causal two-slot tester of the declared commuting phase/dephasing channel.
2. **O162 -- Palindromic Heralded-Code Global Upper Certificate.** A rigorous one-dimensional global certificate, with the three-dimensional full simplex explicitly left open.
3. **O163 -- Certified Detector-Allocation Tube.** Exact propagation of all P429 atom uncertainty through the conditional detector minimax law.

## Scientific boundary

The reconstructed historical kernel

$$
K^*_{\mathrm{legacy}}(d)
=A\frac{\cos(\pi d/4+\pi/6)}{1+\beta d}
$$

remains a fixed-phase subclass of the intermediate legacy kernel and is not the strict kernel. P445--P447 neither derive a kernel nor complete the legacy-to-strict bridge. They do not type the coupling intuition in `DIAGRAMS_KERNEL_TRANSFORMATION.md` as a self-evolution, transfer any legacy physical role, source a selector or orientation, generate dimensional units, or provide laboratory evidence. No Standard Model, gravity theory, total action, or Theory-of-Everything closure is claimed.

## Included files

- `FIN_Local_Research_Checkpoint_P445_P447_EN.pdf`
- `FIN_Local_Research_Checkpoint_P445_P447_EN.md`
- `FIN_Local_Research_Checkpoint_P445_P447_EN.tex`
- `fin_programs_445_447.py`
- `fin_programs_445_447_to_latex.py`
- `test_fin_programs_445_447.py`
- `test_fin_checkpoint_p445_p447.py`
- `FIN_Programs_445_447_Results.json`
- `FIN_Programs_445_447_Summary.csv`
- `FIN_Programs_445_447_P445_Two_Slot_Reduction.csv`
- `FIN_Programs_445_447_P446_Palindromic_Certificate.csv`
- `FIN_Programs_445_447_P446_Full_Simplex_Audit.csv`
- `FIN_Programs_445_447_P447_Atom_Box_Propagation.csv`
- `FIN_Programs_445_447_Figures/p445_p447_exact_reduction_and_intervals.png`
- `FIN_Local_Research_Checkpoint_P445_P447_State.json`
- `FIN_Checkpoint_P445_P447_AGENTS_Guardrail.txt`
- `FIN_PROGRAMS_445_447_RELEASE_MANIFEST.sha256`

## Reproduction

```bash
MPLCONFIGDIR=/tmp/fin-mpl-1038 python3 fin_programs_445_447.py
python3 fin_programs_445_447_to_latex.py
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P445_P447_EN.tex
lualatex -interaction=nonstopmode -halt-on-error FIN_Local_Research_Checkpoint_P445_P447_EN.tex
python3 -m unittest -v test_fin_programs_445_447.py test_fin_checkpoint_p445_p447.py
sha256sum -c FIN_PROGRAMS_445_447_RELEASE_MANIFEST.sha256
```

## Recommended continuation

The preferred next bounded batch is:

1. **P448:** prove full-simplex concavity for the P446 code objective or construct a counterexample;
2. **P449:** extend O161 to three slots or isolate the first genuinely new causal coherence obstruction;
3. **P450:** determine whether O163 is invariant under every admitted signed-atom representation or representation-dependent.

## Keywords

FIN; quantum comb; process tester; channel discrimination; dephasing; GHZ; heralded erasure; representation theory; interval branch-and-bound; trace norm; detector allocation; interval propagation; mathematical physics.
