# Release 10.21 — FIN Laboratory Transfer Package

## P240 Optimal Tomography + P241 Blind-Custody Validator + P242 One-Shot Pipeline

### A scientific transfer package and executable specification for theoretical and experimental physicists

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint and executable specification
- **Version:** 1.0.0
- **Publication date:** 2026-07-27
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory
- **Related project DOI:** https://doi.org/10.5281/zenodo.21435332

---

## Abstract

This release transfers the experimentally actionable part of the FIN research programme to an independent team of theoretical and experimental physicists. It begins from the finite strict twelve-node kernel $K_{\mathrm{strict}}$ on $C_{12}$, not from cosmological or Theory-of-Everything claims. With $W = K_{\mathrm{strict}}$, $s \approx 1.660307$, and $A = sI - W$, the exact Dirichlet identity 

$$\langle f, Af \rangle = \frac{1}{2} \sum_{x,y} W_{xy} |f_x - f_y|^2 \geq 0$$

places the same self-adjoint generator under both the coherent group $U_t = e^{-itA}$ and the reversible heat/Markov semigroup $P_t = e^{-tA}$.

**Program 240** delivers the optimal spectral tomography protocol. The inverse-log noise factor $e^x/x$ yields the unique time minimizer at $x = \tau \lambda_{\max} = 1$. Equal allocation across all twelve basis preparations is proved minimax. Matrix-Bernstein spectral bounds supply a distribution-free concentration guarantee, accompanied by synthetic shot-plan analyses ($50\,000$ shots per preparation).

**Program 241** supplies an executable blind-custody validator checking eleven mandatory measurement, data, calibration, control, chronology, hash, and role-separation fields. It enforces strict provider–registrar–analyst role separation with detached GPG signatures.

**Program 242** is a fail-closed one-shot analysis pipeline. It evaluates holdout prediction $P_{2\tau}^{\mathrm{pred}} = \widehat{P}_\tau^2$ against the sealed $\widehat{P}_{2\tau}$ data and reports binary pass/fail without model repair. No external physical bundle is supplied in this release, keeping Program 242 locked until genuine laboratory data are admitted.

---

## Key Features & Main Results

1. **Dual Semantics of One Generator:** Proves that coherent quantum walk dynamics $U_t = e^{-itA}$ and diffusive Markov heat flow $P_t = e^{-tA}$ share the same spatial generator $A = sI - W$ on $C_{12}$, while clarifying that temporal semantics require operational instrument definition.
2. **Exact Tomographic Time Minimizer (P240):** Proves $x = \tau \lambda_{\max} = 1$ is the unique minimizer of local absolute-noise amplification $e^x/x$ for matrix logarithm rate extraction.
3. **Minimax Equal Allocation (P240):** Proves cyclic symmetrization makes equal shot allocation across all 12 basis preparations minimax within the declared convex design class.
4. **Production Custody Validator (P241):** Provides an 11-field schema and hash validator (`fin_lab_p241_validator.py`) enforcing provider $\neq$ registrar $\neq$ analyst role separation.
5. **Fail-Closed One-Shot Pipeline (P242):** Implements an atomic ledger-locked analysis pipeline (`fin_lab_p242_pipeline.py`) that evaluates time-homogeneity $P_{2\tau} = P_\tau^2$ and the $0.02$ projective spectral fingerprint without post-hoc threshold adjustment or model repair.
6. **Fully Validated Test Suite:** Includes an 18-test regression suite (`test_fin_lab_p240_242.py`), all passing cleanly.

---

## Included Files in `FIN_Laboratory_Transfer_Package_P240_242.zip`

- `FIN_Laboratory_Transfer_Package_P240_242.pdf` — Archival English preprint and executable laboratory transfer specification.
- `fin_lab_p240_optimal_tomography.py` — Target construction, time/allocation optimization, matrix Bernstein bounds, and synthetic shot planner.
- `FIN_Lab_P240_Design_Lock.json` — Canonical registered design lock for Program 240.
- `FIN_Lab_P240_Optimal_Tomography_Results.json` — Executed numerical results for P240.
- `fin_lab_p241_validator.py` — Production 11-field schema, hash, chronology, and role-separation validator.
- `FIN_Lab_P241_Transfer_Template/` — Complete directory of raw data schemas, manifest templates, and calibration headers for laboratories.
- `fin_lab_p242_pipeline.py` — Fail-closed one-run external analysis pipeline.
- `FIN_Lab_P242_Analysis_Lock.json` — Canonical frozen analysis lock for Program 242.
- `test_fin_lab_p240_242.py` — Complete test suite (18 unit/falsification tests).
- `FIN_Lab_P240_242_MANIFEST.sha256` — SHA-256 integrity manifest of all package files.

---

## Keywords

quantum process tomography; Markov semigroups; continuous-time quantum walk; matrix logarithm; matrix Bernstein bounds; blind custody; preregistration; fail-closed pipeline; operational physics; spectral graph theory; experimental quantum optics; integrated photonics.

---

## Suggested Citation

Żuchowski, K. (2026). *FIN Laboratory Transfer Package — Release 10.21: P240 Optimal Tomography + P241 Blind-Custody Validator + P242 One-Shot Pipeline* (FIN Research Monograph & Transfer Package; Version 1.0.0) [Preprint]. Zenodo. https://doi.org/10.5281/zenodo.21435332
