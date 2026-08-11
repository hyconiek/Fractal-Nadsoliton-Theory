# Release 10.68 — FIN ST166–ST177: Carrier Invariants, Exact Uniform Folds, and the Executable Physical Boundary

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint and reproducibility package
- **Version:** 1.0.0
- **Publication date:** 2026-08-11
- **Language:** English
- **License:** CC BY 4.0
- **Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Abstract

Release 10.68 executes FIN research programs ST166–ST177. It asks whether FIN information has a rigorous carrier-independent content, while continuing exact fold classification, recovery, selector geometry, nonlinear-source analysis, interval continuation, noisy refinement, robust tomography, hidden-process inference, bath construction, and operational validation.

ST166 proves the strongest presently justified carrier statement. Let a logical state space be embedded in different stochastic carriers by maps (E_i) with stochastic left inverses (D_i), and transport every logical effect (m) as (mD_i). Then the entire preparation–effect probability table is invariant:

\[
mD_iE_ip=mp.
\]

Two inequivalent finite carriers verify the identity exactly. A same-Shannon-entropy counterexample proves that Shannon entropy alone is not a complete operational invariant. This is a theorem about representation equivalence, not evidence that physical information propagates without any carrier.

ST167 proves that environment couplings restricted to (L_k=f_k(A)) cannot generate a finer vertex pointer algebra. Their joint commutant contains the 22-dimensional strict commutant. A state-dependent, nonlinear, vertex-typed, or symmetry-breaking coupling would therefore be a genuinely new object.

ST168 completely classifies folds on the uniform reflection-even branch. Besides (c=0), stationary amplitudes satisfy

\[
3\rho^3+18\rho^2-284\rho+24=0,
\qquad \rho=c^2.
\]

The cubic has exactly two positive roots. An exact nonzero resultant excludes a constant-mode fold. Five real amplitudes and six positive reflection-even Fourier classes produce exactly 30 geometric fold null lines, or 60 normalized augmented roots after eigenvector sign. This is complete only on the declared symmetry-reduced branch.

ST169 solves optimal recovery for two nonorthogonal qubit phase errors with noisy syndrome. ST170 proves that a symmetry-preserving semigroup cannot start a selector from the uniform state and gives the sharp conditional Fisher–Rao action to reach a supplied probability gap. ST171 proves that the degree-twelve anisotropy coupling is invisible to the strict operator, functional calculus, symmetry, and every field derivative through order eleven.

ST172 certifies all 125 cells of a nuisance cube of halfwidth (7.5\times10^{-5}). ST173 gives exact noisy-refinement parity and majority-decoding laws. ST174 proves a robust four-dimensional commutant subspace under the analytical threshold (16\eta+8\eta^2<2). ST175 sharpens the hidden-HMM pair-transfer spectral-radius bracket to

\[
[0.7211542214950951,\;0.7211542214950967].
\]

ST176 separates the minimal four-qubit bath dimension, depth-one register-SWAP circuit, norm-budget time, and eight-body manifest energy-conserving generator. ST177 implements an eleven-atom operational validator. The local mathematical packet passes structurally, while physical execution correctly remains blocked.

## Program outcomes

- **ST166 — Proven for reversible stochastic code embeddings:** the complete operational probability table is carrier invariant when instruments are transported. Shannon entropy alone is insufficient.
- **ST167 — Proven no-go in strict functional calculus:** (f(A))-couplings cannot generate the 12-dimensional vertex pointer algebra.
- **ST168 — Proven complete on the uniform reflection-even branch:** 30 geometric folds and 60 normalized signed augmented roots.
- **ST169 — Proven:** exact optimal recovery for binary nonorthogonal qubit phase errors with noisy syndrome.
- **ST170 — Proven/conditional:** covariant no-start theorem and sharp Fisher–Rao selector-control cost for supplied control data.
- **ST171 — Proven:** a degree-twelve coupling cannot be inferred from strict linear and order-11 field data.
- **ST172 — Proven:** all 125 interval Krawczyk cells pass over the enlarged nuisance cube.
- **ST173 — Proven for the declared noise model:** exact path-parity and redundant-majority laws.
- **ST174 — Proven analytically with numerical audit:** a unique four-dimensional low commutant subspace is separated under the stated noise threshold.
- **ST175 — Proven:** outward Collatz–Wielandt spectral-radius bracket and improved pair-path exponent.
- **ST176 — Proven construction/tradeoff:** minimal bath dimension and explicit locality, time, and energy-conservation accounting.
- **ST177 — Proven executable boundary:** schema, hashes, role declaration, and all eleven deletion failures are checked; physical execution is false.

## Central conclusion

FIN now has a precise notion of information that can remain invariant under a change of carrier: an equivalence class of tomographically complete preparation–transformation–measurement statistics. This is stronger than equality of Shannon entropy and weaker than a claim of physical existence without a carrier.

The strict operator remains insufficient to source its own fine pointer algebra. Symmetry supplies equivalent possibilities but is not positive feedback by itself; a selector requires instability, state dependence, control, boundary data, or another new source. The exact uniform fold atlas supplies rigorous seeds for a wider validated continuation, but the degree-twelve nonlinear coefficient is independent of all lower-order strict data.

No strict selector, calibrated dimension, laboratory apparatus or record, physical projection theorem, legacy-to-strict bridge completion, legacy physical-role transfer, Standard Model, gravity, (L_{\mathrm{total}}), or Theory-of-Everything closure is claimed.

## Included artifacts

- `FIN_ST166_ST177_Carrier_Invariants_Exact_Folds_and_Executable_Physical_Boundary_Report_EN.pdf`
- `FIN_ST166_ST177_Carrier_Invariants_Exact_Folds_and_Executable_Physical_Boundary_Report_EN.tex`
- `fin_st166_st177_research.py`
- `fin_st177_operational_validator.py`
- `test_fin_st166_st177.py`
- `FIN_ST166_ST177_Results.json`
- `FIN_ST166_ST177_Summary.csv`
- `FIN_ST166_Carrier_Operational_Invariant.json`
- `FIN_ST168_Uniform_Branch_Fold_Classification.json`
- `FIN_ST169_Nonorthogonal_Unitary_Recovery.json`
- `FIN_ST172_Extended_Nuisance_Cover.json`
- `FIN_ST173_Noisy_Refinement_Compression.json`
- `FIN_ST174_Robust_Tomographic_Net.json`
- `FIN_ST175_HMM_Collatz_Certificate.json`
- `FIN_ST176_Microscopic_Local_Bath_Tradeoff.json`
- `FIN_ST177_Operational_Validator_Study.json`
- `FIN_ST177_Operational_Validator_Record.json`
- `FIN_ST166_ST177_Figures/` — nine generated figures
- `FIN_ST166_ST177_SHA256SUMS.txt`

## Validation

The complete local suite passes **28/28 tests**. Live checks cover reversible stochastic encodings, exact operational tables, entropy counterexamples, functional-calculus commutants, exact cubic isolation and resultant, uniform-fold counts, binary recovery formulas, Fisher–Rao costs, the order-11 source obstruction, all 125 nuisance cells, path parity and majority decoding, commutator-Laplacian spectrum and perturbation threshold, interval Collatz bounds, bath dimensions and circuit depth, operational hashes, role separation, and all eleven deletion witnesses.

## Recommended continuation

Programs ST178–ST189 are specified in the report. Highest priority is ST178: replace exact reversible carrier embeddings by irreversible encodings and characterize the sufficient statistic that survives. ST179 should construct and attempt to falsify one strict nonlinear or state-dependent environment coupling outside (f(A)). ST180 should launch validated nonuniform continuation from the 30 exact uniform folds. ST181–ST189 extend nonorthogonal recovery, controlled selector dynamics, twelfth-order source audits, nuisance continuation, correlated refinement noise, factor-algebra recovery, exact hidden-HMM exponents, energy-conserving local resets, and an external immutable-record interface.

## Keywords

FIN; operational information; carrier invariance; stochastic embedding; commutant; functional calculus; fold bifurcation; exact root isolation; quantum recovery; Fisher–Rao geometry; nonlinear identifiability; interval Krawczyk method; noisy refinement; commutator Laplacian; hidden Markov model; Collatz–Wielandt theorem; microscopic bath; operational validator.
