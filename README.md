# Fractal Information Nadsoliton (FIN)

**Author:** Krzysztof Żuchowski  
**Repository:** https://github.com/hyconiek/Fractal-Nadsoliton-Theory  
**DOI (Release 10.88 Compendium):** https://doi.org/10.5281/zenodo.17645737  
**Current Milestone:** Release 11.09+ / Programs ST2352–ST3251 (Current as of September 2026)  

---

## 1. Executive Summary & Epistemic Verdict (September 2026)

This repository documents the mathematical reconstruction, operational testing, and research frontier of the **Fractal Information Nadsoliton (FIN)** theory from **Release 10.88 upwards**. All obsolete early-2026 exploratory tracks have been superseded by the claim-graded compendium and subsequent program suites:

1. **Mathematical Rigor & Internal Closure:**  
   The finite, dimensionless spectral-information framework on the cyclic group \(\mathbb{Z}/12\mathbb{Z}\)—governed by the Dirichlet generator \(A = sI - W\), exact functional calculi, rational Krawczyk inclusions, blockwise algebras, and exact prism renormalization group equations—is **mathematically closed, verified by standard-library reproducible codes, and rigorously audited**.
2. **Full Fundamental Theory of Everything (ToE) Closure:**  
   **NOT YET.** The theory does not derive a continuum spacetime, an action-level Standard Model Lagrangian, Einstein tensor gravity, a microscopic physical time arrow, or canonical SI physical units (seconds, metres, kilograms).
3. **Empirical & Community Confirmation:**  
   **NOT YET.** Laboratory transfer protocols and operational process-tester architectures are formally defined, but held-out external multiteam experimental replication remains pending.

---

## 2. Founding Ontology & Architectural Axioms

The foundational ontology of FIN establishes a strict hierarchy:

$$\boxed{\text{Nadsoliton} \longrightarrow \text{Light} \longrightarrow \text{Matter} \longrightarrow \text{Emergent Observer}}$$

- **The Nadsoliton as Primordial Information:**  
  The nadsoliton itself is the primordial self-similar information of the universe in a persistent solitonic state. There is **no separate informational substrate underneath the nadsoliton**.
- **The Master Separation Theorem (Release 10.88):**  
  A single self-adjoint generator \(A\) generates multiple mathematical shadows (diffusive heat, unitary oscillations, wave dynamics, Green resolvents, quadratic actions, spectral filtering, and contextual Schur reductions). These constructions share spectral data, but **they do not automatically inherit physical reality without explicit coupling and measurement axioms**.
- **Historical Coupling Corrections:**  
  Release 10.88 formally audited the early coupling diagrams:
  - Multiplying an exponential attenuation after a hyperbolic path-sum transformation double-counts damping.
  - The algebraic relation \(d^{1.6} d^{-0.6} = d\) (not \(d^{-1}\)).
  - The integer phase nodes claim was refuted because \(\cos\left(\frac{\pi d}{4} + \frac{\pi}{6}\right) = 0 \iff d = \frac{4}{3} + 4n \notin \mathbb{Z}\).

---

## 3. Authoritative Project Guardrails (from `AGENTS.md`)

All theoretical and automated operations in this repository are strictly bound by the following non-negotiable guardrails:

### 3.1. Kernel Split & Bridge-Completion Guardrail
- **Legacy Kernel:**  
  $$K_{\text{legacy\_ont}}(d) = \alpha_{\text{geo}} \frac{\cos(\omega d + \phi)}{1 + \beta_{\text{tors}} d}$$  
  with \((\alpha_{\text{geo}}, \omega, \phi, \beta_{\text{tors}}) = (4\ln 2, \pi/4, \pi/6, 0.01)\), is restored solely as an **intermediate bridge kernel**. It omits strict-side characteristics of the nadsoliton (nonlinear compression, spatial damping exponent \(\eta\), and certified topological data).
- **Strict Working Kernel:**  
  $$K_{\text{strict\_gate}}(d) = \frac{\cos(\omega d + \phi)}{1 + \beta d^{\eta}}$$  
  with \((\omega, \phi, \beta, \eta) = (0.18575, 0.1625, 1.0, 1.8)\), is the completed, enriched strict working kernel.
- **No Silent Substitution & Role-Transfer Ban:**  
  Do not substitute \(K_{\text{strict\_gate}}\) for \(K_{\text{legacy\_ont}}\). Furthermore, legacy physical-role claims (\(\sin^2\theta_W = \alpha_{\text{geo}}/12\), \(\alpha_{\text{EM}}^{-1} = \frac{\alpha_{\text{geo}}}{2\beta_{\text{tors}}}(1-\beta_{\text{tors}})\), \(\beta^N\) gravity hierarchy) must **never** be transferred onto \(K_{\text{strict\_gate}}\) without an explicit bridge map **and** a separately proven role-transfer theorem.

### 3.2. Selector & Uniqueness Guardrail (`QW-2191`)
- Obstruction `QW-2191` is an active, proven strict-core uniqueness obstruction.
- No strict-core selector closure may be claimed without an explicit symmetry-breaking premise or a newly exported internal selector source. Axiom-augmented closures must remain explicitly labeled as non-strict.

### 3.3. Noncyclic & Strategic Steering Guardrail
- **L5/L12 Noncyclic Guardrail:** Repeated cyclic gate generation under identical blocker-cuts is prohibited (`QW-2381/QW-2382/QW-2383`). Loop expansion requires a noncyclic anchor or genuinely new provider class.
- **State-Map-First Steering:** The theoretical priority is governed by the live proof frontier rather than circular generic bridge repetitions.

---

## 4. Minimal Strict Spectral Core (12-Cycle Dirichlet Generator)

On the discrete 12-cycle \(V = \mathbb{Z}/12\mathbb{Z}\), with cyclic distance \(d(x, y) = \min(|x-y|, 12-|x-y|)\), the strict interaction matrix \(W\) is defined by:

$$W_{xx} = 0, \qquad W_{xy} = K_{\text{strict}}(d(x, y)) \quad (x \neq y)$$

The strict weights evaluate to:
- \(k_1 = K_{\text{strict}}(1) \approx 0.471676\)
- \(k_2 = K_{\text{strict}}(2) \approx 0.182414\)
- \(k_3 = K_{\text{strict}}(3) \approx 0.091176\)
- \(k_4 = K_{\text{strict}}(4) \approx 0.051515\)
- \(k_5 = K_{\text{strict}}(5) \approx 0.030588\)
- \(k_6 = K_{\text{strict}}(6) \approx 0.017700\)

Exact row sum:
$$s = 2 \sum_{d=1}^5 k_d + k_6 \approx 1.660307278766$$

The self-adjoint positive semi-definite Dirichlet generator is:
$$A = sI - W$$

Its eigenvalues \(\lambda_k = s - \hat{W}_k\) satisfy \(\lambda_0 = 0\) (ground state on the constant mode \(\mathbf{1}\)) and \(\lambda_k > 0\) for all \(k \in \{1, \dots, 11\}\), generating:
- **Unitary Dynamics:** \(U(t) = \exp(-iAt)\)
- **Heat Diffusion:** \(H(t) = \exp(-At)\)
- **Resolvent / Green Operator:** \(G(z) = (zI - A)^{-1}\) on the mean-zero subspace
- **Discrete Quadratic Action:** \(S[\psi] = \frac{1}{2} \langle \psi, A\psi \rangle = \frac{1}{4} \sum_{x, y} W_{xy} (\psi_x - \psi_y)^2\)

---

## 5. Major Research Releases (Release 10.88 through Release 11.09)

### Release 10.88 — The FIN Current-State Compendium (August 2026)
- Published as the authoritative state monograph: [`ZENODO_RELEASE_10_88_FIN_CURRENT_STATE_COMPENDIUM.md`](ZENODO_RELEASE_10_88_FIN_CURRENT_STATE_COMPENDIUM.md) / [`FIN_Release_10_88_Current_State_Compendium_EN.pdf`](FIN_Release_10_88_Current_State_Compendium_EN.pdf).
- Unified kernel genealogy, spectral core, internal observer relative scaling, operational testing criteria, and claim-graded status classification.

### Releases 10.89–10.97 — Blockwise Algebra, Causal Cones & Persistence (August 2026)
- **Release 10.89 (ST462–ST551):** Symmetric vacuum, blockwise algebra, and operational duality under \(D_{12}\) symmetry.
- **Release 10.90 (ST552–ST641):** Refinement speed no-go bounds, causal light-cone speed limits, and two-clock continuums.
- **Release 10.91 (ST642–ST731):** Conditional \(3+1\) Maxwell field emergence and strict declared-class no-go boundaries.
- **Release 10.92–10.93 (ST732–ST911):** Current-strict total-no-go theorems across declared ansatz classes; exhaustion of heuristic routes.
- **Release 10.94–10.95 (ST912–ST1091):** Ontological persistence of the nadsoliton under collision and heat bath; distinction between topological protection and pattern/coherence annihilation.
- **Release 10.96–10.97 (ST1092–ST1271):** Operator-to-total-state nonuniqueness theorem; strategic three-package architecture separating mathematical kernels, operational testers, and laboratory records.

### Releases 10.98–11.03 — Operational Process-Testing & Action Falsification (August–September 2026)
- **Release 10.98–10.99 (ST1272–ST1451):** Frozen operational apparatus (OA) dual-dynamics discrimination; interval-certified sequential inference and symmetric grid decisions.
- **Release 11.00–11.01 (ST1452–ST1631):** Structured uncertainty calibration; noncirculant continuums and fiber identifiability.
- **Release 11.02 (ST1632–ST1721):** Lindblad generator KKT systems and quantum clock no-go theorems.
- **Release 11.03 (ST1722–ST1811):** Repository-wide falsification of reverse closure from strict kernel to \(L_{\text{total}}\).

### Releases 11.04–11.09 — Cell Complexes, Hodge Energy & Fixed Points (September 2026)
- **Release 11.04–11.05 (ST1812–ST1991):** Direct gauge-Dirichlet obstructions; proof that automorphism-natural cell complexes cannot close required physical flux without external topological charges.
- **Release 11.06–11.07 (ST1992–ST2171):** Hodge energy classification on 1-forms and 2-cells; equivariant \(d_1\) boundary operators and 3-point source limits.
- **Release 11.08–11.09 (ST2172–ST2351):** Loop-product hidden synergy in pair kernels; ternary instruments and fractal fixed-point mappings.

---

## 6. Recent Research Suite: Programs ST2352 through ST3251

### 6.1. Triangular Prism Geometry & Ternary Gibbs Dynamics (ST2352–ST2621)
- Formulated the discrete triangular prism cell complex (\(C_3 \times K_2\)) as the minimal non-trivial vertical refinement of the discrete cycle.
- Proved that ternary Gibbs configurations generate stable internal sector orbits, but dynamic RG flow does not eliminate fiber gauge freedom.

### 6.2. Kinetic Rates & Microscopic Time-Arrow Obstruction (ST2622–ST2711)
- Proved that every positive local detailed-balance rate depending on energy difference \(\Delta\) is of the form:
  $$r(\Delta) = \exp(-\Delta / 2) h(\Delta)$$
  where \(h\) is an arbitrary positive even function. Analyticity, locality, and \(D_{12}\) covariance do not select \(h\).
- Proved that on the 12-cube cycle space (dimension 20,481), opposite circular currents produce identical Gibbs stationary states and identical positive entropy production rates. Entropy production magnitude **does not select a time orientation**.

### 6.3. Clock-Source Torsor & Operational Completion (ST2712–ST2891)
- Rescaling the continuous generator \(Q \to cQ\) preserves the stationary distribution, detailed balance, and untimed jump sequences. Generator scale \(c\) is an \(\mathbb{R}_+\) clock torsor invisible without externally calibrated timestamps.
- **Five-Axis Operational Completion (ST2802–ST2891):** Proved that a minimal operational characterization requires five typed axes:
  $$\vec{\theta} = (\theta_{\text{eq}}, q_{\text{refine}}, a_{\text{kinetic}}, c_{\text{clock}}, j_{\text{current}})$$
  The five-record observable map has full local Jacobian rank 5 and determinant \(-3.3901 \times 10^{-4}\). No smooth one-scalar parameterization can cover an open neighborhood of this completion family.

### 6.4. Exact Prism RG Fixed Points (ST2892–ST2981)
- The exact non-perturbative prism renormalization group recursion is:
  $$m' = m \frac{1 + r}{1 + r m^2}, \qquad \text{where } m = \tanh\theta, \quad r = \tanh(q)^3$$
- **Fixed Point Theorem:** For any vertical refinement \(q \neq 0\), there is **no finite non-zero fixed point**.
  - Positive \(q > 0\) strictly amplifies \(|\theta|\) toward complete order (\(\pm 1\)).
  - Negative \(q < 0\) attenuates \(\theta\) toward zero.
  - The map never selects the polarity/sign of \(\theta\).

### 6.5. Summable Refinement Hierarchies (ST2982–ST3071)
- For an infinite sequence of refinement layers \(q_n \ge 0\), the multi-level hierarchy has a finite non-zero limit \(\theta_\infty\) if and only if:
  $$\sum_{n=1}^\infty \tanh(q_n)^3 < \infty$$
- For power-law decay \(q_n = C n^{-\alpha}\), the sharp convergence threshold is **\(\alpha > 1/3\)**.
- Every exponential or dyadic decay \(q_n \sim 2^{-\gamma n}\) (\(\gamma > 0\)) is cubically summable and produces a stable, bounded hierarchy.

### 6.6. Dyadic Kernel Recurrence (ST3072–ST3161)
- Tested the dyadic kernel sampling recurrence \(q_n = K_{\text{strict}}(2^n)\).
- Rigorous envelope bound:
  $$|\tanh K_{\text{strict}}(2^n)|^3 \le 2^{-3\eta n} = 2^{-5.4 n}$$
- Absolute cubic sum converges rapidly to \(0.091068\), yielding \(\theta_\infty \approx 0.326689\) for \(\theta_0 = 0.3\).
- **Domain Boundary:** For \(n \ge 3\), the recurrence evaluates distances \(d \in \{8, 16, 32, \dots\}\) outside the frozen \(Z_{12}\) cycle domain \(\{1, \dots, 6\}\), functioning as an extrapolated metric ansatz rather than a derived operator property.

### 6.7. Layer Distance Embedding vs. Functorial Refinement (ST3162–ST3251)
- Audited the exact product refinement Laplacian:
  $$A_{24} = A_{12} \otimes I_2 + I_{12} \otimes B_q, \qquad B_q = \begin{pmatrix} q & -q \\ -q & q \end{pmatrix}$$
- **Intertwining Theorem:** The coarse isometry \(R = I_{12} \otimes \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix}\) satisfies:
  $$A_{24} R = R A_{12}, \qquad \exp(-A_{24} t) R = R \exp(-A_{12} t)$$
- **Weight Transport vs. Distance Doubling:**  
  Exact product refinement transports horizontal edge conductances functorially (\(W_{\text{fine}}(2d) = W_{\text{coarse}}(d)\)), **not** by re-evaluating the kernel at doubled label distances (\(K(2d)\)).
- **Even-Subspace Parity Leakage No-Go:**  
  A naive circulant \(Z_{24}\) Laplacian with weights \(K(d)\) leaks probability into odd sites: the parity-even injection \(E e_i = e_{2i}\) is **not** an intertwining operator for any connected circulant. The true coarse invariant subspace is **pair-uniform**, not parity-even.

---

## 7. Rigorous Epistemic Boundaries

To maintain scientific integrity, the project explicitly demarcates proven results from open challenges:

| Domain | Proven & Certified in FIN | Open / Obstructed (No-Go) |
| :--- | :--- | :--- |
| **Spectral Algebra** | 12-cycle self-adjoint generator \(A\), spectral measure, discrete Green function, exact resolvent bounds. | Continuum differential forms, infinite-dimensional Sobolev spaces. |
| **Symmetry & Selection** | \(D_{12}\) dihedral group symmetries, rational Krawczyk certification, nested dual SDPs. | Unique vacuum selector without explicit symmetry-breaking premise (`QW-2191` active). |
| **Dynamics & Time** | Reversible continuous-time Markov generators, unitary phase shifts, explicit jump chains. | Microscopic arrow of time (opposite currents yield equal entropy production); physical SI second. |
| **Multilevel RG** | Exact prism RG map \(m' = m(1+r)/(1+rm^2)\), cubic summability criterion \(\alpha > 1/3\). | Derivation of layer spacing or physical length units; continuous conformal field theory. |
| **Standard Model & GR** | Finite mathematical shadows, 12-state structural representations. | Tensor-resolved Einstein equations, Yang-Mills gauge actions, Dirac fermions, physical Higgs mechanism. |
| **Empirical Status** | Fail-closed validation schemas, interval decision bounds. | Blinded held-out laboratory verification; independent multiteam replication. |

---

## 8. Key Repository Documents & Entry Points

- **Primary Monograph:** [`ZENODO_RELEASE_10_88_FIN_CURRENT_STATE_COMPENDIUM.md`](ZENODO_RELEASE_10_88_FIN_CURRENT_STATE_COMPENDIUM.md) / [`FIN_Release_10_88_Current_State_Compendium_EN.pdf`](FIN_Release_10_88_Current_State_Compendium_EN.pdf)
- **Repository Guardrails Specification:** [`AGENTS.md`](AGENTS.md)
- **Persistence & Architecture Releases:**
  - [`RELEASE_10_95_TOTAL_NADSOLITON_ONTOLOGICAL_PERSISTENCE_EN.md`](RELEASE_10_95_TOTAL_NADSOLITON_ONTOLOGICAL_PERSISTENCE_EN.md)
  - [`RELEASE_10_97_STRATEGIC_THREE_PACKAGE_ARCHITECTURE_EN.md`](RELEASE_10_97_STRATEGIC_THREE_PACKAGE_ARCHITECTURE_EN.md)
  - [`RELEASE_10_98_FROZEN_OA_DUAL_DYNAMICS_DISCRIMINATION_EN.md`](RELEASE_10_98_FROZEN_OA_DUAL_DYNAMICS_DISCRIMINATION_EN.md)
  - [`RELEASE_10_99_INTERVAL_CERTIFIED_OA_AND_SEQUENTIAL_INFERENCE_EN.md`](RELEASE_10_99_INTERVAL_CERTIFIED_OA_AND_SEQUENTIAL_INFERENCE_EN.md)
  - [`RELEASE_11_00_STRUCTURED_UNCERTAINTY_CALIBRATION_REFINEMENT_EN.md`](RELEASE_11_00_STRUCTURED_UNCERTAINTY_CALIBRATION_REFINEMENT_EN.md)
  - [`RELEASE_11_03_REPOSITORY_WIDE_FALSIFICATION_AND_SOURCE_ACTION_GATE_EN.md`](RELEASE_11_03_REPOSITORY_WIDE_FALSIFICATION_AND_SOURCE_ACTION_GATE_EN.md)
  - [`RELEASE_11_09_FRACTAL_FIXED_POINT_AND_TERNARY_INSTRUMENT_EN.md`](RELEASE_11_09_FRACTAL_FIXED_POINT_AND_TERNARY_INSTRUMENT_EN.md)
- **Execution & Reproduction Suites:**
  - `fin_st3072_st3161_dyadic_kernel_recurrence.py`
  - `fin_st3162_st3251_layer_distance_embedding.py`
  - `fin_st2892_st2981_prism_rg_fixed_points.py`
  - `fin_st2982_st3071_summable_refinement_hierarchy.py`

---

## 9. Podsumowanie w Języku Polskim (Stan Badań od Release 10.88 Wzwyż)

### 9.1. Istota i Status Teorii
Fractal Information Nadsoliton (FIN) to bezwymiarowa, dyskretna teoria spektralno-informacyjna badająca pierwotną samoorganizację informacji w stabilnym stanie solitonicznym (nadsoliton).
- **Domknięcie Matematyczne:** Ściśle udowodnione na 12-cyklu \(\mathbb{Z}/12\mathbb{Z}\) przy użyciu generatora Dirichleta \(A = sI - W\), algebry blokowej, inkluzji Krawczyka oraz ścisłych równań grupy renormalizacji.
- **Status Teorii Wszystkiego (ToE):** **Teoria nie stanowi jeszcze fizycznej Teorii Wszystkiego**. Nie wyprowadza ona czasoprzestrzeni ciągłej, lagranżjanu Modelu Standardowego, tensora Einsteina, fizycznego upływu czasu ani jednostek układu SI.
- **Status Empiryczny:** Brak niezależnej, zewnętrznej weryfikacji laboratoryjnej.

### 9.2. Wiążące Zasady Projektowe (Guardrails)
1. **Podział Jądra:** \(K_{\text{legacy\_ont}}\) jest jedynie pomostem pośrednim. Domkniętym profilem roboczym jest \(K_{\text{strict\_gate}}\). Zakaz cichego podstawiania jąder oraz zakaz transferu ról fizycznych bez formalnego dowodu.
2. **Ontologia Nadsolitona:** Nadsoliton jest pierwotną informacją. Pod nim nie ma żadnej innej warstwy. Obowiązuje ścisła relacja: `nadsoliton -> światło -> materia -> wyłaniający się obserwator`.
3. **Przeszkoda Selektora (`QW-2191`):** Ścisły rdzeń nie posiada unikalnego mechanizmu wyboru próżni bez zewnętrznego założenia o łamaniu symetrii.

### 9.3. Najważniejsze Wyniki Cyklu ST2352–ST3251
- **Brak Mikroskopowej Strzałki Czasu (ST2622–ST2711):** Przeciwbieżne prądy kołowe generują identyczny stan Gibbsa i identyczną dodatnią produkcję entropii.
- **Ranga Operacyjna (ST2802–ST2891):** Pełny opis operacyjny wymaga 5 niezależnych osi (równowaga, refinacja, kinetyka, zegar, prąd).
- **Ścisła Renormalizacja Pryzmatu (ST2892–ST2981):** Równanie \(m' = m \frac{1+r}{1+rm^2}\) dowodzi, że dla \(q \neq 0\) nie ma nietrywialnego punktu stałego; \(q > 0\) wzmacnia uporządkowanie, a \(q < 0\) tłumi je do zera.
- **Sumowalne Hierarchie (ST2982–ST3071):** Krytyczny próg zbieżności wynosi \(\alpha > 1/3\) dla \(q_n = C n^{-\alpha}\).
- **Transfer Wag w Refinacji Produktowej (ST3162–ST3251):** Refinacja \(A_{24} = A_{12} \otimes I_2 + I_{12} \otimes B_q\) splata dynamikę gruboziarnistą w podprzestrzeni parzysto-jednorodnej poprzez funktorowy transfer przewodności, a nie poprzez naiwne ponowne próbkowanie jądra na podwojonych odległościach.
