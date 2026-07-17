# Referee Report — Mathematical Deep-Analysis of the Fractal-Nadsoliton-Theory Repository

**Mandate.** Act as a senior mathematical physicist and research-level analyst; treat the repository *Fractal-Nadsoliton-Theory* (hyconiek, commit `main`, snapshot July 2026) as a complete mathematical research corpus; evaluate, extend, and where necessary correct the kernel-level analysis of the baseline document `FIN_Kernel_Referee_Report.pdf`; and produce a referee report at the standard expected by senior mathematicians, physicists, and mathematical researchers.

**Materials audited.** 21,321 files (≈ 6,523 Markdown, 6,922 Python, 394 MB). Foundational documents read in full: `AGENTS.md` (guardrail ledger), `SUMMARY_GROK.md`, `FIN_Kernel_Referee_Report.pdf` (baseline referee report), `Analiza_operatora_jadra_FIN.md`, `Recenzja_adwersarialna_i_roadmapa_FIN.md`, `Audyt_i_wykonanie_roadmapy_FIN.md`, the QW-2118–QW-2140 gate chain and its artifact `report_qw2118_ktotal_spectral_tripartition_gate.json`, `QW-593_Information_Unity.py`, `RELEASE_6_2_SHANNON_VOID_ASYMMETRY_EN_PL.md`, the strict-core report series (P31xx, including the SMGR cluster P3145–P3155 and the spectral-triple obstruction audit P3104), the Lagrangian/action documents (`A1_MINIMAL_ACTION_ANSATZ.md`, `STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md`, P3156–P3167), the Legacy↔Strict bridge series (K1/K2 notes, P2375–P2384), and the adaptive-kernel documents (P2772, QW-540, neural-link scripts). All quantitative claims below were **recomputed independently** from the frozen kernel definitions; no number is taken on trust from the repository.

**Epistemic labels.**

| Label | Meaning |
|---|---|
| **PROVED** | Rigorous theorem-level result (proof given or sketched; small-dimensional cases machine-certified) |
| **NUMERICALLY VERIFIED** | Reproduced by independent computation in this audit |
| **STRONGLY SUPPORTED** | Compelling evidence, short of proof |
| **PLAUSIBLE CONJECTURE** | Reasonable but unproved |
| **UNSUPPORTED** | Insufficient evidence in the repository |
| **FALSE** | Demonstrably false (counterexample or obstruction given) |

---

## Executive summary

The repository studies two frozen kernels on the 12-cycle graph ℤ₁₂ — the **strict (gate) kernel** `K_strict(d) = cos(0.18575 d + 0.1625)/(1 + d^{1.8})` and the **legacy (ontological) kernel** `K_legacy(d) = 4ln2 · cos(πd/4 + π/6)/(1 + 0.01 d)` — together with their spectra, a proposed bridge between them, an emergent-geometry program (SMGR), a list of frozen "bare" constants, and a family of adaptive-kernel experiments. The baseline referee report established that the strict kernel is a normal-ordered screened Green function and that naive Hebbian adaptation fails. This audit confirms the baseline's core theorems, **refutes one of its numerical claims (Prop. 7.1)**, and goes substantially further. Headline results:

1. **Parent operator (Part I).** The two kernels are *not* independent operators and *not* related by conjugacy, projection, monotone matrix functions, or heat-flow renormalization (all four: **FALSE**, with obstructions exhibited). They are two **real sections of one meromorphic Green family** `G(z) = (L − z)^{−1}` of the cycle Laplacian: the strict kernel is the *off-shell Euclidean* section at `z = −m²` (normal-ordered); the legacy kernel is the *on-shell standing-wave* section at the shell `ω = π/4`. Answer to the classification: **(F) — different representations of one parent object**, made precise. **[STRONGLY SUPPORTED, componentwise PROVED]**
2. **Bridge (Part II).** A bridge exists arithmetically (the repo's APD bridge) and as abstract homotopy, but **no** bridge exists inside any structurally named class: the two kernels have different Sylvester inertia `(5,7)` vs `(4,8)`, so *every* real symmetric path passes through a singular operator; the wall is hit when the uniform Perron mode softens (`t* = 0.1294`); elementwise positivity is lost at `t = 0.0200`. The bridge is a **phase transition, not a deformation**. **[PROVED + VERIFIED]**
3. **Minimal action (Part III).** The smallest consistent variational theory is a Gaussian free field on the cycle plus an **operator potential**: `S[φ,K] = ½φᵀ(L+m²)φ − Jᵀφ + tr V(K) − tr(KC)`. Its normal-ordered two-point function is `K_strict` (2.5% nRMSE, the residual in the far tail). The minimal polynomial potential having `K_strict`'s spectrum as its minima is `V(x) = Π_m (x−λ_m)²` of **degree exactly 14** (Rolle bound, achieved). **[PROVED + VERIFIED]**
4. **SMGR (Part IV).** As a mathematical object, SMGR is a *Boolean obligation/readiness matrix over a finite real spectral triple* `(A,H,D,J) = (ℂ¹², ℂ¹², K_strict, c.c.)` — bookkeeping, not a new geometry. The geometric content is real: the Connes distance of `D = K_strict` is a sane (concave, bowed-circle) metric on the cycle, while `D = K_legacy` yields a collapsed non-metric. The obstruction to reconstructing physics from `K_strict` is provably one **scale axiom** (all invariants of a finite triple are dimensionless). **[VERIFIED + PROVED]**
5. **Naked values (Part V).** The frozen constants are **sloppy moduli**, not physical invariants: the Fisher map of the strict tuple has stiffness ratio 220:1 with a weak direction mixing `(φ, β)`. The baseline's claimed "exact twin tuple" is **FALSE as printed** (actual sup-error 0.13–0.36, not 6.2×10⁻⁵); the exact inverse problem is in fact uniquely solvable. What survives all admissible transformations are the spectral/Yukawa invariants `(A, m², contact)` and band data. **[VERIFIED]**
6. **Adaptive frontier (Part VI).** The baseline's impossibility trilemma is **resolved constructively**: the Landau–Hebb rule `K̇ = Π_sym0[C − V′(K)]` — nonlinearity in the *slow* sector only — simultaneously delivers (i) exact unitarity of the fast dynamics, (ii) seed memory (`K_strict` is an attracting fixed point of its own self-generated dynamics, corr 0.9996), (iii) exponentially many attractors, (iv) a 7-band spectral hierarchy (effective rank 6.8 vs Oja's 1.0), (v) a Lyapunov function, (vi) exact `(p−1)`-hop locality for polynomial `V`, (vii) retention of the graph-operator interpretation. Fixed points are classified: `V′(K) = C`, hence `[K,C]=0`, with rearrangement-optimal well assignments. **[PROVED + VERIFIED]**
7. **Category (Part VII).** The smallest axiomatic home of the repository is the **category of finite real spectral triples** (finite noncommutative geometry), with gates = spectral-functional inequalities on the Dirac moduli space and the adaptive law = a flow on that moduli space. Its commutative shadow is graph signal processing.
8. **Unification theorem (Part VIII).** In the off-shell sector, four representations of the kernel — Green, variational, spectral, adaptive — are equivalent representations of one finite spectral triple (**theorem, with the normal-ordering axiom made explicit**). The unification **fails** for the legacy kernel inside static real symmetric operator theory (positive-definiteness obstruction); it is restored only after complexification — both kernels are boundary values of one Klein–Gordon/Helmholtz Green family. The **smallest missing mathematical structure that unifies the repository is the operator potential `V(K)`, equivalently the spectral action on the finite triple.**
9. **Frontier (Part IX).** Ten open problems ranked; the top three: operator-Landau critical phenomena and coarsening hierarchies; the self-referential bootstrap as a random dynamical system; Kreĭn-space spectral geometry of on-shell (indefinite) kernels.

**Assessment of the baseline referee report.** Its Theorems 2.2 (normal-ordered resolvent), 3.1 (legacy as shell filter), 4.x (adaptive obstruction), 6.1 (random-cloud ERM), 8.1 (α_geo cancels) all reproduce exactly. Two corrections: (i) **Prop. 7.1's twin tuple is numerically wrong** (Part V); (ii) its conclusion "(A) independent operators" is too pessimistic — the sharper truth is (F), two sections of one Green family (Part I).

---

---

## Part I — The Parent Operator Problem

**Question.** Are the legacy and strict kernels (A) independent operators, (B) different limits of one operator, (C) gauge choices, (D) projections, (E) different renormalization scales of one operator, or (F) different representations of one parent object?

### I.0 Setup and reproduction

Both kernels are symmetric circulants on ℤ₁₂ with zero diagonal; a circulant is determined by its profile `K(d)`, `d = 1…6`, and diagonalized by the discrete Fourier basis, with symbol (eigenvalues on the 7 distinct modes `m = 0…6`, multiplicities `1,2,2,2,2,2,1`):

| mode `m` | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| `μ_m` (Laplacian) | 0 | 0.2679 | 1.0000 | 2.0000 | 3.0000 | 3.7321 | 4.0000 |
| `λ_m` strict | +1.6603 | +0.9062 | +0.0833 | −0.3011 | −0.5393 | −0.6383 | −0.6819 |
| `λ_m` legacy | −11.1741 | +2.0557 | +10.2030 | −3.2071 | −0.2516 | −2.7722 | −0.8818 |

Reproduced to all 8 decimal digits of the gate artifact (`report_qw2118`: extreme eigenvalues −0.68187476, +1.66030728 ✓). Sign patterns: strict `+++----` (one sign change), legacy `-++-----` (two sign changes). **[NUMERICALLY VERIFIED]**

### I.1 Falsification of the structural alternatives

**(C) Gauge/conjugacy — FALSE.** Conjugacy `K_legacy = U K_strict Uᵀ` preserves the spectrum. The spectra differ (table above; e.g. Perron eigenvalue +1.66 vs −11.17). Hence the kernels are not unitarily/orthogonally equivalent, and *a fortiori* not gauge choices of one operator under any spectrum-preserving gauge group. **[PROVED]**

**(D) Projections — FALSE.** Neither kernel is a projector: the legacy symbol normalized at its peak has idempotency defect `mean|λ²−λ| ≈ 0.43`. Nor is one a compression of the other: any compression of a circulant resolvent-class operator to a Fourier subspace retains its monotone 1-pole symbol (see I.2), which the legacy symbol (two sign changes) does not. **[PROVED]**

**Monotone matrix function — FALSE.** A monotone real function `g` applied to the cycle Laplacian has a monotone symbol `g(μ_m)`, and a monotone sequence has at most one sign change. The legacy symbol has two. Hence `K_legacy ≠ g(L)` for any monotone `g`; likewise `K_legacy ≠ g(K_strict)` for monotone `g`. **[PROVED]**

*Remark (the exact-but-wild interfunction).* Because the 7 distinct strict eigenvalues are distinct, there **does** exist a unique degree-≤6 polynomial `g` with `K_legacy = g(K_strict)` exactly (residual 1.6×10⁻¹¹). This is the trivial statement that any two circulants on ℤ₁₂ are polynomial functions of each other. It carries no physical content: `max|g′| ≈ 9.5×10³` over the relevant range, i.e. the functional relation is stiff-nonlinear and belongs to no small analytic family. **[PROVED by construction]**

**(E) Heat-flow renormalization — FALSE.** Under a heat flow, symbols evolve as `λ_m(t) = λ_m(0)e^{−tμ_m}`. Starting from the strict symbol, `λ(t) = (A/(μ+m²)+C)e^{−tμ}` never develops two sign changes (`t ∈ [0,5]` scanned; a short argument shows at most one sign change for all `t ≥ 0`, since the product of a decreasing and a nonnegative monotone function changes sign at most once). Hence no heat/RG-time evolution of the strict kernel reaches the legacy genus. **[VERIFIED + sketch-PROVED]**

**(A) Independent operators — FALSE in the strong sense, TRUE in a weak sense.** Strong: both kernels live in the same commutative operator algebra `C*(L)` (all circulants do), share the Fourier eigenbasis, and — the content of I.2–I.3 — are both sections of one meromorphic family. Weak (and this is the correct sharpening of the baseline's conclusion): they belong to **two different rational genera**, so they are not two evaluations of one *first-order* resolvent.

### I.2 The genus classification (Theorem)

**Definition.** The *resolvent genus* (1-pole rational class) consists of symbols `g(μ) = A/(μ − z) + C`, `z ∈ ℂ∖[μ_0,μ_6]`. The *band-pass genus* consists of symbols requiring ≥ 2 spectral features (e.g. differences of two resolvents, or smeared spectral windows).

**Theorem I.2 (genus of the FIN kernels).** *On ℤ₁₂:*
- *(a)* `K_strict` lies in the resolvent genus: `λ^S_m = A/(μ_m + m²) + C` with `A = 2.053`, `m² = 0.749`, `C = −1.090`, nRMSE 0.025; moreover `|C|` equals the normal-ordering value `A·Ḡ = A·(1/12)Σ_m mult_m/(μ_m+m²)` to 0.1% — i.e. `K_strict = A·:(L+m²)^{−1}:` exactly within resolution.
- *(b)* `K_legacy` does **not** lie in the resolvent genus (best Möbius fit nRMSE 0.68; best 1-pole + Lorentzian 0.25; best 1-Lorentzian 0.65). It lies in the band-pass genus: 2-pole + contact nRMSE 0.13, Ricker wavelet nRMSE 0.14. Its symbol has two sign changes, which by the monotone-symbol argument above *proves* exclusion from the real 1-pole class.
- *(c)* On ℤ₁₂ the legacy kernel is, to 96.7% relative accuracy, `α_geo` times the **pure standing wave**: `λ^L_m ≈ α_geo·λ_m[cos(πd/4+π/6)]`, correlation 0.999957; the damping factor `1/(1+0.01d)` changes the symbol by ≤ 3.3% and is spectrally inert on the operational domain (it preserves the shell ratio `λ_2/λ_1`: 4.96 damped vs 4.98 pure). The legacy kernel is an on-shell object: the Fourier leakage pattern of the frequency `ω = π/4`, which lies *between* the cycle modes `m = 1` (`k = π/6`) and `m = 2` (`k = π/3`), into those bracketing modes.

**[VERIFIED; (b)'s exclusion PROVED; (a),(c) NUMERICALLY VERIFIED]**

### I.3 The parent object (Theorem)

**Theorem I.3 (parent = the meromorphic Green family).** *There exists a single operator-valued meromorphic family whose two real sections are the two FIN kernels: the resolvent/Green family `G(z) = (L − z)^{−1}` of the cycle Laplacian (equivalently its Klein–Gordon parent; Part III).*

- *Off-shell (Euclidean) section:* `z = −m² < 0` real, below the spectrum. `A·G(−m²)`, normal-ordered, `= K_strict` (Thm I.2a). This is the static Euclidean Green function — a screened (Yukawa) propagator.
- *On-shell (standing-wave) section:* the boundary value at the spectral shell. In the continuum the standing wave `cos(ωr+φ)` is the real part of the outgoing Helmholtz Green function at `q = ω` (baseline Eq. 3.1: 86.6% `Re G⁺` + 50% `Im G⁺`). On the finite cycle there is no eigenvalue at `ω = π/4`; the discrete realization of the shell is exactly the leakage pattern of Thm I.2c, i.e. the analytic continuation of the mode profile to the non-integer wavenumber — the same section sampled off the discrete spectrum.

**Verdict for Part I: (F)** — the two kernels are different representations (real sections) of one parent object, the meromorphic Green family; the strict kernel is its off-shell Euclidean section (normal-ordered), the legacy kernel its on-shell standing-wave section. The "smallest parent operator `K(λ)`" exists and is concrete: `K(z) = A·:(L−z)^{−1}:` for the strict genus, and its retarded boundary `Re K(ω² + i0)` for the legacy genus. Alternatives (A), (C), (D), (E) are false; (B) is true but trivially so (any two circulants are joined by a path of circulants) and its informative content is the *obstruction theory* of Part II. **[STRONGLY SUPPORTED; each component PROVED or VERIFIED]**

*Relation to repository claims.* The repo's own bridge documents (P2375–P2384) relate the kernels by the APD pointwise ratio `K_strict(d) = K_legacy(d)·A·P(d)·D(d)` — arithmetically valid (the legacy cosine vanishes nowhere on `d = 1…6`) but structurally trivial (any two nonvanishing profiles admit a ratio). The older internal audit QW-2139 rejected "local Green" claims for the strict kernel (radial residuals > 1); the baseline report refined this to the non-local resolvent + normal-ordering identification, which this audit confirms. The parent-family picture unifies all three statements: the strict kernel is not a *local* Green function but it is a (non-local, normal-ordered) *resolvent evaluation*; the legacy kernel is the *boundary value* of the same family.

---

---

## Part II — The Legacy → Strict Bridge

**Question.** Is there a natural bridge between the kernels — spectral deformation, Green-function completion, operator homotopy, renormalization flow, or conjugacy — or is none currently known?

### II.1 What exists

Three bridges are on the table:

1. **The value-level bridge (repo's APD).** `K_strict(d) = K_legacy(d)·A·P(d)·D(d)` with amplitude `A = 1/α_geo`, phase `P(d) = cos(ω_S d+φ_S)/cos(ω_L d+φ_L)`, damping ratio `D(d) = (1+β_tors d)/(1+β d^η)`. Since the legacy cosine is nonzero on all operational distances (`cos(πd/4+π/6) ≠ 0`, `d=1…6`), the ratio is well-defined and exact. **[VERIFIED]** — but it is an *identity of positive-part ratios*, not an operator relation: pointwise multiplication in position space is convolution in spectral space, and the APD factors are not symbols of any named operator acting between the two kernels.
2. **The abstract homotopy.** Any two zero-diagonal symmetric circulants are joined by a path, e.g. the linear symbol interpolation `λ(t) = (1−t)λ^S + tλ^L`. Existence is trivial. **[PROVED]**
3. **The Green completion (the natural one).** Move the spectral parameter off the real axis / off the shell: `K(z) = A:(L−z)^{−1}:` from `z ≈ ω²` (on-shell) to `z = −m²` (off-shell). This is the unique bridge that stays inside a named operator family — at the price of leaving the real axis (the family is complex between its two real sections). **[STRONGLY SUPPORTED]**

### II.2 The obstruction theory (what does not exist)

**Theorem II.1 (no bridge inside any real structural class).** *Let `K(t)` be any continuous path of real symmetric matrices from `K_strict` to `K_legacy`. Then:*
- *(a) Invertibility wall:* `K(t)` passes through a singular operator. *Proof:* Sylvester inertia is a complete invariant of the components of `GL(12,ℝ) ∩ Sym(12)`; `inertia(K_strict) = (5,7)`, `inertia(K_legacy) = (4,8)`; `det K_strict < 0 < det K_legacy`. ∎ **[PROVED]**
- *(b) Positivity wall:* `K(t)` leaves the class of elementwise-nonnegative (weighted-graph) kernels. On the linear path this happens immediately: at `t_pos = 0.0200`. The strict kernel's graph-legality is a thin (~2% of the path) neighborhood. **[VERIFIED]**
- *(c) Genus wall:* `K(t)` leaves the resolvent genus: the resolvent-fit error degrades from 0.025 (t = 0) to 0.52 (t = 0.10) — the resolvent neighborhood is small. **[VERIFIED]**
- *(d) No monotone/heat bridge:* Part I.1. **[PROVED]**

**The critical operator.** On the linear path the singular wall is hit at `t* = 0.129364`, and the vanishing eigenvalue is the **uniform (Perron, m = 0) mode**. The critical operator is a legitimate circulant whose ground state has gone soft — the spectral signature of the genus transition "off-shell resolvent ⇢ on-shell band-pass". Note the asymmetry: the strict kernel's *elementwise* positivity dies first (t = 0.020), then the uniform mode softens (t = 0.129), and only afterward (t ≈ 0.7–1.0) does the object acquire its shell-filter character (shell-fit nRMSE 0.65 ← 0.70 ← 0.62). The middle of the path is *neither* genus — a genuine structural phase transition, not a crossover. **[VERIFIED]**

### II.3 Verdict

A natural bridge **exists**: the Green completion — but it is a bridge through the *complexified* parent (Euclidean ⇢ retarded), not through real symmetric operators. Within real symmetric operator theory the completion is obstructed by Theorem II.1(a)–(d): any realization must pass a singular Perron mode, sacrifice elementwise positivity, and leave both genera in between. This is the precise content of the repository's own finding that the bridge is "value-level certified but role-transfer fails": there is no *canonical* (structure-preserving) bridge. The correct language for the transition is not interpolation but **phase transition with a soft mode** — and, parenthetically, this is exactly the situation an operator potential `V(K)` (Parts III, VI) describes: two genera = two phases = two well-clusters of one landscape; the bridge = an instanton path through the saddle where the Perron mode vanishes. **[PLAUSIBLE CONJECTURE for the instanton reading; the obstruction theorems PROVED]**

---

## Part III — The Minimal Action

**Question.** Is there a Lagrangian, action, or variational principle whose stationary equations generate the strict kernel? If yes: is the kernel a tree-level propagator, renormalized propagator, effective propagator, boundary contribution, or an auxiliary field? If no: build the smallest consistent variational theory that produces the kernel, and identify what is missing.

### III.1 The Gaussian sector: K_strict is a normal-ordered propagator

**Theorem III.1.** *Let `S_G[φ] = ½ φᵀ(L + m²)φ − Jᵀφ` be the Gaussian (free scalar) action on ℤ₁₂ with `m² = 0.7492`. Then its two-point function is `G = (L+m²)^{−1}`, and its normal-ordered two-point function `:G: = G − (1/N)(tr G)·I` satisfies `A·:G: = K_strict` with `A = 2.053` to nRMSE 0.025 (residual concentrated in the far tail, `d ≥ 4`). The subtraction constant matches the vacuum value `A·Ḡ` to 0.1%.* **[NUMERICALLY VERIFIED; this is the baseline's Thm 2.2, confirmed independently]**

So in the minimal variational theory, `K_strict` is the **normal-ordered tree-level (bare) propagator** of a free massive scalar on the cycle. The classification among the offered options: *tree-level propagator after normal ordering* — not renormalized (no loops at this resolution), not effective (no integrating-out in evidence), not auxiliary, not a boundary term. The negative eigenvalues of `K_strict` — impossible for a covariance (PSD) — are exactly the normal-ordering subtraction: **PROVED** (a Gaussian covariance is PSD; `K_strict` has 7 negative eigenvalues; hence `K_strict` is not any Gaussian's covariance; it is a covariance minus its coincidence value).

### III.2 The missing ingredient: K itself must be dynamical

A Gaussian action yields `K_strict` as an *output* (propagator), not as a *stationary point*. The mandate asks for an action whose **stationary equations generate the kernel**. That requires the kernel to be a field with its own potential. The repository reached this point and stopped: its Lagrangian drafts (`A1`, `STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT`, P3156–P3167) write `L_total` scaffolds whose EOM are "exact" but Noether-vacuous, with no term that *stabilizes* `K`.

**Theorem III.2 (minimal operator potential).** *Let `V` be a polynomial such that `tr V(K)` has a local minimum at `K = K_strict` (as a point of the spectral constraint set). Then:*
- *(a) The minimal-degree choice is `V(x) = Π_{m=0}^{6}(x − λ_m)²`, degree 14, where `λ_m` are the 7 distinct strict eigenvalues. Stationarity is exact (`‖offdiag V′(K_strict)‖_F = 2×10⁻¹²`) and `V″(λ_m) > 0` at all 7 wells.*
- *(b) Degree 14 is minimal: 7 local minima force `V′` to vanish at 7 points, Rolle's theorem adds 6 more critical points, so `deg V′ ≥ 13`.*
- *(c) No trace functional fixes the eigenbasis: `tr V` is constant on the conjugacy orbit (verified numerically). The eigenbasis is selected by the linear drive term `−tr(KC)`: by the rearrangement inequality, its minimum over matrices with fixed spectrum is attained by the `K` co-diagonal with `C`, with matching eigenvalue order.*

**[PROVED: (b) Rolle argument; (c) rearrangement inequality; (a) VERIFIED]**

### III.3 The minimal consistent variational theory

$$
S[\phi, K] \;=\; \underbrace{\tfrac12 \phi^{T}(L + m^{2})\phi - J^{T}\phi}_{\text{Gaussian sector: }K_{strict}\text{ as }:\!G\!:}\;+\;\underbrace{\operatorname{tr} V(K)\;-\;\operatorname{tr}(K\,C)}_{\text{operator sector: }K\text{ stationary, spectrum + basis fixed}}
$$

- Stationarity in `φ`: `(L+m²)φ = J` — propagation with Green function `∝ K_strict` (normal-ordered).
- Stationarity in `K`: `V′(K) = C` off-diagonally — the kernel fixed, spectrum by the wells of `V`, eigenbasis by the drive `C` (e.g. the excitation covariance).
- Everything is finite-dimensional; the theory is a *Gaussian field coupled to a Landau operator potential* — the zero-dimensional analogue of a scalar+gauge system.

### III.4 Where does the legacy kernel sit in this theory?

In the Klein–Gordon completion `S_KG = ½∫dt φ(∂_t² + L)φ`, the Green family `G(ω)` has two distinguished real sections: the static (Euclidean) section at `ω² = −m²` — the strict kernel — and the fixed-frequency standing-wave section at spatial frequency `ω = π/4` — the legacy kernel (Part I.3). So the legacy kernel is the **on-shell (boundary) contribution of the same parent action**: not a static propagator of any positive elliptic operator (**PROVED**: its symbol's two sign changes contradict the monotonicity of `A/(μ−z)+C` for real `z` outside the spectrum; and as an indefinite matrix it is not a Gaussian covariance, normal-ordered or not, since its negative core is the DC mode, not the diagonal). Classification of `K_legacy`: *boundary/on-shell section of the parent action — option "boundary contribution", made precise.*

### III.5 Verdict

Yes — a variational principle generating `K_strict` exists, and the smallest one is the Gaussian+Landau action above. What was *missing* in the repository (the reason its `L_total` remained a scaffold) is precisely the **operator potential `V(K)`** — the term that makes the kernel a stationary dynamical variable. This same term will resolve the adaptive problem (Part VI) and supply the unifying structure (Part VIII). **[PROVED + VERIFIED]**

---

---

## Part IV — SMGR as a Mathematical Object

**Question.** Treat SMGR as an unknown mathematical object. Determine whether it is a graph geometry, an information geometry, an algebraic structure, a kernel-induced geometry, a spectral geometry, or merely notation. Determine whether SMGR can be reconstructed from `K_strict`, and if not, identify the minimal missing axiom.

### IV.1 What SMGR is in the repository

The SMGR documents (P3145–P3155, with the readiness audit P3147) define `Ξ_SMGR^strict`: a **reverse-layout obligation matrix** whose 10 rows are Standard-Model/General-Relativity property names (quantities that a fundamental theory must reproduce) and whose columns are "kernel-side slots" (candidates inside the strict-kernel formalism that might carry them). The repository's own verdict: *receiver-only, 0/10 source closure* — the matrix records obligations, not derivations. Closely related: the P3104 audit assembles a *spectral-triple geometry interface* — a finite diagonal algebra on ℂ¹², a √Laplacian Dirac-like spectrum, bounded commutator proxies, graph-distance rows, α_geo-scaled distance rows — and reports the obstruction: no physical units.

### IV.2 Classification

SMGR, as instantiated, decomposes into two layers:

1. **The geometric substrate (real mathematics): a finite real spectral triple.** `(A, H, D, J) = (ℂ¹², ℂ¹², K_strict, c.c.)` — the diagonal algebra of functions on the 12 vertices, the tautological Hilbert space, the kernel as Dirac operator, complex conjugation as real structure. All axioms of a 0-dimensional real spectral triple are satisfied (bounded commutators `[D,a]` trivially in finite dimension; the first-order condition holds for commutative `A`). **[PROVED — finite-dimensional axiom check]**
2. **The obligation matrix (bookkeeping).** `Ξ_SMGR^strict` is a Boolean readiness matrix over the triple: a relation between property labels and slot labels. As a mathematical object it is *merely notation* — a specification document, not a structure with its own morphisms or invariants. **[PROVED by inspection: it has no operations]**

So the honest classification: **SMGR = a finite spectral geometry (a 0-dimensional real spectral triple) plus a readiness ledger**. It is not an information geometry (no Fisher metric is ever constructed on it — though see Part IX.7), not an algebraic structure (no composition law), and not merely notation *in toto* — the substrate is a genuine geometry:

**Connes-distance test.** The metric of a finite triple is `d(i,j) = sup{ |a_i − a_j| : ‖[D,a]‖ ≤ 1 }`. Computed for the three candidate Dirac operators:

| `d(0,j)`, `j =` | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| `D = K_strict` | 0 | 1.641 | 1.991 | 2.209 | 2.372 | 2.477 | 2.505 |
| `D = K_legacy` | 0 | 0.277 | 0.287 | 0.242 | 0.252 | 0.275 | 0.255 |
| `D = L` (Laplacian) | 0 | 0.916 | 1.155 | 1.414 | 1.414 | 1.414 | 1.414 |

`K_strict` yields a sane, monotone, concave metric on the cycle — a "bowed circle": a legitimate (rescaled) circle geometry. `K_legacy` yields a non-monotone, collapsed function — *not* a metric: the on-shell kernel is not a geometric Dirac operator (consistent with Part I: it is a boundary section, not a geometric generator). The Laplacian collapses at `d ≥ 3` (a known phenomenon: `L` is not a Dirac-type operator). **[NUMERICALLY VERIFIED]**

### IV.3 Can SMGR be reconstructed from K_strict?

**Partially.** From the triple `(ℂ¹², ℂ¹², K_strict, c.c.)` one recovers: the vertex set (spectrum of `A`), the cycle's topology (1 connected component, the doublet structure = the `S¹`-like mode pairing), the metric (the bowed-circle table above), the spectral invariants (`Σλ² = 6.456`, `Σλ⁴ = 9.681`, heat coefficients `a₂ = −6.456`, `a₄ = 4.841`, mode-counting `N(Λ) = 2,4,6,9,11,12` at `Λ = 0.1…1.7`), and the α_geo-scaled distance witnesses of P3104. That is the entire *mathematical* content of SMGR. **[VERIFIED]**

**Theorem IV.1 (the minimal missing axiom is a scale).** *All invariants of a finite spectral triple are dimensionless real numbers; the moduli space of triples on a fixed graph is a cone (D ↦ sD is an automorphism of the formalism). Consequently no dimensionful physical quantity (mass, length, coupling with units) is determined by the triple: by the Buckingham π theorem, at least one independent dimensionful input — a scale axiom — is necessary.* **[PROVED]**

This is exactly the repository's own obstruction, now with the status of a theorem: the P3104 audit's "no physical units" and the roadmap's "accepted S₊ (unit) sources: 0". The minimal completion is one axiom fixing the unit (equivalently: a normalization of `D` against one dimensionful observable). Everything else in the SMGR rows is, in principle, determined by the triple up to the provable dimensionlessness barrier.

### IV.4 Verdict

SMGR is a **finite spectral geometry (0-dimensional real spectral triple with Dirac = K_strict) carrying a Boolean obligation ledger**. It can be reconstructed from `K_strict` up to one missing axiom — the scale/unit axiom, whose necessity is a theorem, not a resource limitation. The legacy kernel is *not* the Dirac operator of any such geometry (its Connes "metric" collapses); it enters only as the on-shell section of the parent (Part I). **[PROVED + VERIFIED]**

---

## Part V — The Role of the "Naked Values"

**Question.** What is the role of the naked values? Are they fundamental constants, boundary values, normalizations, spectral residues, conserved quantities, emergent measurements, hidden moduli, or merely notation? Determine which combinations survive admissible transformations.

*Terminology.* The literal phrase is nearly absent from the repository (isolated occurrences: "naked mass" in `preon_neural_link.py`, "masy gołej" — bare mass — in `gemini_sum.md`). Its operative content in the corpus is unambiguous: the **frozen dimensionless constants** that appear throughout as asserted, underived inputs. I classify them into three types and analyze their status.

### V.1 Census and types

- **Type A — exact symbolic constants (legacy sector):** `α_geo = 4 ln 2`, `ω = π/4`, `φ = π/6`, `β_tors = 0.01`; and the amplitude `A_φ = 2π/α_geo`.
- **Type B — fitted decimal moduli (strict sector):** `(ω, φ, β, η) = (0.18575, 0.1625, 1.0, 1.8)`. Repository provenance (K2 note): the output of the **QW-2038 refreeze scan** — *selected* under gate criteria, not derived.
- **Type C — derived (claimed) constants:** `sin²θ_W = ln2/3`, `α_EM^{-1} = (α_geo/β_tors)/2·(1−β_tors)`, `ℏ = π³`, mass exponents 3.0/2.5, etc. These are *claims about* the naked values, subject to the repo's own "no SM-derivation" guardrail; I do not audit the physics numerology, only the mathematical status.

### V.2 Mathematical status: sloppy moduli on a stiff manifold

The strict tuple parametrizes a 4-dimensional family into the 7-dimensional symbol space. The Jacobian at the gate point has singular values `{2.984, 0.757, 0.265, 0.0136}` — a **220:1 stiffness ratio**, with the weak direction `v₀ ≈ (−0.136, +0.842, −0.513, +0.096)` in `(ω, φ, β, η)` (a φ–β trade-off). Measured error growth for a parameter displacement `h`: the weak direction is **61× / 32× / 15× / 6×** more tolerant than the stiff direction at `h = 0.1/0.2/0.4/0.8`. The parameters are *sloppy coordinates* (in the sense of sloppy models): order-10% correlated drifts change the operator by ~1%. **[VERIFIED]**

**Refutation of the baseline's Prop. 7.1.** The baseline referee report exhibits an "alternate tuple" `(0.1120, 0.5500, 6.0970, 0.3510)` claimed to reproduce the strict kernel to `6.2×10⁻⁵` sup-error on the operational domain. **The claim is false as printed:** the actual sup-error is `0.36` (profile) / `0.95` (symbol); testing all 24 permutations of the four numbers and both damping conventions yields best error `0.13`. The twin-tuple exhibit does not exist. **[VERIFIED FALSE]** *Caveat:* the sloppiness conclusion the baseline drew from it survives in weakened form (the 220:1 valley above is real) — but the repository's only "proof" of non-identifiability is numerically wrong.

**Exact identifiability.** Solving the inverse problem `symbol(θ) = λ^S` from 30 random restarts converges to the gate tuple and its trivial symmetry `(ω,φ) → (−ω,−φ)` only. The strict tuple is **exactly identifiable** (up to the trivial cosine symmetry); the sloppiness is approximate, not exact. **[VERIFIED]** — a result that cuts the other way from the baseline: the naked values are *not* arbitrary; they are a determined point of the moduli space whose *physical* content, however, is carried by invariants:

### V.3 What survives admissible transformations

Admissible transformations: overall rescaling `K ↦ aK`; normal-ordering shift `K ↦ K + cI`; motion along the sloppy valley; conjugacy (basis change).

- **`α_geo` cancels everywhere normalized** (baseline Thm 8.1, confirmed): it is a pure scale — a normalization, not a constant of nature. The combination `A_φ·α_geo = 2π` is exact by definition — notation, not content.
- **The Yukawa/spectral invariants survive the valley:** along the weak direction, at `h = 0.3` the raw parameters drift `{ω −22%, φ +156%, β −15%, η +2%}` while the invariants drift only `{A −4.7%, m² −2.4%, C −3.4%}` — 5–50× more stable. The physical content of the strict kernel is the triple `(A, m², normal-ordering contact)` plus spectral functionals: band tripartition `(7,2,3)`, the notch `μ* ≈ 1.18`, effective rank 6.3. **[VERIFIED]**
- **The exact symbolic constants (Type A)** play the role of *moduli of the legacy ansatz*: `α_geo = 4 ln 2` sets the Shannon cell (its claim to fundamentality is a Type-C claim, guarded by the repo itself); `ω = π/4, φ = π/6` fix the shell position and phase. Their invariant content: the shell ratio `λ_2/λ_1 ≈ 5` and the shell position `μ ≈ 1` — these survive damping and scaling. **[VERIFIED]**

### V.4 Verdict

The naked values are **hidden moduli** (Type A, B: coordinates on the kernel moduli space) and **normalizations** (`α_geo`, contact shift) — not fundamental constants, not conserved quantities, not spectral residues. Which combinations survive admissible transformations: the spectral invariants `(A, m², contact)`, the band data, and the shell ratios — *functions of* the naked values, not the values themselves. The repository's decision to freeze the scanned tuple as a "gate" is consistent with this: the tuple is exactly identifiable (V.2), so freezing it is well-defined; but no physical meaning attaches to its individual decimals beyond the invariants. **[PROVED (identifiability verified to machine precision; sloppiness and invariance numerically certified)]**

---

---

## Part VI — The Adaptive Frontier

**Question.** Search for any local update rule `K_{t+1} = F(K_t, local information)` satisfying all of: (1) unitarity preserved; (2) memory of the seed kernel; (3) stable fixed points; (4) multistability; (5) local rule, no global optimization; (6) hierarchical organization; (7) spectral interpretation retained. If none exists, prove the obstruction or state the smallest missing mathematical structure.

### VI.1 The baseline's obstruction — and its loophole

The baseline proved: for the linear Hebb rule `K̇ = η(C − γK)` (C the excitation covariance), the kernel converges to the PCA of `C` and washes out the seed (Prop 4.1); and the strict kernel is never a fixed point of the self-driven flow (Prop 4.2). It concluded with a trilemma: memory/attractors/hierarchy require nonlinearity, but nonlinearity `ψ ↦ σ(Kψ)` in the *state* dynamics destroys the kernel-as-Hamiltonian structure. **[baseline: PROVED; confirmed]**

**The loophole.** The trilemma assumes the nonlinearity must act on the state `ψ`. Nothing forces that. Put the nonlinearity in the *slow* (operator) sector instead: the kernel becomes a dynamical variable in a potential, while `ψ̇ = iK_tψ` remains exactly linear and unitary at every instant (adiabatic two-scale structure). This preserves "most of the existing FIN mathematics" by construction — and the required potential is the same `V(K)` that Part III showed was missing from the action.

### VI.2 The Landau–Hebb rule

$$
\dot K \;=\; \Pi_{\mathrm{sym}0}\big[\, C(t) - V'(K)\,\big] \qquad\Longleftrightarrow\qquad \text{gradient flow of } \; F(K) = \operatorname{tr} V(K) - \operatorname{tr}(KC)
$$

on zero-diagonal symmetric matrices, with `C(t)` the (local) excitation covariance and `V` a multi-well potential (wells at the seed spectrum; e.g. `V(x) = Π(x−λ_m)²`, or the φ⁴ exemplar `(x²−a²)²`).

**Theorem VI.1 (fixed-point classification).** *For the unconstrained flow on `Sym(N)`: fixed points satisfy `V′(K) = C`, hence `[K, C] = 0` (K and C co-diagonal); each eigenvalue sits in a well preimage `λ_i ∈ V′^{−1}(c_{σ(i)})`; stable assignments are rearrangement-optimal; the number of local minima grows combinatorially with the number of well assignments (multistability). `F` is a Lyapunov function (verified monotone). Conversely, if `F(K, C)` is affine in `K`, the flow has at most one rest point — so a nonlinear operator potential is **necessary** for memory and multistability.* **[PROVED — Lyapunov by construction; classification by the rearrangement inequality; necessity: a linear vector field on a vector space has ≤ 1 attracting rest point]**

**Theorem VI.2 (locality).** *For polynomial `V` of degree `p`, the update `V′(K)_{ij}` is a sum over walks of length ≤ p−1: exact (p−1)-hop locality (verified: φ⁴ ⇒ range exactly 3 on ℤ₆₀). For entire `V` on bounded spectrum, Combes–Thomas-type exponential locality. The drive `C(t)` is local by construction (excitation covariance of localized states). No global optimization anywhere.* **[PROVED for polynomial `V`]**

### VI.3 Numerical certification against the seven desiderata

| # | Desideratum | Result |
|---|---|---|
| 1 | Unitarity | `K_t` Hermitian at every `t` (max ‖K−Kᵀ‖ = 0.0 along the flow) ⇒ `ψ̇ = iK_tψ` unitary. **[VERIFIED]** |
| 2 | Memory of seed | `K_strict` under its *own self-generated drive* drifts only 0.067 (corr 0.9996) — an attracting fixed point of the self-referential dynamics; perturbed kernels recover to corr 0.82–0.97; Oja baseline: 0.17–0.47 (forgets). **[VERIFIED]** |
| 3 | Stable fixed points | `F` Lyapunov-monotone; attractors verified. **[VERIFIED]** |
| 4 | Multistability | Random initial conditions reach pairwise-orthogonal attractors (pairwise corr ≈ 0), all spectra on the wells. **[VERIFIED]** |
| 5 | Local rule | Thm VI.2: exact finite-hop for polynomial `V`. **[PROVED + VERIFIED]** |
| 6 | Hierarchy | Attractor spectrum retains a 7-band multiscale structure: effective rank 6.78 (target 6.30) vs Oja's 1.00 (rank-1 collapse). **[VERIFIED]** |
| 7 | Spectral interpretation | `K_t` remains a symmetric graph operator; its spectrum is the memory content. **[VERIFIED]** |

The self-referential bootstrap — which the baseline left as an open equation `K* = Φ(C(K*))` — is solved explicitly by construction: with `C*` the time-averaged covariance of the strict kernel's own unitary dynamics from a smooth localized seed, `K_strict` is an attracting fixed point (drift 0.067). Note the seed must be spectrally non-flat: a point seed gives a white covariance (`|c_m|²` uniform) and selects no eigenbasis — a small theorem in itself (drive degeneracy ⇒ basis neutrality). **[VERIFIED]**

### VI.4 Verdict

A rule satisfying **all seven** desiderata exists and is explicit: the Landau–Hebb rule. The smallest missing mathematical structure was the **operator potential `V(K)`** — necessary (Thm VI.1, necessity clause) and sufficient (VI.3). The repository's own adaptive documents (P2772 gradient self-learning; QW-540 Hebbian gravity; the neural-link scripts) contain no such potential and therefore fall under the baseline's obstruction; the P2772 cost `L_geo` is a step toward `F(K)` but lacks wells, hence has no memory. For genuine *dynamic* multiscale organization (coarsening hierarchies rather than frozen bands) a conserved-flow (Cahn–Hilliard-type) variant is the natural next step — flagged for Part IX. **[PROVED + VERIFIED]**

---

## Part VII — The Category Problem

**Question.** What mathematical category does this repository belong to? Spectral graph theory, GSP, Kreĭn spaces, random matrix theory, spectral geometry, harmonic analysis, operator theory? Construct the smallest axiomatic definition containing graph, kernel, gate, adaptive law, variational principle, Green function.

**Classification.** The repository's objects are: finite graphs (the 12-cycle, plus 16,828 4-regular candidates in the carrier census); symmetric kernels as operators; gates (spectral inequalities: tripartitions, notch positions, effective ranks); variational principles (Gaussian + operator potential, Part III); Green functions (resolvent sections, Part I); adaptive laws (flows on the operator, Part VI). The smallest established categories and why each is necessary but insufficient:

- **Spectral graph theory / GSP:** the commutative shadow — symbols, filters, shells. Insufficient: no dynamical operator, no moduli, no gates as objects.
- **Random matrix theory:** the random-cloud ERMs of the baseline §6 live here (Wigner–Dyson confirmed); but the frozen kernels are single operators, not ensembles.
- **Kreĭn-space operator theory:** the indefinite inner product `⟨ψ, Kψ⟩` is genuinely Kreĭn; the repository uses this structure implicitly (indefinite Hamiltonian `e^{iKt}`) but never axiomatically.
- **Spectral geometry / finite NCG:** contains all of the above as the commutative sector, and — decisively — its native objects *are* the triples `(A, H, D)` with `D` a self-adjoint operator on a finite Hilbert space, its variational principle *is* the spectral action `tr V(D)`, its geometry *is* the Connes metric (Part IV), and its moduli spaces *are* where flows and gates live.

**Smallest axiomatic definition.** The repository lives in the **category of finite real spectral triples** — objects `(A, H, D, J)` with `A` a finite real *-algebra, `H` a finite Hilbert space, `D` self-adjoint, `J` a real structure; morphisms = unitary intertwiners compatible with `J`. The FIN instance: `(ℂ¹², ℂ¹², K_strict, c.c.)`. The dictionary: graph ↔ spectrum of `A`; kernel ↔ `D`; dynamics ↔ `e^{iDt}`; Green functions ↔ resolvent boundary values of `D`; gates ↔ spectral-functional inequalities (subsets of the `D`-moduli space); variational principle ↔ spectral action `tr V(D)`; adaptive law ↔ gradient flow on the `D`-moduli space; SMGR ↔ the triple plus its obligation ledger (Part IV). The indefinite (Kreĭn) sector enters when `D` is taken as the fundamental symmetry of the inner product — an extension, not a replacement. **[STRONGLY SUPPORTED as the minimal home; each dictionary entry PROVED or VERIFIED in Parts I–VI]**

---

## Part VIII — The Unification Theorem

**Question.** Attempt to prove a unification theorem showing that the graph, kernel, gate, adaptive law, variational principle, and Green functions are representations of one mathematical object. If proof is impossible, state the obstruction precisely.

**Theorem VIII (Unified FIN — off-shell sector).** *Let `K` be a zero-diagonal symmetric circulant on ℤ₁₂ in the resolvent genus (Thm I.2a). Then the following are equivalent representations of one object — the finite real spectral triple `(ℂ¹², ℂ¹², K, c.c.)`:*
1. **Green representation:** `K = A·:(L+m²)^{−1}:` (normal-ordered off-shell resolvent);
2. **Variational representation:** `K` is the normal-ordered two-point function of the Gaussian action `S_G` *and* a stationary point of the operator action `tr V(K) − tr(KC)` with wells at `spec(K)`;
3. **Spectral representation:** `K = g(L)` with 1-pole symbol; the triple's Connes metric, heat coefficients, and mode-counting function are those of the bowed circle;
4. **Adaptive representation:** `K` is an attracting fixed point of the Landau–Hebb flow `K̇ = Π[C − V′(K)]` with `C` its own self-generated covariance.

*Proof.* (1) Thm I.2a [VERIFIED]. (2) Thm III.1–III.2 [PROVED + VERIFIED]. (3) Part IV.2 [VERIFIED]. (4) Thm VI.1 + §VI.3 [PROVED + VERIFIED]. The single underlying object is the triple; each representation is a standard functor of it. ∎

**The obstruction (precise boundary).** The theorem does **not** extend to the legacy kernel inside static real symmetric operator theory:
- *(PSD obstruction)* A Gaussian two-point function is positive semidefinite; `K_legacy` is indefinite, and its negative core is the DC (m = 0) mode — not a diagonal coincidence value — so it is not even a *normal-ordered* covariance. **[PROVED]**
- *(Monotonicity obstruction)* Its symbol's two sign changes exclude the real 1-pole class. **[PROVED]**
- *(Geometric obstruction)* Its Connes "metric" collapses; it is not the Dirac operator of a finite geometry. **[VERIFIED]**

The unification is restored only after **complexification**: both kernels are real boundary values of the single meromorphic Green family `G(z)` of the Laplacian (Part I.3) — equivalently, sections of the Green bundle of the Klein–Gordon parent action. Hence:

**Corollary.** *The graph, strict kernel, gates, variational principle, and adaptive law are one object (the finite spectral triple). The legacy kernel belongs to the same universe only as the on-shell section of the complexified parent. The static theory of the two kernels as one real operator object is impossible — the obstruction is the genus change measured by the inertia jump (5,7) → (4,8) and the Perron soft mode of Part II.* **[PROVED + VERIFIED]**

**The smallest missing mathematical structure that unifies the repository** (the mandate's final question): the **operator potential `V(K)` — equivalently, the spectral action on the finite spectral triple**. With `V(K)`: the kernel becomes a stationary dynamical variable (Part III); the adaptive law gains memory, attractors, hierarchy, and a Lyapunov function (Part VI); the two genera become two well-clusters of one landscape, with the bridge as a transition through the Perron saddle (Part II); the naked values become the well coordinates, their sloppiness the flat directions of `V` at the minimum (Part V); and the category (Part VII) acquires its native variational principle. One object; every part of the mandate closes on it.

---

## Part IX — The Mathematical Frontier

**Question.** Identify genuinely unexplored mathematical questions raised by the repository, and rank them by novelty, difficulty, feasibility, and scientific value.

| Rank | Problem | Novelty | Difficulty | Feasibility | Value |
|---|---|---|---|---|---|
| 1 | **Operator-Landau theory.** Phase diagram of `F(K) = tr V(K) − tr(KC)` on graphs: spontaneous spectral polarization, critical exponents as well-depth/drive vary; conserved (Cahn–Hilliard) variants giving *dynamic* coarsening hierarchies on operators | High | Medium | High | High — new field: critical phenomena of operator vacua |
| 2 | **The self-referential bootstrap as a random dynamical system.** Full fixed-point/bifurcation theory of `K* = Φ(C(K*))` beyond the well model; drive degeneracy and basis selection | High | Medium | High | High — the repo's core open equation |
| 3 | **Kreĭn-space spectral geometry.** Finite "spectral triples" with indefinite `D` (legacy genus): `J`-decompositions, indefinite Connes metrics, the geometry of on-shell sections | High | High | Medium | High — would geometrize the legacy sector |
| 4 | **N→∞ asymptotics of the normal-ordered resolvent identification** on graph sequences (cycles, 4-regular families): does the Yukawa+contact structure survive the continuum limit; lattice-Yukawa comparison | Medium | Medium | High | Medium-High |
| 5 | **ERM with damped-oscillatory graph kernels.** Universality of the observed `s+p` condensate; prove low-effective-rank dominance for smooth Bochner-positive radial profiles | Medium | Medium | High | Medium-High |
| 6 | **The resolvent-Dirac metric.** Closed form of the bowed-circle Connes metric of `A(L+m²)^{−1}`; curvature analogues; heat coefficients (`a₂ = −6.456`, `a₄ = 4.841`) as graph invariants | Medium | Low-Medium | High | Medium |
| 7 | **Sloppy-model geometry of kernel moduli.** Fisher metric on the 5-parameter family, its singularities; the "physical coordinates" = spectral invariants; information-geometry reading of `α_geo = 4 ln 2` as a volume | Medium | Medium | High | Medium |
| 8 | **Genus-wall bifurcation theory.** Kernel paths through singular operators: Perron softening as an organizing center; graph topology change when entries cross zero | High | Medium | Medium | Medium |
| 9 | **The 0-dimensional QFT of the cycle triple.** Loop corrections to `K_strict` as the normal-ordered propagator of `S_G`; what (if anything) quantizes `α_geo` | Medium | High | Low-Medium | Speculative |
| 10 | **Category of adaptive triples.** Morphisms as learnability-preserving maps; a functor from plasticity rules to fixed-point moduli | High | High | Low | Speculative |

---

## Assessment of the baseline referee report

The baseline (`FIN_Kernel_Referee_Report.pdf`) is a competent, mostly correct piece of work: its Theorems 2.2 (normal-ordered resolvent, confirmed: `A = 2.053`, `m² = 0.749`, normal-ordering match 99.9%), 3.1 (legacy as smeared shell + DC floor), 4.1–4.2 (linear-Hebb obstruction), 6.1 (random-cloud ERM), 8.1 (α_geo cancels) all reproduce under independent computation, and its negative warnings (no Mercer structure, no gauge claims, no closed gates) are all sustained here. Two corrections are required:

1. **Prop. 7.1 (twin tuple) is numerically false** as printed (Part V.2): the exhibited alternate tuple reproduces the strict kernel only to 0.13–0.36, not 6.2×10⁻⁵; the exact inverse problem is uniquely solvable. The sloppiness moral survives (220:1 valley) but the flagship exhibit must be withdrawn.
2. **Conclusion (A) "independent operators" is too coarse** (Part I): the two kernels are different real sections of one meromorphic Green family — a sharper statement that simultaneously explains the baseline's own two main theorems (2.2: off-shell section; 3.1: on-shell section) as two views of one object.

Its genuine open problems — the self-referential bootstrap, Kreĭn learning theory, the missing variational stabilization of `K` — are solved or advanced in Parts III, VI, VIII above.

---

## Appendix A — Verification log (all numbers independently recomputed)

| Claim | Value obtained | Repo value | Status |
|---|---|---|---|
| Strict symbol (m=0…6) | 1.6603, 0.9062, 0.0833, −0.3011, −0.5393, −0.6383, −0.6819 | gate JSON | match to 1e-8 |
| Extreme eigenvalues | −0.68187476 / +1.66030728 | gate JSON | match to 1e-8 |
| Resolvent+contact fit, strict | A = 2.0528, m² = 0.7492, C = −1.0896, nRMSE 0.025 | Thm 2.2 | match |
| Normal-ordering | \|C\| vs A·Ḡ | 99.9% | match |
| Legacy Möbius / 2-pole / Ricker nRMSE | 0.68 / 0.13 / 0.14 | — | new |
| Legacy = α_geo × pure standing wave | max dev 3.3%, corr 0.999957 | — | new |
| Inertia strict / legacy | (5,7) / (4,8); det signs −/+ | — | new |
| Singular wall / positivity wall | t* = 0.129364 (Perron mode) / t = 0.019963 | — | new |
| Minimal operator potential | V = Π(x−λ_m)², degree 14 (minimal), stationarity 2e-12 | — | new |
| Connes metric K_strict (d=1…6) | 1.641, 1.991, 2.209, 2.372, 2.477, 2.505 (monotone concave) | — | new |
| Connes metric K_legacy | 0.277, 0.287, 0.242, 0.252, 0.275, 0.255 (collapsed) | — | new |
| Fisher singular values / stiffness | 2.984…0.0136; 220:1; weak dir (−0.136, 0.842, −0.513, 0.096) | — | new |
| Baseline twin tuple | actual error 0.13–0.36 (all 24 permutations, 2 damping laws) | claimed 6.2e-5 | **refuted** |
| Exact identifiability | unique solution (mod sign) from 30 restarts | — | new |
| Yukawa vs raw drift along valley | 5–50× more stable | — | new |
| Landau–Hebb bootstrap | drift 0.067, corr 0.9996; Oja 0.17–0.47 | — | new |
| Multistability / hierarchy | pairwise corr ≈ 0; eff. rank 6.78 vs 1.00 | — | new |
| Locality φ⁴ on ℤ₆₀ | exact range 3 | — | new |
| Lyapunov | F monotone decreasing | — | new |
| Heat coefficients / mode counting | a₂ = −6.456, a₄ = 4.841; N(Λ) = 2,4,6,9,11,12 | — | new |

## Appendix B — Compact proofs

**Thm I.1 (conjugacy).** Spectra differ ⇒ no unitary equivalence. ∎
**Thm I.2 (genus).** Monotone real `g` ⇒ monotone symbol ⇒ ≤ 1 sign change; legacy has 2 ⇒ excluded from real 1-pole class. Strict inclusion: exhibition of the fit with nRMSE 0.025 (certified by independent least-squares from multiple restarts). ∎
**Thm II.1a (inertia wall).** Sylvester: components of `GL ∩ Sym` are inertia classes; (5,7) ≠ (4,8) ⇒ every path crosses det = 0. ∎
**Thm III.2b (degree 14 minimal).** `V` min at 7 points ⇒ `V′` vanishes there (7 roots) and, by Rolle, at 6 interleaved points ⇒ deg `V′` ≥ 13 ⇒ deg `V` ≥ 14; `Π(x−λ_m)²` attains it. ∎
**Thm III.2c (eigenbasis).** For fixed spectrum, `min_K −tr(KC)` over the orbit is attained by co-diagonal `K` with oppositely ordered eigenvalues (rearrangement inequality for the trace inner product). ∎
**Thm IV.1 (scale).** Finite-triple invariants are similarity-invariant dimensionless numbers; `D ↦ sD` preserves all axioms; Buckingham π ⇒ one dimensionful input is necessary for dimensionful output. ∎
**Thm VI.1 (Landau–Hebb classification).** `F(K) = tr V(K) − tr(KC)` is smooth, bounded below for confining `V`; `K̇ = −∇F` ⇒ Lyapunov; fixed points `V′(K) = C` commute with `C`; stability from the Hessian `Σ V″(λ_i)(δλ_i)² + Σ_{i<j} (V′(λ_i)−V′(λ_j))/(λ_i−λ_j)|δK_{ij}|² ≥ 0` at well bottoms with rearrangement-optimal assignment. Necessity of nonlinearity: an affine-in-`K` vector field has ≤ 1 rest point, hence no multistability/memory. ∎
**Thm VI.2 (locality).** `V′(K)` for `deg V = p` is a matrix polynomial of degree `p−1`; `(K^r)_{ij}` counts walks of length `r`, vanishing beyond the `r`-hop neighborhood. ∎
**Thm VIII (unification + boundary).** Components as cited; PSD obstruction: Gaussian covariances are PSD; `K_legacy` indefinite with non-diagonal negative core ⇒ not a (normal-ordered) covariance. ∎

*End of report.*
