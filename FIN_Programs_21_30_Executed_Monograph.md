# FIN Beyond the Spectral Core

## Executed Programs 21–30, operational falsification, and inverse reconstruction of an action from the strict kernel

**FIN Research Supplement — Release 10.2**  
**Author:** Krzysztof Żuchowski  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Publication type:** Preprint  
**Version:** 1.0.0  
**Date:** 19 July 2026  
**Language:** English  
**License:** CC BY 4.0

## Confidence convention

Every principal statement is marked as follows.

- **[Proven]** means an analytic proof is supplied, or an exact finite theorem is reduced to auditable algebra.
- **[Strong evidence]** means a reproducible numerical result is stable under the declared perturbations, but is not promoted to a general theorem.
- **[Moderate evidence]** means a result is model-dependent or has not survived a broad enough alternative class.
- **[Speculative]** means a mathematically admissible research direction lacks a source theorem or decisive computation.
- **[Refuted]** means the declared claim has a proof-level obstruction or an explicit counterexample.

All times, rates, masses, actions, and distances in the finite calculations are dimensionless. No SI unit, physical clock, experimental apparatus, or physical interpretation is inferred from the strict operator alone.

---

# Executive summary

This monograph executes the ten studies that were previously recommended as Programs 21–30. It also performs a separate inverse-action audit requested after the roadmap was written. The audit begins from the empirical observation

\[
K_{\rm strict}\approx a_Y\,\mathcal N\!\left[(L+m_Y^2I)^{-1}\right],
\qquad
\mathcal N(X)=X-\operatorname{diag}X,
\]

and asks whether a smallest Lagrangian can be reconstructed from that approximate propagator relation.

The work produces five constructive results.

1. A targeted intervention set distinguishes four declared temporal models on the twelve-mode system with nonzero finite-shot margins. **[Strong evidence]** It is model discrimination, not unrestricted process-tensor tomography.
2. Two different positivity-preserving refinements reproduce the twelve-mode strict weights and converge to inequivalent limits: one nonlocal fixed-circumference jump operator and one local diffusive limit obtained from an exponentially completed tail. **[Proven]** Their coexistence proves that the twelve-mode data do not select a unique continuum.
3. A doubled incidence Dirac operator gives an exact local \(z=1\) refinement of the cycle Laplacian, while a weighted incidence operator gives an exact square root of the strict generator but is complete-graph nonlocal. **[Proven]** No tested construction has strict matching, locality, and a FIN-selected Lorentzian limit simultaneously.
4. A free massive Euclidean extension of the strict generator is reflection positive, with exact factorization of the reflection form. **[Proven]** A positive higher-derivative quadratic action supplies a counterexample: Euclidean positivity alone does not imply Osterwalder–Schrader positivity.
5. The inverse Gaussian action exists conditionally. The exact variational object is most cleanly written as a covariance-level functional containing a log determinant. **[Proven]** The displayed Gaussian-plus-operator-potential action does not by itself couple its field and kernel sectors, and a post-hoc potential does not uniquely or noncircularly select the strict kernel. **[Refuted as an unconditional uniqueness claim]**

The work also produces five limiting results.

1. Full two-slot process-tensor tomography on a generic \(d=12\) system has an operation-space scale \(d^8=429,981,696\). A modest instrument set is possible only after symmetry, low-rank, finite-memory, or candidate-model restrictions. **[Proven]**
2. Population-only records do not determine a bath: two continuous quantum semigroups have the same complete FIN population process but one-time Choi ranks 12 and 144. **[Proven]**
3. The declared finite-time self-covariance one-form is strongly nonintegrable. Its first odd instability is complex and is not a pitchfork certificate. **[Strong evidence]**
4. The exact no-go library passes five executable certificates, but independent Lean kernel replay was not completed because no Lean toolchain was configured. **[Partially completed]**
5. The blinded synthetic challenge reached 80% accuracy. It separates the declared FIN profile from many positive-circulant alternatives but loses specificity for several alternatives at the 5–20% perturbation scales. It validates a protocol, not FIN as physics. **[Moderate evidence]**

The deepest conclusion is therefore precise:

> The strict spectral object supports several mathematically consistent operational, variational, stochastic, hyperbolic, and Euclidean completions. Programs 21–30 narrow those completions and construct decisive tests, but they do not derive the missing state space, clock, instrument, environment, selector, dimensional scale, or empirical domain from the operator itself.

---

# Part I — Frozen mathematical input and research discipline

## 1. The strict finite object

The executed calculations use the twelve-mode radial matrix

\[
(K_s)_{xy}=k(d_{12}(x,y)),
\]

with zero diagonal and radial row

\[
(0,
0.469985673,
0.192043552,
0.091428614,
0.047029169,
0.024131223,
0.011070817).
\]

The declared row sum is

\[
s=1.660307278766099,
\]

while the nine-decimal row gives \(1.660307279\), a difference of \(2.34\times10^{-10}\). Define

\[
A_M=sI-K_s.
\]

With positive off-diagonal weights,

\[
\langle f,A_Mf\rangle
=\frac12\sum_{x,y}(K_s)_{xy}|f_x-f_y|^2\ge0.
\]

Thus the same finite generator supports at least

\[
U_t=e^{-itA_M}
\quad\text{and}\quad
P_t=e^{-tA_M}.
\]

The first is unitary on amplitudes; the second is a symmetric stochastic semigroup on populations. Nothing in the real spectral data selects one temporal category. **[Proven]**

## 2. Kernel-split boundary

This monograph studies the strict operational kernel only. It does not identify the strict kernel with the legacy ontological kernel, transfer legacy physical roles, or close the legacy-to-strict bridge. The two kernels may be related through an explicit future completion theorem, but no such theorem is assumed here.

## 3. Selector and unit boundaries

All operators used here are reflection even unless an odd diagnostic or apparatus polarity is explicitly introduced. Therefore no calculation below discharges the strict selector obstruction. Likewise, every finite number computed from the dimensionless matrices is scale weight zero. No physical length, time, mass, energy, action unit, or value of \(\hbar\) follows.

## 4. Reproducibility

The primary computation is

*fin_programs_21_30_experiments.py*.

It writes

*FIN_Programs_21_30_Results.json*

and generates the figures used in this preprint. NumPy and SciPy perform the matrix calculations. Randomized studies use the fixed seed 20260719 or a program-specific deterministic offset. Synthetic results are therefore exactly replayable up to platform-level floating-point differences.

---

# Part II — Program 21: intervention-aware process-tensor tomography

## 1. Question

Can intermediate interventions distinguish coherent, classical Markov, dephasing, and finite-memory dynamics that share the strict spatial generator?

## 2. Declared candidate family

Four finite models were tested.

- **Coherent:** \(\rho\mapsto U_t\rho U_t^\dagger\).
- **Markov:** \(p\mapsto P_tp\).
- **Site dephasing:** the Lindblad equation

\[
\dot\rho=-i[A_M,\rho]+\gamma(\operatorname{diag}\rho-\rho),
\qquad \gamma=0.55.
\]

- **Finite classical memory bath:** a hidden branch retains one of two dephasing rates, \(\gamma=0.08\) or \(1.25\), through all time slots.

The intervention signature contains: no intervention, an unread position measurement, the full two-time position record, the full three-time record, a terminal record, and three reset preparations. This is an instrument-restricted process signature.

## 3. Optimized finite design

A bounded scan selected

\[
(t_1,t_2,t_3)=(1.16875,2.125,3.50625).
\]

Eight normalized probability blocks produced a 1,944-component signature. After centering the four model rows, the signature matrix had rank three, the maximum possible affine rank for four distinct models. The smallest pairwise Jensen–Shannon divergence was

\[
4.17677\times10^{-3}\ \text{nats},
\]

and the smallest total-variation distance was

\[
0.0697765.
\]

In 80 simulated trials per model with 400 shots per probability block, maximum likelihood classified all 320 trials correctly. **[Strong evidence within the declared four-model family]**

FINFIGUREP21

## 4. Why this is not full tomography

An unrestricted linear operation on a \(d\)-level density matrix lives in a real space of dimension \(d^4\). For \(d=12\),

\[
d^4=20,736.
\]

Two unrestricted intervention slots scale as

\[
d^8=429,981,696,
\]

before a final informationally complete measurement is counted. Therefore the word “tomography” must be qualified. The executed protocol identifies a small candidate family under a restricted instrument menu; it does not reconstruct a generic twelve-mode process tensor. **[Proven]**

## 5. Falsification tests

- Removing interventions collapses coherent and final-dephased population models in several preparations.
- Replacing the hidden branch by a reset branch removes the memory signature.
- Allowing arbitrary readout maps can send every model to the uniform distribution.
- Increasing the candidate family without increasing the instrument span can produce likelihood-equivalent processes.

## 6. Verdict

- Targeted intervention-aware discrimination: **[Strong evidence]**.
- Modest unrestricted informationally complete tomography: **[Refuted]**.
- A physical process tensor derived from \(A_M\) alone: **[Refuted]**.
- Recommended next move: reduce process-tensor complexity only through an independently justified symmetry, memory-order, or low-rank theorem.

The operational Markov condition and process-tensor framework follow the intervention-based formulation of [Pollock et al.](https://doi.org/10.1103/PhysRevLett.120.040405) and the complete framework of [Pollock et al.](https://doi.org/10.1103/PhysRevA.97.012127).

---

# Part III — Program 22: positivity-preserving refinement

## 1. Obstruction in the naive continuation

The frozen analytic profile

\[
k(r)=\frac{\cos(0.18575r+0.16250)}{1+r^{1.8}}
\]

is positive for integer distances \(1\le r\le6\), but first becomes negative at \(r=8\). Consequently, the naive \(N=16\) extension is not an off-diagonally nonnegative Markov kernel. **[Proven]**

## 2. Construction A: fixed-circumference nonlocal refinement

Define

\[
(W_N^{\rm circ})_{xy}
=\frac{12}{N}
k\!\left(\frac{12d_N(x,y)}{N}\right),
\qquad x\ne y.
\]

Then \(W_{12}^{\rm circ}=K_s\) exactly. Since \(0<12d_N/N\le6\), all off-diagonal entries are positive at every \(N\). The associated form converges by Riemann sums to

\[
\mathcal E_{\rm circ}(f)
=6\int_0^1\!\int_0^1
k(12d_{\mathbb T}(x,y))
|f(x)-f(y)|^2\,dx\,dy.
\]

The limiting jump kernel is positive at every separation. This is a well-defined nonlocal Dirichlet form, not a local Laplacian. **[Proven for Fourier modes; strong evidence for full form convergence]**

## 3. Construction B: anchored exponential tail

Retain the six strict anchor values exactly and define

\[
w(d)=
\begin{cases}
k(d),&1\le d\le6,\\
k(6)q^{d-6},&d>6,
\end{cases}
\qquad
q=\frac{k(6)}{k(5)}=0.458774\ldots .
\]

This family is positive and summable. Its second moment is

\[
D_2=\sum_{d\ge1}w(d)d^2=4.40795255\ldots .
\]

For fixed low Fourier mode,

\[
N^2\lambda_{N,k}
\longrightarrow
4\pi^2D_2k^2.
\]

At \(N=512\), the scaled first gap differed from the predicted limit by \(2.70\times10^{-4}\). This construction therefore gives a local diffusive continuum after \(N^2\) scaling. **[Proven by dominated expansion; numerically verified]**

FINFIGUREP22

## 4. Smooth completely monotone approximation

A nonnegative exponential-mixture fit

\[
w(d)=\sum_jc_je^{-\lambda_jd},\qquad c_j\ge0,
\]

used two active terms and reproduced the six anchors with relative residual \(2.805\times10^{-3}\). This supplies a smooth positivity-preserving family, but it no longer matches the frozen row exactly.

## 5. Nonuniqueness theorem

Constructions A and B agree with the complete \(N=12\) strict row and preserve positivity, but their continuum limits are inequivalent: A is nonlocal at fixed circumference; B is locally diffusive after \(N^2\) scaling. Therefore finite matching plus positivity does not select refinement, locality, scaling exponent, or continuum geometry. **[Proven by counterexample]**

## 6. Verdict

- Existence of a positivity-preserving refinement: **[Proven]**.
- Uniqueness or FIN-internal selection: **[Refuted]**.
- A local continuum without an added scaling rule: **[Refuted]**.
- Best next theorem: add an independently motivated refinement functor, locality criterion, and normalization, then test Mosco convergence rather than only low modes.

Mosco convergence is the correct operator-form topology because it controls the associated semigroups and resolvents; see, for example, [Kuroda's graph-limit study](https://arxiv.org/abs/1106.5476).

---

# Part IV — Program 23: doubled Dirac and hyperbolic refinement

## 1. Local incidence construction

Let \(B_N\) be the oriented vertex-edge incidence matrix of the \(N\)-cycle and define

\[
D_N=
\begin{pmatrix}
0&B_N^\top\\
B_N&0
\end{pmatrix}.
\]

Then

\[
D_N^2=
\begin{pmatrix}
B_N^\top B_N&0\\
0&B_NB_N^\top
\end{pmatrix}.
\]

The square residual was exactly zero at every tested size. The positive dispersion is

\[
\omega_N(k)=2\sin\frac{\pi|k|}{N},
\]

and the fitted finite-size exponent was \(z=0.999082\). **[Proven analytically; numerically verified]**

FINFIGUREP23

The probability outside the conservative radius \(\lceil2t+2\rceil\) was below \(6.4\times10^{-12}\) in the tested local propagations. This is a finite locality diagnostic, not an exact relativistic light cone.

## 2. Exact strict square root

For every unordered pair \(i<j\), introduce an edge with incidence entries \(\pm\sqrt{(K_s)_{ij}}\). The resulting weighted incidence matrix \(B_W\) satisfies

\[
B_W^\top B_W=A_M.
\]

Hence

\[
D_W=
\begin{pmatrix}
0&B_W^\top\\
B_W&0
\end{pmatrix}
\]

is an exact doubled square root of the strict generator. **[Proven]**

Because every off-diagonal strict weight is positive, \(B_W\) contains all 66 unordered pairs. It is a complete-graph incidence operator and therefore nonlocal in the cycle metric. Its effective square-root exponents from the first mode pairs are approximately 0.532 and 0.435, not a stable \(z=1\) law. **[Strong evidence]**

## 3. Orientation obstruction

Changing edge orientation multiplies rows of \(B\) by signs. This conjugates the doubled Dirac operator by a block sign unitary and leaves \(D^2\) invariant. Likewise, \(D\) and \(-D\) have the same square. The construction therefore supplies no canonical spatial or temporal orientation. **[Proven]**

## 4. Verdict

- Local \(z=1\) cycle Dirac refinement: **[Proven]**.
- Exact square root of the strict generator: **[Proven]**.
- Exact strict matching plus cycle locality plus stable \(z=1\): **[Not obtained]**.
- Lorentz symmetry: **[Speculative]**; it would additionally require a clock/length calibration, causal cone, continuum representation, spin structure, and symmetry theorem.

The incidence Dirac identity is standard graph operator theory; a modern treatment is [Casiday et al.](https://doi.org/10.1080/03081087.2022.2158297). The FIN-specific result is the explicit locality-versus-exactness obstruction for this frozen weighted operator.

---

# Part V — Program 24: Osterwalder–Schrader positivity

## 1. Free massive temporal completion

Introduce an explicit Euclidean time coordinate and the free Hessian

\[
H_T=L_T\otimes I+I\otimes(A_M+m^2I),
\qquad m^2=0.3.
\]

For even temporal size \(T\), reflect \(t\mapsto-t\bmod T\), and restrict test functions to \(t=1,\ldots,T/2-1\). The reflection matrix is the corresponding block of \(C_T=H_T^{-1}\).

For \(T=8,12,16,24\), no eigenvalue below numerical tolerance was found. At \(T=16\), the smallest value was of order \(10^{-16}\). This agrees with the exact transfer-kernel factorization of a free massive covariance. **[Proven in the declared finite free theory]**

The infinite-time sampled covariance

\[
C(\tau)=\frac{e^{-|\tau|\Omega}}{2\Omega},
\qquad
\Omega=(A_M+m^2I)^{1/2},
\]

has reflection form

\[
C(t_i+t_j)
=e^{-t_i\Omega}(2\Omega)^{-1}e^{-t_j\Omega},
\]

which is a Gram matrix and hence positive. **[Proven]**

## 2. Higher-derivative counterexample

Replace the time Hessian by

\[
L_T+\alpha L_T^2,
\qquad \alpha>0.
\]

The full quadratic Hessian remains positive, yet at \(T=16\) the reflection form has minimum eigenvalues approximately

\[
-3.70\times10^{-7},\quad
-2.23\times10^{-4},\quad
-1.83\times10^{-2}
\]

for \(\alpha=0.01,0.1,1\). Thus a stable Euclidean quadratic action need not reconstruct a positive Hilbert-space theory. **[Proven by finite counterexample]**

## 3. Scope

The result constructs a conditional free temporal extension. It does not derive the time lattice, select a physical mass, establish an interacting continuum, or show that FIN itself satisfies all Osterwalder–Schrader axioms. The original reconstruction theorem requires a full family of Euclidean functions, covariance, symmetry, regularity, and positivity; see [Osterwalder and Schrader](https://projecteuclid.org/euclid.cmp/1103858969).

## 4. Verdict

- Free reflection-positive completion: **[Proven]**.
- Positivity of the Euclidean Hessian as sufficient: **[Refuted]**.
- FIN-selected temporal/QFT completion: **[Speculative]**.
- Physical Hamiltonian from analytic continuation alone: **[Refuted]**.

---

# Part VI — Program 25: SPAM- and clock-robust Bayesian identification

## 1. Nuisance model

For temporal model \(m\), use

\[
p_m(x|t,\gamma,\varepsilon)
=(1-\varepsilon)p_m^0(x|\gamma t)+\frac{\varepsilon}{12},
\]

with calibrated ranges

\[
0.95\le\gamma\le1.05,
\qquad
0\le\varepsilon\le0.06.
\]

The candidate temporal categories are strict coherent evolution and strict diffusion.

## 2. Maximin design

A greedy nuisance-profiled design selected

\[
t=(0.615385,0.737179,0.493590).
\]

The worst profiled Jensen–Shannon separation remained above \(7.7\times10^{-2}\) nats. With 600 shots per time and 120 trials per true model, profile likelihood classified all 240 trials correctly. **[Strong evidence inside the declared nuisance box]**

The optimized single snapshot remained distinguishable even under a broader clock/SPAM grid in this restricted two-model comparison. Therefore “three times are universally necessary” is false. Multiple times instead protect against misspecification and enlarge the alternative class that can be rejected. **[Refuted as a universal necessity claim]**

## 3. Exact nonidentifiability examples

Final dephasing is population invisible:

\[
\operatorname{diag}(U_t\rho U_t^\dagger)
=\operatorname{diag}\Delta(U_t\rho U_t^\dagger).
\]

No number of terminal population snapshots separates those two channels. Likewise, an unrestricted readout channel may map both temporal models to the uniform distribution. Calibration assumptions are therefore logically necessary. **[Proven]**

## 4. Verdict

- Coherent-versus-diffusive identification under bounded nuisance: **[Strong evidence]**.
- Unrestricted temporal-category identification from populations: **[Refuted]**.
- Universal three-time necessity: **[Refuted]**.
- Hardware performance: **[Not tested]**.

The treatment of SPAM parameters as inferred quantities follows the logic of self-consistent Bayesian device characterization; compare [Landa et al.](https://doi.org/10.1103/PhysRevResearch.4.013199). Imperfect clocks are themselves physical resources rather than harmless labels; see [Xuereb et al.](https://doi.org/10.1103/PhysRevLett.131.160204).

---

# Part VII — Program 26: minimal environment and dilation identifiability

## 1. Two population-equivalent quantum semigroups

Let \(p_d(t)=(P_t)_{d0}\) and let \(S\) be the cyclic shift. Define

\[
\mathcal R_t(\rho)
=\sum_{d=0}^{11}p_d(t)S^d\rho S^{-d}.
\]

Define a second semigroup with edge jumps

\[
L_{ij}=\sqrt{(K_s)_{ij}}|i\rangle\langle j|.
\]

It has the closed form

\[
\mathcal E_t(\rho)
=e^{-st}(\rho-\operatorname{diag}\rho)
+\operatorname{diag}(P_t\operatorname{diag}\rho).
\]

For every input state,

\[
\operatorname{diag}\mathcal R_t(\rho)
=\operatorname{diag}\mathcal E_t(\rho)
=P_t\operatorname{diag}\rho.
\]

Every population preparation, position instrument, and position record is therefore identical in the two models. **[Proven]**

At \(t=1\), the population residual was below \(4\times10^{-16}\), while the Choi ranks were 12 and 144. Since minimal one-time Stinespring dimension equals Choi rank, the same complete population process admits radically different environmental sizes. For the coherent uniform input, the channel-output trace distance was

\[
\frac{11}{12}(1-e^{-s})=0.742426\ldots .
\]

One coherent preparation separates the models, but populations never do. **[Proven]**

## 2. Finite realization rank

A three-frequency finite bath generated a coherence signal whose noiseless Hankel matrix had exactly three nonzero singular values. With complex noise \(2\times10^{-3}\), the remaining singular values rose to order \(10^{-2}\), making model order threshold-dependent. **[Strong evidence]**

## 3. Finite-window nonuniqueness

Two positive spectral measures were constructed to agree at seven training times to \(8.3\times10^{-17}\) while differing later by 0.230456. Thus even exact finite-window coherence data do not select a unique spectral density. **[Proven by null-space construction]**

## 4. Verdict

- One-time Choi rank under informationally complete control: **[Identifiable]**.
- Unique autonomous bath from population records: **[Refuted]**.
- Finite-window realization rank under a declared exponential model: **[Strong evidence]**.
- Unique infinite dilation: **[Refuted]**.

This distinction is consistent with realization-theoretic dimension estimation: dimension bounds require a fixed interaction model and sufficient probe access, as emphasized by [Sone and Cappellaro](https://arxiv.org/abs/1702.03280).

---

# Part VIII — Program 27: relational-observable rewrite audit

## 1. Complete frame invariant

A frame is \(f=(r,\lambda)\in\mathbb Z_{12}\times\{\pm1\}\). Under

\[
(a,\varepsilon)\cdot(r,\lambda)
=(a+\varepsilon r,\varepsilon\lambda),
\]

the system–apparatus pair has invariant

\[
(\Delta,c)
=\bigl(\lambda_a(r_s-r_a),\lambda_s\lambda_a\bigr).
\]

All 576 frame pairs and 24 dihedral transformations were enumerated. The diagonal action has 24 orbits, exactly the number of \((\Delta,c)\) values, and zero invariant failures occurred across 13,824 transformed cases. Therefore \((\Delta,c)\) is complete. **[Proven]**

## 2. Dynamical covariance

For

\[
p_{\rm obs}(\delta|f_s,f_a)
=p_{\rm base}(r_a+\lambda_a\delta-r_s),
\]

maximum covariance residuals were below \(2.5\times10^{-16}\) for coherent dynamics and below \(4.2\times10^{-17}\) for diffusion. **[Strong numerical certificate]**

## 3. What cannot be rewritten away

Absolute site-zero probability varied by 0.192 in a generic transformed-state audit. A signed current was reflection odd: treating it as invariant produced residual 0.870, whereas its odd covariance residual was \(5.6\times10^{-17}\). Multiplying by apparatus polarity makes it relational, but the polarity is then a physical resource of the apparatus.

The audit classified twelve operational objects:

- eight admit invariant or relational representations;
- absolute origin requires a reference;
- physical clock reading requires calibration;
- absolute orientation requires an odd resource;
- persistent record identity requires a memory-bearing apparatus.

## 4. Verdict

- Completeness of \((\Delta,c)\) for frame pairs: **[Proven]**.
- Relational rewrite of tested operational observables: **[Proven in the declared registry]**.
- Elimination of an absolute strict-core selector requirement: **[Refuted]**.
- Repository-wide historical role-transfer audit: **[Not attempted]**.

The mathematical setting is the resource theory of reference frames and asymmetry; compare [Gour and Spekkens](https://arxiv.org/abs/0711.0043) and the mode decomposition of [Marvian and Spekkens](https://arxiv.org/abs/1312.0680).

---

# Part IX — Program 28: adaptive integrability and odd-sector bifurcation

## 1. Concrete covariance map

For the localized state \(|0\rangle\), define

\[
C_T(K)=\frac1T\int_0^T
\operatorname{Re}\left(e^{itK}|0\rangle\langle0|e^{-itK}\right)dt.
\]

The one-form

\[
\alpha_K(H)=\operatorname{Tr}(C_T(K)H)
\]

is locally exact only if its Jacobian is self-adjoint in the Frobenius pairing.

## 2. Jacobian calculation

The real symmetric zero-diagonal tangent space has dimension 66. Reflection decomposes it into a 36-dimensional even block and a 30-dimensional odd block. A centered finite-difference derivative at \(K_s\), with \(T=4\), gave

\[
\frac{\|DC_T-DC_T^*\|_F}{\|DC_T\|_F}=1.66281.
\]

The even–odd mixing blocks were below \(1.5\times10^{-9}\), confirming equivariance. The leading odd eigenvalue was

\[
0.161035+0.413525i.
\]

Thus a centered linear family \(\mu DC_T-I\) reaches zero real part near \(\mu=6.21\), but the nonzero imaginary part signals an oscillatory instability, not a real pitchfork. **[Strong evidence]**

FINFIGUREP28

Independent scans at \(T=0.5,2,5,10\) kept the integrability defect near 1.69–1.73, and finite-difference refinement left it stable. The nonintegrability is not a step-size artifact. **[Strong evidence]**

## 3. Fixed-point and degeneracy obstructions

For unitary self-driving,

\[
\operatorname{Tr}(K C_T(K))
=\langle0,K|0\rangle.
\]

Since \((K_s)_{00}=0\), the best scalar relation \(C_T(K_s)\approx\gamma K_s\) has \(\gamma\approx0\) and relative residual one. The strict kernel is not a fixed point of this scalar-damped self-covariance law. **[Proven for the declared update class]**

The infinite-time dephasing map is discontinuous at the degenerate circulant spectrum along a tested odd perturbation, with a nonvanishing limiting jump near 0.2635. Therefore a Fréchet derivative for a center-manifold calculation is ill posed there unless degeneracies are resolved or finite averaging is retained. **[Strong evidence]**

## 4. Verdict

- Automatic feedback-gradient structure: **[Refuted for \(C_T\)]**.
- Strict scalar-damped self-covariance fixed point: **[Refuted]**.
- Nonzero internal odd source: **[Not found]**.
- FIN-specific pitchfork selector: **[Not established]**.
- A post-hoc centered flow can be engineered, but it is not a derivation.

---

# Part X — Program 29: machine-checked operational no-go library

## 1. Executed certificates

Five minimal theorem schemas were encoded as exact or finite executable checks.

1. **Selector no-section.** A reflection-fixed input cannot map equivariantly to a sign torsor flipped by reflection.
2. **Positive-scale no-section.** Scale-weight-zero data cannot select a nonzero element of a free weight-one positive torsor.
3. **Unitary–stochastic intersection.** A nonnegative unitary stochastic matrix is a permutation matrix.
4. **Finite-bath obstruction.** A finite autonomous unitary bath cannot equal strictly exponential irreversible decay for all positive times.
5. **Dynamical nonselection.** The same positive graph generator admits both \(e^{-itA}\) and \(e^{-tA}\) when no temporal-category axiom is supplied.

All five executable assertions passed. The certificate payload has SHA-256

**964822117316b107419908cfc2373bedb53119678714f1c0848ce0d6a1b54e29**.

## 2. Hidden-premise audit

Each theorem fails if its defining premise is removed. An odd source evades the selector theorem; a scale-charged source evades the scale theorem; negative or complex entries evade the stochastic intersection; an infinite or reset bath evades finite recurrence; an operational temporal axiom evades dynamical nonselection. These are new premises, not contradictions.

## 3. Formalization boundary

The mathematical proofs and exact Python certificates pass. The installed *lean* command was only an *elan* launcher without a configured toolchain, so an independent Lean or Isabelle kernel replay was not completed. Calling the present library “proof-assistant certified” would be false.

## 4. Verdict

- Mathematical theorem content: **[Proven]**.
- Executable finite certificates: **[Proven]**.
- Independent proof-kernel replay: **[Incomplete]**.
- Next step: port one theorem at a time to a pinned Lean version and archive the build manifest.

---

# Part XI — Program 30: blinded twelve-mode held-out challenge

## 1. Preregistered synthetic protocol

The analysis used a positive-circulant alternative library containing 100 profiles with log-weight perturbation scales

\[
0.05,0.10,0.20,0.40,0.70.
\]

Only an overall rate \(\alpha\) and uniform SPAM parameter \(\varepsilon\) were fitted.

Training data contained spectral modes 1 and 2, diffusion at \(t=0.35\), and unitary populations at \(t=0.55\). Held-out data contained modes 3–6, diffusion at \(t=0.9,1.7\), and unitary populations at \(t=1.1,2.2\). Each population distribution used 4,000 shots; spectral Gaussian noise was 0.012.

The hidden label sequence was committed by SHA-256 before scoring. This is a deterministic software commitment, not an externally timestamped preregistration.

## 2. Results

Across 30 balanced trials, the combined spectral-and-multitime protocol classified 24 correctly:

\[
\text{accuracy}=0.80,
\qquad
\text{FIN sensitivity}=0.80,
\qquad
\text{alternative specificity}=0.80.
\]

The confusion matrix, with true classes in rows, was

\[
\begin{pmatrix}
12&3\\
3&12
\end{pmatrix}.
\]

Median held-out \(\Delta\mathrm{NLL}=\mathrm{NLL}_{\rm generic}-\mathrm{NLL}_{\rm FIN}\) was \(+20.409\) on FIN data and \(-72.724\) on alternative data. Specificity was \(2/3\) at perturbation scales 0.05, 0.10, and 0.20, and \(1\) at 0.40 and 0.70. Each scale had only three alternative trials, so these scale-resolved estimates have large sampling uncertainty. **[Moderate evidence]**

An exploratory spectral-only affine challenge against a particularly close two-exponential positive alternative was near chance. It was not the preregistered primary analysis. This negative control demonstrates that success depends materially on the declared alternative library and on multitime data.

## 3. Falsification result

The correct conclusion is a resolution boundary:

> The present synthetic data reject many moderate alternatives but do not establish uniqueness against near-FIN positive circulant profiles.

The test also does not validate hardware, dimensional calibration, or fundamental physics. Its purpose is to freeze a workflow that can later be applied to externally generated data.

## 4. Verdict

- Synthetic held-out pipeline: **[Strong evidence]**.
- Robustness against moderate alternatives: **[Moderate evidence]**.
- Uniqueness against near profiles: **[Refuted at the present resolution]**.
- Hardware challenge: **[Not executed]**.
- Fundamental-physics inference: **[None]**.

The next legitimate experiment is an independently fabricated twelve-mode device with raw counts, fabrication tolerances, randomized acquisition order, independent calibration, and a generic-model comparison frozen before unblinding.

---

# Part XII — Cross-cutting study: from the strict kernel to an inverse action

## 1. Three operators that must not be conflated

Use distinct notation:

- \(K_s\): the zero-diagonal strict kernel used as a Hamiltonian or weighted adjacency;
- \(A_M=sI-K_s\): the positive Markov/Dirichlet generator;
- \(Q_Y=L+m_Y^2I\): a candidate Gaussian precision.

These operators may share a Fourier eigenbasis, but they have different mathematical roles.

## 2. Exact representability criterion

Let \(L\) be the cycle Laplacian,

\[
G_m=(L+m^2I)^{-1},
\qquad
\mathcal N(G_m)=G_m-\operatorname{diag}G_m.
\]

On a vertex-transitive graph, \(\operatorname{diag}G_m=g_0I\). Therefore

\[
K=a\mathcal N(G_m)
\]

holds exactly if and only if \(K\) is scalar on every \(L\)-eigenspace and, in a common spectral basis,

\[
\kappa_r=\frac{a}{\mu_r+m^2}-ag_0.
\]

Equivalently, after a contact value is chosen, \((\kappa_r-c_0)^{-1}\) must be affine in \(\mu_r\). **[Proven]**

The strict kernel commutes with \(L\) to machine precision, so the functional-calculus prerequisite is exact. The one-pole spectral law is not exact.

## 3. Recomputed Yukawa fit

Equal weighting of the seven distinct Fourier symbols gives approximately

\[
a_Y=2.05275945,
\qquad
m_Y^2=0.74918143,
\qquad
c_0=-1.08957768.
\]

The relative symbol residual is 0.024916, the full-matrix relative residual is 0.028675, and the off-diagonal correlation is 0.999339. Relative real-space errors grow through the tail and reach approximately 27% at distance six.

A Frobenius-optimal contact-subtracted fit instead gives

\[
a_Y=2.04199291,
\qquad
m_Y^2=0.74899211,
\]

with relative matrix residual 0.028187. The fitted parameters therefore depend slightly on the declared loss. **[Strong numerical evidence]**

FINFIGUREACTION

The scientifically correct statement is:

> \(K_s\) is well approximated by a contact-subtracted discrete screened-Poisson resolvent on \(C_{12}\).

Exact equality is refuted. “Yukawa” here refers to a discrete screened Green form; it is not a three-dimensional Yukawa tail and not a fermion–Higgs Yukawa coupling.

## 4. The conditional Gaussian action

Let

\[
Q_Y=L+m_Y^2I.
\]

Then

\[
S_G[\phi;J]
=\frac{1}{2a_Y}\phi^\top Q_Y\phi-J^\top\phi
\]

has Euler–Lagrange equation

\[
Q_Y\phi=a_YJ
\]

and, after specifying a unit Gaussian measure, covariance

\[
G=a_YQ_Y^{-1}.
\]

Its contact-subtracted covariance approximates \(K_s\). If the coefficient is written as \(1/2\) rather than \(1/(2a_Y)\), the amplitude must be absorbed into a field or source normalization. **[Proven]**

Classical stationarity produces a response \(\phi=GJ\). Calling \(G\) a two-point function additionally assumes a Gaussian measure and a normalization convention. No physical \(\hbar\) is generated.

## 5. Why the zero-diagonal kernel is not a covariance

The contact-subtracted kernel has trace zero. A nonzero positive-semidefinite Hermitian matrix has positive trace, so a nonzero zero-trace kernel must be indefinite. The strict matrix has seven negative eigenvalues. **[Proven]**

Therefore the stable Gaussian covariance is \(G\), not \(K_s=\mathcal N(G)+R\). The term “normal ordering” should be read here as contact/diagonal subtraction unless a vacuum and Wick ordering are independently defined.

## 6. The proposed two-sector action

Consider

\[
S[\phi,K]
=\frac12\phi^\top(L+m^2I)\phi-J^\top\phi
+\operatorname{tr}V(K)-\operatorname{tr}(KC_d).
\]

For fixed drive \(C_d\), stationarity gives

\[
(L+m^2I)\phi=J,
\qquad
V'(K)=C_d
\]

on the unconstrained symmetric space. The mixed \(\phi\)-\(K\) Hessian is exactly zero. Thus the written functional is a direct sum: the Gaussian sector and the operator sector are juxtaposed, not dynamically linked. **[Proven]**

If \(C=C(\phi)\), the field equation acquires an adjoint derivative. If \(C=C(K)\), differentiating \(\operatorname{tr}[KC(K)]\) adds the omitted term \(DC(K)^*[K]\). The simple equation \(V'(K)=C(K)\) is not automatically variational. **[Proven]**

## 7. Degree fourteen: what is and is not minimal

Let \(\kappa_0,\ldots,\kappa_6\) be the seven distinct strict eigenvalues and

\[
p(x)=\prod_{r=0}^6(x-\kappa_r),
\qquad V(x)=p(x)^2.
\]

Every \(\kappa_r\) is a nondegenerate scalar minimum. Seven prescribed scalar wells require a polynomial of degree at least fourteen: between seven minima there must be at least six additional critical points, so \(\deg V'\ge13\). **[Proven under the seven-well axiom]**

This does not make \(K_s\) the unique matrix minimum. \(\operatorname{tr}V(K)\) is constant on conjugacy orbits and admits multiple root assignments. With \(C_d=0\), it cannot select an eigenbasis. With nonzero \(C_d\), \(V'(K_s)=0\) means the unconstrained gradient at \(K_s\) is \(-C_d\), not zero.

Unconditionally, a degree-two action already selects the known target:

\[
F_*(K)=\frac12\|K-K_s\|_F^2
=\frac12\operatorname{tr}K^2-\operatorname{tr}(KK_s)+\text{const}.
\]

It is unique but completely tautological because the drive stores the answer. Therefore “degree fourteen is the smallest possible action” is false unless the seven-well restriction is stated explicitly. **[Refuted as an unconditional claim]**

## 8. Exact covariance-level stationary reconstruction

The inverse propagator problem has a cleaner variational solution. For \(G>0\), define

\[
\Gamma_Q(G)
=\frac12\operatorname{tr}(QG)-\frac12\log\det G.
\]

Then

\[
\delta\Gamma_Q=0
\quad\Longleftrightarrow\quad
Q-G^{-1}=0
\quad\Longleftrightarrow\quad
G=Q^{-1}.
\]

**[Proven]** This is the finite Gaussian/2PI variational structure behind the slogan “reconstruct the action from the propagator.”

To bind the normal-ordered kernel explicitly, introduce a multiplier:

\[
\Gamma[G,K,\Lambda]
=\frac12\operatorname{tr}(QG)
-\frac12\log\det G
+\langle\Lambda,K-\mathcal N(G)\rangle.
\]

Its stationary equations give

\[
G=Q^{-1},
\qquad
K=\mathcal N(G),
\qquad
\Lambda=0.
\]

For the empirical strict kernel, the equality must be replaced by a residual or penalty because the screened one-pole fit has about 2.8% matrix error.

## 9. Exact but nonlocal inverse actions

Because the minimum strict eigenvalue is approximately \(-0.681875\), every \(q>0.681875\) makes

\[
G_q=K_s+qI
\]

positive. The Gaussian precision \(Q_q=G_q^{-1}\) then reproduces \(K_s\) exactly after contact subtraction. This gives infinitely many exact inverse actions, parameterized by the lost coincidence value \(q\), and they are generally nonlocal. **[Proven]**

At a propagator-optimal contact choice, the exact precision differs from the nearest-neighbour form \(L+m^2I\) by about 3.4% in Frobenius norm. Optimizing precision locality lowers that residual but worsens the forward kernel fit. The “best inverse action” therefore depends on which space and loss function are declared. **[Strong evidence]**

## 10. Final verdict on inverse reconstruction

| Claim | Verdict |
|---|---|
| exact nearest-neighbour Yukawa identity | **[Refuted]** |
| high-correlation discrete screened-Green approximation | **[Strong evidence]** |
| conditional local Gaussian action for the fitted component | **[Proven]** |
| exact inverse action after contact completion | **[Proven but nonunique]** |
| strict kernel itself as a stable Gaussian covariance | **[Refuted]** |
| degree-14 minimality for seven scalar wells | **[Proven conditionally]** |
| unconditional degree-14 minimality | **[Refuted]** |
| noncircular selection by fitted \(V,C_d\) | **[Refuted]** |
| physical Lagrangian, units, clock, or selector closure | **[Unsupported]** |

The central distinction is between two inverse problems:

1. **Propagator to field action.** After a contact completion, positivity, and a locality choice, invert the covariance. This is standard and conditionally rigorous.
2. **Observed kernel to a law selecting that kernel.** This is infinitely nonunique. A potential and drive fitted after seeing \(K_s\) encode the target rather than predict it.

The finite covariance functional is the cleanest mathematical answer to the first problem. No result in this monograph solves the second without an additional independent principle.

The covariance-level variational construction is the finite analogue of the two-particle-irreducible effective-action formalism introduced by [Cornwall, Jackiw and Tomboulis](https://doi.org/10.1103/PhysRevD.10.2428).

---

# Part XIII — Integrated theorem and obstruction inventory

## 1. Newly established theorems

### Theorem A — positive-refinement nonuniqueness

The same positive twelve-mode strict row admits a nonlocal fixed-circumference continuum and a local diffusive continuum. Finite matching and positivity cannot select the continuum category. **[Proven]**

### Theorem B — Dirac locality–exactness trade-off

The nearest-neighbour incidence Dirac operator is local and \(z=1\) but does not square to the strict generator. The exact strict weighted incidence square root uses all 66 edges and is nonlocal in the cycle metric. **[Proven for the two declared constructions]**

### Theorem C — free OS-positive completion

For \(m^2>0\), the covariance \(e^{-|t|\sqrt{A_M+m^2I}}/[2\sqrt{A_M+m^2I}]\) is reflection positive. **[Proven]**

### Theorem D — population-dilation nonidentifiability

Complete population records of the strict Markov semigroup do not determine a quantum dilation or one-time Choi rank. **[Proven by ranks 12 versus 144]**

### Theorem E — complete relational frame invariant

The pair \((\Delta,c)\) completely classifies the diagonal dihedral action on system–apparatus frames. **[Proven]**

### Theorem F — contact-subtracted covariance obstruction

A nonzero zero-trace Hermitian matrix cannot be positive semidefinite. The strict normal-ordered kernel is not a Gaussian covariance. **[Proven]**

### Theorem G — covariance-level inverse action

The stationary point of \(\Gamma_Q(G)=\tfrac12\operatorname{tr}(QG)-\tfrac12\log\det G\) is exactly \(G=Q^{-1}\). **[Proven]**

## 2. Refuted claims

- The twelve-mode spectrum determines a unique continuum. **[Refuted]**
- Unitarity alone implies \(z=1\) or Lorentz symmetry. **[Refuted]**
- Positive Euclidean Hessian implies reflection positivity. **[Refuted]**
- Population records determine a unique bath. **[Refuted]**
- A fitted feedback covariance is automatically a gradient. **[Refuted]**
- A generic two-slot process tensor can be reconstructed with a modest unrestricted instrument set. **[Refuted]**
- The strict kernel is exactly a nearest-neighbour Yukawa propagator. **[Refuted]**
- The direct-sum action couples the field and kernel sectors. **[Refuted]**
- Degree fourteen is unconditionally the smallest action selecting \(K_s\). **[Refuted]**

## 3. Unresolved physical bridges

None of Programs 21–30 derives the following from the strict spectral core:

- a state/probability axiom;
- a temporal category;
- a dimensionful clock;
- a physical preparation and instrument;
- an apparatus and persistent record;
- a bath or environmental law;
- a selector or nonzero odd resource;
- a scale-charged unit source;
- a local Lorentzian refinement;
- an interacting reflection-positive continuum;
- an independently motivated empirical domain.

These are typed missing structures, not one unnamed scalar.

---

# Part XIV — Research recommendations after execution

## 1. Highest-probability next studies

1. **Externally preregistered multitime emulator challenge.** Freeze raw data format, nuisance priors, generic alternatives, and acceptance thresholds before receiving the data. Probability of a technically decisive emulator result: 0.70.
2. **Low-rank process-tensor theorem.** Prove or refute a symmetry/Markov-order reduction that lowers the \(d^8\) two-slot burden without assuming the answer. Probability: 0.55.
3. **Pinned proof-assistant release.** Formalize the five no-go theorems in a fixed Lean toolchain and archive the build hash. Probability: 0.90.
4. **Mosco competition between refinements.** Compare the local exponential-tail and nonlocal fixed-circumference constructions using the same topology and held-out sizes. Probability: 0.80.
5. **Interacting OS gate.** Add one declared onsite interaction and prove reflection positivity or exhibit the first negative form. Probability: 0.35.

## 2. High-risk, high-impact studies

6. **Strict-local Dirac compromise theorem.** Optimize sparsity under controlled operator error and determine whether a stable \(z=1\) window exists. Probability: 0.35.
7. **Independent covariance law.** Specify \(C(K)\) before seeing \(K_s\), test its Helmholtz integrability conditions, and demand held-out-size prediction. Probability: 0.25.
8. **Exact inverse-action locality theorem.** Bound the decay of \((K_s+qI)^{-1}\) and identify the contact value minimizing nonlocal precision tails. Probability: 0.75.
9. **Operational odd-resource experiment.** Introduce an explicitly measured chiral apparatus variable and verify relational covariance without claiming an absolute selector. Probability: 0.65.
10. **Scale-bearing conversion experiment.** Treat the conversion scales as explicit calibration axioms and test conditional predictions; do not relabel them strict outputs. Probability: 0.60.

## 3. Stop rules

- Do not repeat invariant scalar searches for an absolute orientation.
- Do not infer units from dimensionless monomials.
- Do not call a synthetic challenge a hardware validation.
- Do not call candidate-model discrimination full process tomography.
- Do not infer Lorentz symmetry from \(z=1\) alone.
- Do not infer OS reconstruction from a positive quadratic matrix alone.
- Do not call a post-hoc potential predictive.
- Do not transfer legacy physical roles onto the strict kernel without a bridge and role-transfer theorem.

---

# Final scientific judgment

Programs 21–30 materially improve the theory's mathematical and operational map. They construct two positive refinements, an exact local Dirac square root of the cycle Laplacian, an exact nonlocal square root of the strict generator, a free reflection-positive temporal extension, robust bounded-nuisance discriminators, explicit dilation nonidentifiability, a complete relational-frame theorem, a concrete adaptive nonintegrability certificate, an executable no-go library, and a blinded synthetic test with a measured resolution boundary.

They do not turn the strict operator into physics.

The inverse-action study sharpens the reason. Reconstructing a Gaussian precision from a positive covariance is mathematically standard. Reconstructing a unique law that selects an observed contact-subtracted kernel is not. The first problem is solved conditionally by inversion or by the covariance-level log-determinant functional. The second remains nonunique unless the potential, drive, locality, contact term, units, temporal category, and empirical interface are specified independently.

The deepest surviving interpretation is therefore:

> FIN currently supplies a distinguished finite spectral/Dirichlet object together with a rich family of mathematically admissible shadows. Programs 21–30 identify which operational experiments can separate those shadows and which inverse constructions are circular. The missing physical principle is not another spectral identity; it is an independently sourced operational and dimensional law that selects one completion and survives held-out experiment.

This conclusion is **[Proven]** as a nonselection statement over the explicitly constructed countermodels, and **[Moderate evidence]** as a judgment about the future theory.

---

# Reproducibility statement

Run:

```text
MPLCONFIGDIR=/tmp/mplconfig python3 fin_programs_21_30_experiments.py \
  --output FIN_Programs_21_30_Results.json \
  --figures-dir FIN_Programs_21_30_Figures
```

The archival release should contain:

- this Markdown source;
- the English PDF;
- the generated LaTeX source;
- *fin_programs_21_30_experiments.py*;
- *FIN_Programs_21_30_Results.json*;
- the figure directory;
- a SHA-256 manifest.

## Suggested citation before DOI assignment

Żuchowski, K. (2026). *FIN Beyond the Spectral Core: Executed Programs 21–30, Operational Falsification, and Inverse Reconstruction of an Action from the Strict Kernel* (FIN Research Supplement, Release 10.2; Version 1.0.0) [Preprint]. Zenodo. DOI to be assigned.
