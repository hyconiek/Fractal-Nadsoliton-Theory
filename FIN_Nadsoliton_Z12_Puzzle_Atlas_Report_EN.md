# Nadsoliton Puzzle Atlas

## Z12 Simulator Audit, Controlled Cross-Disciplinary Intuition, and a New Memory Object after Reduction

**Version:** Release 10.22 — English research edition  
**Author:** Krzysztof Żuchowski  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Date:** 28 July 2026  
**Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Confidence convention

- **[Proven]** — an algebraic or numerical result with an explicit test and
  tolerance.
- **[Strong evidence]** — stable in the declared test class, but not a general
  theorem.
- **[Moderate evidence]** — a meaningful mechanism with unresolved
  identifiability.
- **[Speculative]** — a research hypothesis, not a claim about nature.
- **[Refuted]** — a claim falsified in the declared model.

## 1. Executive verdict

The most productive “assembly of the puzzles” does not yield a new scalar
parameter or a hidden selector. It yields the operator-valued object

\[
\Sigma_E(z)=A_{EH}(zI+A_{HH})^{-1}A_{HE},
\]

the frequency-dependent **self-energy of the reduced process**, obtained after
partitioning the nodes into observed \(E\) and hidden \(H\) sectors.
**[Proven]**

This object unifies five previously separate pictures:

1. Green function/resolvent, through the exact block-resolvent formula;
2. Schur compression, as the static limit \(\Sigma_E(0)\);
3. wave and diffusion dynamics, through the same analytic resolvent family;
4. memory and environment, through an exact time-domain memory kernel;
5. the observer, as an explicitly supplied accessible/hidden partition rather
   than a metaphysical component of the equation.

It does not solve the three principal omissions of FIN: it does not select an
orientation, generate a physical unit, or provide an experimental preparation
and apparatus. **[Proven]**

The deepest intuition supplied by the atlas is therefore:

> Fractal reduction should not be modeled solely as a static transformation
> of one kernel into another. Exact elimination of an information level
> generates a process with memory. The static Schur complement is only its
> zero-frequency shadow.

This is a known reduction mechanism in mathematical physics, but its concrete
form for the strict FIN operator and its relation to wave–diffusion duality
provide a new testable direction inside the project. It is not evidence of new
physics.

## 2. Method: intuition with safeguards

“Intuitive” comparison can easily turn graphical resemblance into false
equivalence. The atlas therefore uses four layers:

1. **Shape layer:** identify similar plots, symmetries, flows, and diagrams.
2. **Mechanism layer:** determine whether both structures use the same
   algebraic operation.
3. **Prediction layer:** ask whether the analogy produces a new computable and
   falsifiable quantity.
4. **Assumption-debt layer:** list everything that must be added before the
   analogy becomes a physical model.

Each match receives one of five classifications:

- **exact equivalence** — an explicit operation-preserving map exists;
- **common mechanism** — the operation is shared but the objects or
  interpretations differ;
- **structural analogy** — similar geometry or spectrum without an
  equivalence theorem;
- **metaphor** — visually useful, with no evidential force;
- **anti-analogy** — the resemblance exposes the precise point at which FIN
  differs from the comparison theory.

The script `fin_nadsoliton_puzzle_atlas.py` scanned:

- 12 individual characteristics against 15 fields: **180 comparisons**;
- all unordered pairs of characteristics: \(\binom{12}{2}=66\);
- every pair against 15 fields: **990 comparisons**.

This yields **1170 explicit similarity indices**. Tag-based scoring is only a
candidate-search device. It is not evidence and never raises the confidence
status of a hypothesis.

![Heuristic matching matrix of individual puzzles against scientific fields.](FIN_Nadsoliton_Puzzle_Atlas_Figures/single_puzzle_domain_atlas.png)

## 3. What `Z12 sim.html` actually does

The code implements a finite thought experiment on \(Z_{12}\):

\[
W_{xy}=K(d(x,y)),\qquad
A=sI-W,\qquad
U_t=e^{-itA}.
\]

The initial state is a partially coherent preparation on two nodes:

\[
\rho_0(\phi,\eta)
=\frac12\Big(
|a\rangle\langle a|+|b\rangle\langle b|
+\eta e^{-i\phi}|a\rangle\langle b|
+\eta e^{i\phi}|b\rangle\langle a|
\Big).
\]

After the evolution \(\rho_t=U_t\rho_0U_t^\dagger\), the code evaluates

\[
p_x(t)=(\rho_t)_{xx},
\qquad
J(x\to y)=2W_{xy}\operatorname{Im}(\rho_{yx}),
\]

and two different aggregate current observables:

\[
C_1=\sum_xJ(x\to x+1),
\]

\[
C_\chi=\sum_x\sum_{d=1}^{5}d\,J(x\to x+d).
\]

\(C_1\) sees only the first hop harmonic. \(C_\chi\) is an oriented group
current of the full kernel, excluding the antipodal edge \(d=6\), because a
sign cannot be assigned to the shortest displacement there without an
additional convention.

The simulator’s intended distinction is sound: a **sign receiver** is not a
**sign source**. A radial operator can transport prepared chirality into a
measurable current, but it does not choose the sign of the phase.

## 4. Methodological errors and corrections

### 4.1. Inconsistent current sign

Two diagonal contributions used the opposite sign of the imaginary part from
the interference contributions, violating the single definition of
\(J(x\to y)\).

**Correction:** every contribution is now computed as

\[
2W_{xy}\operatorname{Im}(\rho_{yx}).
\]

**[Proven]**

### 4.2. False claim for \(\eta=0\)

Removing coherence eliminates cross terms but need not eliminate the local
currents of the two independently propagated packets. The audit gives

\[
\max_{x,y}|J_{xy}|_{\eta=0}=0.0995050916,
\qquad
C_\chi|_{\eta=0}\approx0.
\]

**Correction:** the interface now states that interference circulation
vanishes, not that all local currents vanish. **[Proven]**

### 4.3. Confusing one harmonic with the full current

For the preparation \(a=10,b=2\),

\[
C_1(\phi)\approx0
\]

for every tested phase, even though the full operator contains long-range
hops. At \(\phi=0.7\),

\[
C_\chi=+0.1211880893,\qquad
C_\chi(-0.7)=-0.1211880893,
\]

while

\[
C_1(0.7)=-2.36\times10^{-16}.
\]

**Correction:** the simulator displays \(C_\chi\) as the principal counter and
\(C_1\) as a blind-harmonic control. **[Proven]**

### 4.4. Noncanonical “circulation over all edges”

The original sum weighted distant edges by an implicitly chosen cyclic
displacement without declaring the branch choice or antipodal convention.

**Correction:** the observable is now explicitly defined as \(C_\chi\), with a
statement that it depends on a chosen orientation and lift of distance. It is
a chirality detector, not its source.

### 4.5. False “asymmetry” from changing separation

The pair \(a=10,b=3\) still has a reflection interchanging the two equal
sources. The test at \(\phi=0\) gives \(C_\chi\approx0\).

**Correction:** the switch is labeled “different separation,” not symmetry
breaking. **[Proven]**

### 4.6. False monotonicity in \(|\phi|\)

\(C_\chi\) is odd and periodic:

\[
C_\chi(-\phi)=-C_\chi(\phi),\qquad
C_\chi(0)=C_\chi(\pi)=0.
\]

It may grow linearly near zero, but not globally with \(|\phi|\).

**Correction:** the narration and plot now describe periodic behavior.
**[Proven]**

### 4.7. Spectrum–kernel mismatch

The spectrum of \(A\) was pasted as a fixed array, so a change of \(K\) could
leave stale eigenvalues.

**Correction:** the spectrum is now computed from the kernel by Fourier
transform, and the page verifies

\[
s=1.660307278766099,\qquad
\lambda_{\min}(A)\approx0,\qquad A\succeq0.
\]

### 4.8. Missing continuity-equation control

**Correction:** the page checks normalization and

\[
\dot p_x+\sum_yJ(x\to y)=0.
\]

The independent NumPy audit gives a maximum residual of
\(5.55\times10^{-17}\).

### 4.9. Misleading visual elements

The central arrow rotated even at zero current, while the gauge scale depended
on interaction history.

**Correction:** zero current is now static, and the scale is determined
reproducibly from the current \(C_\chi(\phi)\) curve.

### 4.10. Audit verdict

`audit_z12_sim.py` passes **15/15 tests**. **[Proven]**

This does not show that FIN is physically realized. It establishes only the
internal consistency of the declared finite model.

## 5. Twelve individual puzzles

| ID | FIN characteristic | Nearest known structures | Shared feature | Principal difference |
|---|---|---|---|---|
| C01 | strict oscillatory-damped kernel | spectral graph theory, graph signal processing, band-pass kernels | radial filter and Fourier spectrum | no self-generated continuum or physical interpretation |
| C02 | legacy intermediate bridge profile | screened Green functions, oscillatory propagators | oscillation, damping, amplitude | legacy roles cannot be transferred to strict without a theorem |
| C03 | positive Laplacian \(A\) | resistor networks, Markov processes, Dirichlet forms | positivity, constant preservation, gradient energy | time and units are dimensionless |
| C04 | duality \(e^{-itA}\) and \(e^{-tA}\) | quantum walks, heat kernels, Wick-like continuation | common functional calculus | common generator does not imply operationally identical processes |
| C05 | Green function, resolvent, action | Gaussian fields, inverse field theory, Yukawa-like propagation | operator inverse and quadratic action | interactions, measure, and action scale are additional data |
| C06 | orientation cocycle and chiral receivers | cohomology, holonomy, topological phases, asymmetry resources | sign carrier and reflection transformation | the pair \(\pm\) exists, but no sign is selected |
| C07 | fractal Schur compression | Kron reduction, multigrid, RG, tensor networks | elimination of degrees of freedom | the strict family is not closed under the audited Schur map |
| C08 | adaptive operator and memory | neural operators, reservoir computing, adaptive control | feedback and memory state | learning law and target are not derived from strict |
| C09 | information and geometry | diffusion maps, Fisher geometry, optimal transport | metrics from distributions and semigroups | internal metric is not an SI metre |
| C10 | state–clock–instrument–record | process tensors, causal models, system identification | interventions and multitime records | concrete apparatus remains external |
| C11 | missing dimensional scale | dimensional analysis, renormalization, calibration | scale orbit and need for a standard | no mass, energy, or SI time without an anchor |
| C12 | missing selector | spontaneous symmetry breaking, reference frames, resource theory | multiple branches and asymmetry cost | bifurcation gives branches, not a canonical choice |

## 6. Multi-puzzle assemblies

### 6.1. Spectral core: C01 + C03 + C04

This is the assembly closest to a finite graph transport theory.
**[Proven]**

- \(W\) defines connectivity geometry.
- \(A=sI-W\) defines Dirichlet energy.
- \(e^{-itA}\) produces coherent propagation.
- \(e^{-tA}\) produces relaxation/diffusion.

The equivalence is exact at the level of functional calculus, but not at the
level of experimental records.

### 6.2. Propagator inverted into an action: C01 + C03 + C05

If \(G=(L+m^2I)^{-1}\), then the minimal quadratic action

\[
S[\phi]=\frac12\phi^\top(L+m^2I)\phi-J^\top\phi
\]

has stationary equation \((L+m^2I)\phi=J\) and response \(\phi=GJ\).
**[Proven]**

This correctly reconstructs an action from a propagator. It is not a unique
fundamental law, because interactions, measures, and sectors invisible in the
two-point Green function may be added.

### 6.3. Chirality as a resource: C06 + C12

The closest framework is the resource theory of asymmetry. Reflection-covariant
operations cannot create asymmetry from a reflection-symmetric state. The
prepared phase \(\phi\) is a resource; \(C_\chi\) is a receiver.
**[Proven]**

This explains why the simulator can read a sign but cannot discharge
`QW-2191`.

### 6.4. Information as geometry: C04 + C09

The heat semigroup defines diffusion distances and multiscale spectral
coordinates, closely resembling diffusion maps. **[Strong evidence]**

The anti-analogy is essential: \(P_t\) can generate a **dimensionless**
geometry, but not a metre, second, or energy without calibration.

### 6.5. Operational observer: C04 + C10

The same \(A\) does not determine the observed process. One also needs

\[
(\mathcal H,\rho,\mathcal E_t,\tau,\mathcal P,\mathcal I,\mathcal E,
\mathcal A,\mathcal R,\mathcal C,\mathcal D).
\]

State, intervention, environment, and recording determine whether the record
appears wave-like, diffusive, or history dependent. **[Proven]**

### 6.6. Static compression: C05 + C07

Schur/Kron reduction yields the exact boundary response at a fixed resolvent
parameter. This is an equivalence of block algebra, not proof of strict
self-similarity. **[Proven]**

The previously audited failure of the strict family to close under the Schur
map remains in force.

### 6.7. Dynamic compression: C04 + C05 + C07 + C10

This is the richest new assembly. When some nodes are invisible, the exact
observed process is not generated by one static Schur complement. The
frequency-dependent \(\Sigma_E(z)\) and a memory kernel appear instead.
**[Proven]**

### 6.8. Fractality plus memory: C07 + C08

If successive layers are eliminated, each contributes its own self-energy and
memory kernel, producing a hierarchy

\[
\Sigma^{(1)}(z),\Sigma^{(2)}(z),\ldots.
\]

This may resemble continued fractions, multiscale impedance, or hierarchical
memory models. **[Moderate evidence]**

No result proves that the hierarchy reproduces the strict kernel or exponent
\(9/5\).

### 6.9. Adaptation plus a process with memory: C08 + C10

An adaptive operator learning from process records resembles system
identification with memory and neural operators. **[Moderate evidence]**

The model may learn the apparatus rather than FIN. Permutation controls,
hold-out tests, and identifiability analysis are essential.

### 6.10. Two independent obstructions: C11 + C12

The missing scale and missing orientation are different torsors:

- scale carries an action of \(\mathbb R_{>0}\);
- orientation carries an action of \(\mathbb Z_2\) or a relevant automorphism
  group.

A positive scalar cannot select a sign, and a pseudoscalar cannot generate a
unit of length. **[Proven]**

This excludes the shortcut in which one “informational constant” is expected
to generate both dimensional units and a selector.

### 6.11. Legacy + Green + strict: C02 + C05 + C01

Legacy may be viewed as an intermediate propagator profile and strict as a
later enrichment. **[Moderate evidence]**

This does not license silent substitution of strict for legacy, transfer of
legacy physical roles, or treatment of amplitude/damping as a complete
completion map.

### 6.12. Spectral action: C03 + C05 + C09

The trace \(\operatorname{Tr}f(A/\Lambda)\) resembles the spectral-action
principle, but \(A\) alone is not a complete spectral triple and FIN does not
generate \(\Lambda\). **[Moderate evidence]**

This is a useful anti-analogy: the existing theory identifies the additional
objects required—an algebra, representation, Dirac-type operator,
grading/real structure, and scale.

## 7. New working theorem: dynamic Schur reduction generates memory

Partition the space into observed nodes \(E\) and hidden nodes \(H\):

\[
A=
\begin{pmatrix}
A_{EE}&A_{EH}\\
A_{HE}&A_{HH}
\end{pmatrix}.
\]

### Theorem 1 — exact resolvent block

For \(z>0\),

\[
\left[(zI+A)^{-1}\right]_{EE}
=
\left[
zI+A_{EE}
-A_{EH}(zI+A_{HH})^{-1}A_{HE}
\right]^{-1}.
\]

Therefore,

\[
\Sigma_E(z)=A_{EH}(zI+A_{HH})^{-1}A_{HE}.
\]

**Proof.** This is the block-matrix inverse formula using the Schur complement.
**[Proven]**

### Theorem 2 — no single exact static reduced generator

For \(z>0\),

\[
\Sigma_E'(z)
=
-A_{EH}(zI+A_{HH})^{-2}A_{HE}\preceq0.
\]

If coupling to the hidden sector is nonzero in an active direction, then
\(\Sigma_E'(z)\neq0\). Consequently, no constant matrix \(B\) can satisfy

\[
\left[(zI+A)^{-1}\right]_{EE}=(zI+B)^{-1}
\]

for every \(z>0\).

**Proof.** The right-hand side would require

\[
B=A_{EE}-\Sigma_E(z)
\]

to be independent of \(z\), contradicting \(\Sigma_E'(z)\neq0\).
**[Proven]**

### Theorem 3 — exact composition defect after projection

Let \(P\) select sector \(E\), \(Q=I-P^\ast P\), and let \(T(t)\) be either
\(e^{-tA}\) or \(e^{-itA}\). Then

\[
PT(t+s)P^\ast
-PT(t)P^\ast PT(s)P^\ast
=
PT(t)QT(s)P^\ast.
\]

The right-hand side is the amplitude or mass of trajectories that left the
observed sector and returned. **[Proven]**

### Time-domain form

For heat dynamics,

\[
\dot x=-A_{EE}x-A_{EH}h,\qquad
\dot h=-A_{HE}x-A_{HH}h.
\]

Exact elimination of \(h\) gives

\[
\dot x(t)
=-A_{EE}x(t)
+\int_0^t M_E(t-s)x(s)\,ds
-A_{EH}e^{-tA_{HH}}h(0),
\]

\[
M_E(t)=A_{EH}e^{-tA_{HH}}A_{HE}.
\]

This is a memory equation. The “environment” has not been inserted by hand; it
arises mathematically from the declared split between visible and hidden
degrees of freedom. **[Proven]**

### Result for strict \(Z_{12}\)

Choose

\[
E=\{0,2,4,6,8,10\},\qquad
H=\{1,3,5,7,9,11\}.
\]

| Test | Result |
|---|---:|
| static-Schur row-sum residual | \(2.78\times10^{-16}\) |
| maximum resolvent-identity residual | \(3.83\times10^{-15}\) |
| \(\|\Sigma'(0.05)\|_2\) | \(0.919783\) |
| \(\|\Sigma'(1)\|_2\) | \(0.290955\) |
| heat composition defect at \(t=0.5\) | \(0.119015\) |
| wave composition defect at \(t=0.5\) | \(0.305427\) |
| static-Schur heat error at \(t=1\) | \(0.451941\) |
| static-Schur wave error at \(t=1\) | \(0.921176\) |
| exact hidden-excursion identity residual | \(<6.3\times10^{-16}\) |

![Composition defect and the error of replacing dynamic reduction by a static Schur generator.](FIN_Nadsoliton_Puzzle_Atlas_Figures/dynamic_schur_memory.png)

### Significance

- **[Proven]** Fractal elimination of degrees of freedom generates memory even
  when the full dynamics is a homogeneous semigroup or group.
- **[Proven]** Wave and diffusion dynamics share the same block-resolvent
  object but produce different time records.
- **[Moderate evidence]** \(\Sigma_E(z)\) is a strong candidate for the missing
  mathematical interface between compression, environment, and apparatus.
- **[Refuted]** Static \(A_{\rm Schur}\) is not an exact generator of the entire
  reduced dynamics for the audited partition.
- **[Refuted]** Memory after reduction does not establish a fundamental arrow
  of time or a conscious observer.

## 8. “Shadow shapes”: new intuitions

### 8.1. The missing object may be a noncommutativity defect

Many unresolved components share the form

\[
\text{reduce first, then evolve}
\neq
\text{evolve first, then reduce}.
\]

The defect of this commutation is an observable shadow of hidden degrees of
freedom. **[Strong evidence]**

### 8.2. The observer is a section of a process, not an added substance

Mathematically, an observer may mean:

- a selected subalgebra of accessible observables;
- a projection \(P\);
- a set of permitted instruments;
- a rule for recording multitime data.

This does not introduce an “informational layer beneath the nadsoliton.” It is
an operational section of nadsoliton states. **[Moderate evidence]**

### 8.3. \(C_\chi\) and the chiral bispectrum receive the same resource type

Both are odd under inversion and reverse sign with the preparation. They may
belong to a common class of monotone asymmetry receivers.
**[Moderate evidence]** Neither is thereby a sign selector.

### 8.4. Different observables see different harmonics

The blindness of \(C_1\) when \(C_\chi\neq0\) is a general tomographic lesson.
One current or one interference image does not identify the full operator.
**[Proven]**

### 8.5. Static and dynamic self-similarity are distinct hypotheses

Failure of strict to close under a static Schur map does not rule out a simpler
structure for the family \(\Sigma^{(n)}(z)\). This is a new, well-typed
hypothesis. **[Speculative]**

### 8.6. Scale can be relative, but an experiment requires an anchor

Spectral ratios are identifiable without an absolute clock; absolute energies
and times are not. **[Proven]** A projective fingerprint can be tested first,
but a physical theory still requires calibration.

### 8.7. Information geometry may be real without being spacetime

The heat kernel gives a diffusion metric and scale hierarchy but does not
establish Lorentz signature or spacetime geometry. **[Proven]**

### 8.8. The most defensible architecture still has two additional packages

- \(W_0\): the strict dimensionless information-operator core;
- an operational/conversion package: preparation, scale, instrument,
  environment, and record;
- a sector package: an explicit orientation/branch resource.

Dynamic self-energy may enrich the first interface, but it does not remove the
other two packages.

## 9. Falsification of analogies

| Tempting analogy | Destructive test | Verdict |
|---|---|---|
| strict = Yukawa propagator | test whether the inverse and local operator are derived rather than fitted | conditional reconstruction only |
| duality = observer paradox | same generator but different instruments and records | no mathematical paradox |
| Schur = RG fixed point | test closure of the strict family | [Refuted] in the audited realization |
| chiral current = selector | phase inversion gives an equally valid opposite sign | [Refuted] as a source |
| bifurcation = unique choice | symmetric law produces a pair of branches | [Refuted] without bias |
| entropy = energy | no temperature, Hamiltonian, or apparatus reset | [Refuted] without additional axioms |
| diffusion distance = physical length | generator rescaling changes calibration | [Refuted] as an SI metre |
| process tensor = ontology | depends on accessible interventions and system/environment split | operational description only |
| \(\Sigma(z)\) = new fundamental field | changing the projection changes \(\Sigma\) | [Refuted] as projection independent |
| graphical resemblance = theory equivalence | no operation- and observable-preserving maps | [Refuted] methodologically |

## 10. Twelve subsequent research programs

### Program 243 — Dynamic Schur theorem for arbitrary partitions

Prove conditions for positivity and complete monotonicity of \(\Sigma_E(z)\),
minimality of the memory realization, and a necessary-and-sufficient criterion
for exactly Markovian reduction.

**Probability of a valuable result:** 0.90.  
**Success criterion:** an if-and-only-if theorem plus tests of every
nonisomorphic \(Z_{12}\) partition.

### Program 244 — Common analytic self-energy of wave and diffusion

Build one resolvent package and derive its limits

\[
z>0\quad\text{(heat)},\qquad z=\epsilon-i\omega\quad\text{(wave)}.
\]

**Probability:** 0.82.  
**Risk:** confusing analytic continuation with physical Wick rotation without
axioms.

### Program 245 — Identifiability of the observed/hidden partition

Determine whether distinct \(E/H\) splits can generate the same process on
\(E\), and construct equivalence classes of minimal realizations.

**Probability:** 0.80.

### Program 246 — \(C_\chi\) as a twist-response theorem

Define a flux-twisted family \(A(\theta)\) and test

\[
C_\chi\stackrel{?}{=}
\left.\frac{d}{d\theta}
\operatorname{Tr}[\rho A(\theta)]
\right|_{\theta=0}
\]

with correct treatment of \(d=6\), gauge covariance, and orientation reversal.

**Probability:** 0.86.

### Program 247 — Current-harmonic tomography

Measure

\[
C_d=\sum_xJ(x\to x+d),\qquad d=1,\ldots,5
\]

instead of only \(C_1\), and determine the minimum preparations and
observables required to identify the chiral part of a state.

**Probability:** 0.84.

### Program 248 — Asymmetry-resource budget

Compare \(C_\chi\), the chiral bispectrum, and the earlier
\(\Lambda(\rho,A)\) as receivers of one resource. Prove bounds of the form

\[
|\langle C\rangle_\rho|\leq M_R(\rho)\|C\|_\infty.
\]

**Probability:** 0.78.  
**Boundary:** classification of receivers, not a new selector.

### Program 249 — Fractal cascade of memory

Iterate layer elimination and study

\[
\Sigma^{(n+1)}(z)=\mathcal F(\Sigma^{(n)}(z)).
\]

Seek a fixed point in operator-function space rather than in the existing
static-kernel family.

**Probability:** 0.62.  
**Reward:** high; this would be the proper dynamic version of fractal
compression.

### Program 250 — Process-memory test with a complete record

Extend P240–P242 by an intervention at an intermediate time. Distinguish:

- a single semigroup on 12 states;
- a static six-state model;
- the exact six-state process with memory.

**Probability:** 0.75, laboratory dependent.

### Program 251 — False positives of the spectral fingerprint

Generate broad Laplacian classes and determine how many non-FIN operators pass
the projective-spectrum, semigroup, harmonic-current, and prescribed-reduction
memory tests.

**Probability:** 0.88.

### Program 252 — Stability of heat/Fisher geometry under reduction

Compare diffusion distances and the Fisher metric before and after Schur
reduction and dynamic reduction.

**Probability:** 0.73.  
**Boundary:** the result remains dimensionless.

### Program 253 — Can dynamic self-energy generate strict damping?

Treat this as one new bridge atom, not a repetition of a generic
legacy-to-strict audit. Test whether elimination of a declared self-similar
layer can generate

\[
(1+\beta d^{9/5})^{-1}
\]

without encoding \(9/5\) in the input.

**Probability:** 0.38.  
**Stop rule:** one explicit ansatz; do not enlarge the family after a no-go.

### Program 254 — Blind hidden-excursion test

Preregister

\[
T_E(2\tau)-T_E(\tau)^2
=PT(\tau)QT(\tau)P^\ast
\]

and measure both sides from independent records.

**Probability:** 0.70.  
**Value:** directly distinguishes static reduction from a process with memory.

## 11. Research priority

| Priority | Program | Reason |
|---:|---|---|
| 1 | P243 | highest theorem confidence and immediate clarification of reduction |
| 2 | P246 | completes simulator methodology and current definition |
| 3 | P247 | prevents blindness to all but one harmonic |
| 4 | P244 | formally joins wave and diffusion after reduction |
| 5 | P251 | measures specificity of the complete FIN fingerprint |
| 6 | P245 | determines what can actually be inferred about an “environment” |
| 7 | P250 | moves the result into an operational multitime process |
| 8 | P248 | orders chiral receivers without pretending to find a source |
| 9 | P252 | tests information geometry |
| 10 | P254 | supplies a falsifiable memory experiment |
| 11 | P249 | high reward, greater risk |
| 12 | P253 | most difficult and most vulnerable to target coding |

## 12. What this report does not claim

The report does not export:

- a strict selector or discharge of `QW-2191`;
- a canonical scale of length, time, mass, energy, or action;
- closure of the legacy-to-strict bridge;
- transfer of legacy physical roles;
- \(L_{\rm total}\), the Standard Model, GR, or a ToE;
- evidence that any laboratory realizes FIN;
- evidence that projection-generated memory is fundamental memory in nature.

## 13. Shortest interpretation surviving falsification

The strict FIN core is a finite positive spectral generator whose different
functions produce coherent propagation, diffusion, a Green function, and a
quadratic action. Chirality belongs to resources of states or operations, not
to the radial generator. When some degrees of freedom are hidden—as naturally
occurs both in fractal compression and for an observer—the exact shadow of the
full theory is not another static kernel but a frequency-dependent self-energy
and a process with memory.

This is currently the richest mathematical picture obtained by assembling the
puzzles. It is more precise and more testable than the claim that one plot
merely “looks like” a familiar phenomenon.

## 14. Comparative literature

1. Dörfler, F.; Bullo, F., *Kron Reduction of Graphs with Applications to
   Electrical Networks*, <https://arxiv.org/abs/1102.2950>.
2. Lin, Y. T.; Tian, Y.; Anghel, M.; Livescu, D., *Data-driven learning for
   the Mori–Zwanzig formalism*, <https://arxiv.org/abs/2101.05873>.
3. Pollock, F. A. et al., *Non-Markovian quantum processes: complete framework
   and efficient characterisation*, <https://arxiv.org/abs/1512.00589>.
4. Marvian, I.; Spekkens, R. W., *Modes of asymmetry*,
   <https://arxiv.org/abs/1312.0680>.
5. Shuman, D. I. et al., *The Emerging Field of Signal Processing on Graphs*,
   <https://arxiv.org/abs/1211.0053>.
6. Coifman, R. R.; Lafon, S., *Diffusion maps*,
   <https://doi.org/10.1016/j.acha.2006.04.006>.
7. Chamseddine, A. H.; Connes, A., *The Spectral Action Principle*,
   <https://arxiv.org/abs/hep-th/9606001>.

## 15. Reproducible artifacts

- `Z12 sim.html` — corrected interactive simulator.
- `audit_z12_sim.py` — independent mathematical audit.
- `Z12_sim_methodological_audit.json` — results of 15 tests.
- `fin_nadsoliton_puzzle_atlas.py` — atlas and dynamic-Schur study.
- `FIN_Nadsoliton_Puzzle_Atlas_All_Singles.csv` — 180 comparisons.
- `FIN_Nadsoliton_Puzzle_Atlas_All_Pairs.csv` — 990 comparisons.
- `FIN_Nadsoliton_Puzzle_Atlas_Results.json` — values and candidate ranking.
- `FIN_Nadsoliton_Puzzle_Atlas_Figures/` — figures.
