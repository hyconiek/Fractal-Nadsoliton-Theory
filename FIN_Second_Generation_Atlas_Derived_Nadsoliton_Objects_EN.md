# Second-Generation Atlas

## Derived Nadsoliton Objects

**Subtitle:** Systematic search over combinations C01...C12, construction of
objects O01...O15, and third-generation synthesis G01...G10  
**Version:** Release 10.23 — English edition  
**Author:** Krzysztof Żuchowski  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Date:** 28 July 2026  
**Repository:** <https://github.com/hyconiek/Fractal-Nadsoliton-Theory>

## Confidence convention

- **[Proven]** — an analytic result or an explicit numerical identity within a
  stated tolerance.
- **[Strong evidence]** — a stable construction whose general theorem or
  interpretation remains open.
- **[Moderate evidence]** — a typed mechanism whose identifiability or
  uniqueness remains open.
- **[Speculative]** — a constructive hypothesis intended for falsification.
- **[Refuted]** — a connection rejected in the stated class, usually because
  no common typed operation exists.

## 1. Executive summary

The second phase of the atlas no longer asks only what the puzzles C01...C12
resemble. It asks:

\[
\text{which new objects can be constructed by composing their operations?}
\]

A finite, reproducible search was performed:

| Class | Count |
|---|---:|
| individual puzzles | 12 |
| all pairs | 66 |
| all triples | 220 |
| selected couplings of 4–6 puzzles | 7 |
| combinations scanned in total | 305 |
| rows yielding a typed object or mechanism | 109 |
| purely visual combinations rejected | 184 |
| second-generation O-objects | 15 |
| third-generation G-objects | 10 |

![Coverage of the combinatorial search and output of the construction grammar.](FIN_Second_Generation_Atlas_Figures/combination_search_coverage.png)

The strongest derived object is not a single matrix but a contextual family:

\[
\mathfrak S_A:
E\longmapsto
B_E(z)
=\operatorname{Schur}_{V\setminus E}(zI+A),
\qquad z>0.
\]

For every observed context \(E\), the hidden sector generates the self-energy

\[
\Sigma_E(z)
=A_{EH}(zI+A_{HH})^{-1}A_{HE}.
\]

This family has three exact properties:

1. nested reductions compose exactly by associativity of the Schur complement;
2. \(\Sigma_E(z)\) is a positive operator-valued Stieltjes function;
3. its inverse Laplace transform is the memory kernel of the reduced process.

We call the resulting G02 object the **Stieltjes–Schur Context Functor**.

It is **[Proven]** as finite operator algebra and **[Speculative]** as a
foundation for physical ontology.

The second particularly promising object is the chiral memory susceptibility

\[
\Xi_E(z)
=
\left.
\partial_\theta\Sigma_E(z,\theta)
\right|_{\theta=0}.
\]

It is nonzero and odd under reflection, but it is a **receiver** of the
inserted twist \(\theta\), not a source of orientation. **[Strong evidence]**

The third result clarifies “negative information coupling.” Loss of accessible
information can be represented without modifying the kernel by the contraction
ledger

\[
\mathcal L_t(p,q)
=D(p\Vert q)-D(P_tp\Vert P_tq)\geq0.
\]

This is loss of distinguishability in the accessible channel. It is not
automatically destruction of information or thermodynamic energy. **[Proven]**

## 2. Generation rather than analogy matching

Every combination passes through five gates:

1. **Common operation:** can the puzzle operations be composed with compatible
   input and output types?
2. **New object:** is the output more than a list of constituents?
3. **New space:** does the object live in a recognizable space of operators,
   functionals, measures, or functors?
4. **Observable and dynamics:** does the construction provide a testable
   number, record, or evolution?
5. **Falsification gate:** what outcome would refute the construction or its
   stronger interpretation?

Stop rule:

> A combination without a typed operation does not generate an object.
> Verbal or visual resemblance alone is marked [Refuted].

The complete result for all 305 combinations is stored in
`FIN_Second_Generation_Combination_Search.csv`. Each row records the common
operation, candidate object, object space, new observable, new dynamics,
status, and whether the result is a new coupling or merely inherits a
lower-order object.

## 3. Base inventory C01...C12

| ID | Base object | Dominant operation | Critical limitation |
|---|---|---|---|
| C01 | strict kernel | convolution/spectralization on \(Z_{12}\) | no automatic continuum |
| C02 | legacy intermediate bridge | propagator-like damped profile | no complete completion map or role transfer |
| C03 | \(A=sI-W\succeq0\) | Dirichlet form and spectral calculus | no energy unit |
| C04 | \(e^{-itA},e^{-tA}\) | two rays of functional calculus | different operational semantics |
| C05 | Green/resolvent/action | operator inversion and variation | propagator does not determine the full interaction |
| C06 | cocycle/chiral receivers | sign under inversion | no non-premise sign |
| C07 | fractal Schur compression | elimination of degrees of freedom | strict is not a static fixed point |
| C08 | adaptation and memory | operator–record feedback | learning law is additional input |
| C09 | information/geometry | entropy, Fisher metric, diffusion | no physical unit |
| C10 | operational process | preparation–instrument–environment–record | concrete laboratory remains external |
| C11 | scale obstruction | quotient by \(\mathbb R_{>0}\) | no canonical section |
| C12 | selector obstruction | orientation quotient/torsor | `QW-2191` remains open |

## 4. Main table: second-generation objects

| FIN combination | Derived object | Mathematical form | Analogous fields | New intuition | Status |
|---|---|---|---|---|---|
| C03+C05+C07 | O01 Dynamic Spectral Compression Kernel | \(\Sigma_E(z)=A_{EH}(zI+A_{HH})^{-1}A_{HE}\) | Feshbach, Mori–Zwanzig, Kron reduction | dynamic reduction generates memory | [Proven] |
| C03+C05+C07 | O02 Frequency-Dependent Effective Action | \(S_z[x]=\frac12\langle x,(zI+A_{EE}-\Sigma_E(z))x\rangle-\langle j,x\rangle\) | effective action, influence functional | Green, compression, and variation become one object | [Proven] |
| C05+C07+C10 | O03 Resolvent Context Presheaf | \(E\mapsto\operatorname{Schur}_{V\setminus E}(zI+A)\) | category theory, open systems, network reduction | each observer/context has a compatible effective operator | [Proven] |
| C03+C06+C12 | O04 Chiral Spectral Response Bundle | \((M_R,\Lambda,C_\chi,\operatorname{Im}B)\) | asymmetry resources, harmonic analysis | distinct sign receivers detect one resource type | [Strong evidence] |
| C06+C07+C10 | O05 Chiral Memory Susceptibility | \(\Xi_E(z)=\partial_\theta\Sigma_E(z,\theta)|_0\) | Kubo response, flux susceptibility, open systems | hidden degrees of freedom possess chiral memory response | [Strong evidence] |
| C03+C04+C09 | O06 Heat Information Metric Tower | \(D_t^2(x,y)=\sum_z|P_t(x,z)-P_t(y,z)|^2/\pi_z\) | diffusion maps, information geometry | informational distance is multiscale | [Proven] |
| C10+C11 | O07 Operational Calibration Torsor | \((A,\tau)\sim(cA,\tau/c)\) | metrology, gauge fixing, identification | experiment selects scale; algebra supplies an orbit | [Proven] |
| C01+C03+C11 | O08 Projective Spectral Fingerprint | \((\lambda_1/\lambda_{\max},\ldots)\) | inverse spectral theory, control | part of the prediction is testable before absolute calibration | [Proven] |
| C04+C09+C10 | O09 Information Contraction Ledger | \(\mathcal L_t=D(p\Vert q)-D(P_tp\Vert P_tq)\) | data processing, statistics, sensory bottlenecks | “loss” is loss of accessible distinguishability | [Proven] |
| C04+C07+C10 | O10 Dual-Ray Reduced Process Tensor | \(\mathcal T_r[\mathcal M_1,\ldots,\mathcal M_r]\) for wave/heat projections | process tensors, cybernetics, causal models | memory becomes visible through multitime interventions | [Moderate evidence] |
| C08+C09+C10 | O11 Adaptive Memory Geometry Flow | \((\dot A,\dot q)=(-\nabla_A\mathcal L,-\operatorname{grad}_F\mathcal L)\) | neural operators, adaptive networks, plasticity | learning may occur in memory-geometry space | [Speculative] |
| C05+C07+C08 | O12 Multiscale Memory Functor | \(\Sigma^{(0)}(z)\to\Sigma^{(1)}(z)\to\cdots\) | continued fractions, hierarchical reservoirs, RG | seek a fractal fixed point in memory-function space | [Strong evidence] |
| C01+C02+C07 | O13 Legacy-Strict Completion-Defect Tower | \(\Delta_n=K_{\rm strict}^{(n)}-\mathcal C_n[K_{\rm legacy}^{(n)}]\) | error propagation, RG defects, obstruction theory | study bridge failure dynamically instead of asserting similarity | [Moderate evidence] |
| C06+C09+C12 | O14 Spectral Asymmetry Transport Cost | \(\inf_{\Phi\ {\rm covariant}}\operatorname{Cost}(\Phi)\) | optimal transport, resource conversion | orientation has operational cost, but cost does not select a sign | [Moderate evidence] |
| C10+C11+C12 | O15 Two-Torsor Physicalization Bundle | \(\mathcal T_{\rm scale}\times\mathcal T_{\rm orient}\) | principal bundles, frames, calibration | scale and orientation require separate sections | [Strong evidence] |

## 5. Explanatory power and falsification

| Derived object | What it explains | What it does NOT explain | Falsifying test |
|---|---|---|---|
| O01 DSCK | exact memory after elimination | physical environment and scale | constant \(\Sigma_E(z)\) despite nonzero coupling |
| O02 Effective Action | exact contextual source response | nonlinear interactions and quantization | stationary solution disagrees with block Green function |
| O03 Context Presheaf | compatibility of nested reductions | unique observer context | nonzero direct-versus-nested residual |
| O04 Chiral Bundle | common covariance of sign receivers | source of sign | symmetric state with a nonzero odd receiver |
| O05 Chiral Memory | memory response to phase twist | selector for \(\theta\) | loss of oddness under reflection |
| O06 Metric Tower | geometry of distinguishability | SI length or Lorentz signature | distance grows under a contractive semigroup |
| O07 Calibration Torsor | non-identifiability of absolute scale | laboratory reference source | distinguish \((A,\tau)\) from \((cA,\tau/c)\) using the same \(P_\tau\) |
| O08 Fingerprint | scale-free spectral predictions | absolute energy | fingerprint changes under positive rescaling |
| O09 Information Ledger | accessible loss of distinguishability | thermodynamic energy or ontological destruction | violation of data processing |
| O10 Process Tensor | operational memory and interventions | concrete microscopic dilation | a memoryless model reproduces all interventions |
| O11 Adaptive Flow | possible operator-learning mechanism | strict source law or biology | no hold-out advantage over a static model |
| O12 Memory Functor | hierarchical dynamic compression | a strict fixed point | no stabilization of any functional class |
| O13 Defect Tower | where the bridge fails | role transfer and completion theorem | residual converges to zero without target coding |
| O14 Transport Cost | amount of asymmetry resource required | canonical choice of that resource | free operation increases the monotone |
| O15 Two-Torsor Bundle | independence of scale and orientation problems | source of either section | one internal datum canonically closes both torsors |

## 6. Ranking O01...O15

Scores use novelty, mathematical rigor, FIN relevance, and falsifiability from
0 to 5. Overinterpretation risk also ranges from 0 to 5, where 5 is highest.

![Score matrix for second-generation objects.](FIN_Second_Generation_Atlas_Figures/derived_object_score_matrix.png)

| ID | Novelty | Rigor | FIN | Falsifiability | Risk |
|---|---:|---:|---:|---:|---:|
| O01 | 4 | 5 | 5 | 5 | 2 |
| O02 | 3 | 5 | 5 | 5 | 2 |
| O03 | 4 | 5 | 4 | 5 | 2 |
| O04 | 3 | 4 | 5 | 5 | 3 |
| O05 | 4 | 4 | 4 | 5 | 3 |
| O06 | 2 | 5 | 4 | 4 | 3 |
| O07 | 2 | 5 | 5 | 5 | 2 |
| O08 | 2 | 5 | 5 | 5 | 2 |
| O09 | 2 | 5 | 4 | 5 | 4 |
| O10 | 3 | 3 | 5 | 5 | 3 |
| O11 | 4 | 2 | 3 | 3 | 5 |
| O12 | 4 | 4 | 4 | 4 | 3 |
| O13 | 3 | 4 | 5 | 5 | 2 |
| O14 | 3 | 3 | 4 | 4 | 4 |
| O15 | 2 | 4 | 5 | 5 | 2 |

O01 has the best balance. O01+O03+O12 yields the richest higher-order
structure. O11 has the greatest overinterpretation risk because analogies with
learning and biology are not FIN source laws.

## 7. Theorems generated in the second phase

### 7.1. O01 belongs to the operator Stieltjes class

For \(z>0\) and \(A_{HH}\succeq0\),

\[
\Sigma_E(z)=A_{EH}(zI+A_{HH})^{-1}A_{HE}.
\]

For every \(n\geq0\),

\[
(-1)^n\Sigma_E^{(n)}(z)
=n!A_{EH}(zI+A_{HH})^{-(n+1)}A_{HE}\succeq0.
\]

Thus \(\Sigma_E\) is a completely monotone positive operator-valued function.
**[Proven]**

Its poles and spectral measure encode hidden modes; its inverse Laplace
transform is a decaying memory kernel. Consequently, the admissible space for
dynamic compression is far narrower than the set of all matrix-valued
functions. For strict \(Z_{12}\), orders \(0,\ldots,4\) have minimum
eigenvalues above \(-6.4\times10^{-17}\), consistent with positivity up to
floating-point error.

### 7.2. O03 composes functorially

Let \(G\subset F\subset E\subset V\). For \(B(z)=zI+A\),

\[
\operatorname{Schur}_{E\setminus G}
\left(\operatorname{Schur}_{V\setminus E}B\right)
=\operatorname{Schur}_{V\setminus G}B.
\]

This is associativity of Gaussian elimination/the Schur complement.
**[Proven]** The direct-versus-nested residual in selected strict \(Z_{12}\)
contexts is \(2.47\times10^{-16}\).

“Functor” is used in the restricted sense that the context poset and reduction
morphisms form a compatible operator diagram. No sheaf axioms or quantum
contextuality theorem are claimed.

### 7.3. O02 reconstructs the exact contextual Green function

Set

\[
B_E(z)=zI+A_{EE}-\Sigma_E(z).
\]

The functional

\[
S_z[x;j]=\frac12\langle x,B_E(z)x\rangle-\langle j,x\rangle
\]

has stationary solution

\[
x_\star=B_E(z)^{-1}j
=\left[(zI+A)^{-1}\right]_{EE}j.
\]

The audited stationarity residual is \(4.90\times10^{-16}\). **[Proven]**
This exactly reconstructs an effective action from the resolvent. It is not a
complete physical action because \(z\), the action unit, measure, and
interactions remain additional data.

### 7.4. O05: memory has a chiral response direction

Introduce a Hermitian twisted family \(A(\theta)\), in which a hop by \(d\)
acquires phase \(e^{id\theta}\). Then

\[
\Xi_E(z)=
\left.\partial_\theta
\left[
A_{EH}(\theta)(zI+A_{HH}(\theta))^{-1}A_{HE}(\theta)
\right]\right|_{\theta=0}.
\]

For the even/odd split of strict \(Z_{12}\),

\[
\|\Xi_E(0.2)\|_2=1.2175076427.
\]

For reflection \(R\),

\[
R_E\Xi_E(z)R_E=-\Xi_E(z),
\]

with residual \(2.17\times10^{-16}\). **[Proven]** for the stated twisted
family and split.

Likewise,

\[
C_\chi=
\left.\frac{d}{d\theta}\operatorname{Tr}[\rho A(\theta)]
\right|_{\theta=0}
\]

has residual \(3.23\times10^{-13}\), promoting \(C_\chi\) from an ad hoc
indicator to an explicit flux-twist response. The twist and its sign are still
inserted in the test family; O05 does not discharge `QW-2191`.

### 7.5. O09: negative informational coupling as channel contraction

For a Markov channel \(P_t\) and positive distributions \(p,q\),

\[
D(P_tp\Vert P_tq)\leq D(p\Vert q).
\]

Define

\[
\mathcal L_t(p,q)=D(p\Vert q)-D(P_tp\Vert P_tq)\geq0.
\]

Across 1000 deterministically reproducible pairs:

| Quantity | Result |
|---|---:|
| minimum \(\mathcal L_t\) | 0.0485457 |
| median | 0.718427 |
| maximum | 3.15526 |
| violations below \(-10^{-12}\) | 0 |

This is **[Proven]** as an instance of the data-processing theorem. It is
methodologically preferable to inserting “negative information” into a scalar
kernel because it identifies the two hypotheses, the contracting channel, and
the difference between inaccessible and destroyed information. It does not
pretend to be energy without temperature and scale.

### 7.6. O07+O08: exact calibration orbit

\[
e^{-\tau A}=e^{-(\tau/c)(cA)}.
\]

The residual at \(c=7.25\) is \(3.43\times10^{-17}\), and the spectral
fingerprint remains invariant with residual \(6.66\times10^{-16}\).
**[Proven]**

An experiment can therefore test a scale-free fingerprint first, but physical
time requires an external calibration section.

### 7.7. O06: heat geometry contracts

The maximum squared diffusion distance on \(Z_{12}\) is:

| \(t\) | \(\max_{x,y}D_t^2(x,y)\) |
|---:|---:|
| 0.10 | 17.3358 |
| 0.25 | 11.0222 |
| 0.50 | 5.69186 |
| 1.00 | 2.00931 |

**[Proven]** for the tested family. The information geometry exists, but it is
dimensionless and semigroup-time dependent.

## 8. Second-order analogies

| FIN object | First analogy | Derived phenomenon | Returning FIN intuition | Boundary |
|---|---|---|---|---|
| O01 | Feshbach/Mori–Zwanzig | self-energy, generalized Langevin equation | hidden nodes generate memory | split is supplied |
| O02 | effective action | influence functional | contextual resolvent has a variational source | no complete QFT |
| O03 | context categories | compatible restrictions and diagrams | observer as context selection | no contextuality theorem |
| O04 | asymmetry modes | resource monotones, frames | receivers can be classified harmonically | no source resource |
| O05 | Kubo/flux response | susceptibility, pumping | memory has a chiral tangent | no topological quantization |
| O06 | diffusion maps | data geometry, coarse graining | heat generates internal geometry | no spacetime |
| O07 | gauge fixing/metrology | scale section | calibration belongs to operational theory | not a strict source |
| O08 | inverse spectra | system identification | test ratios before units | false positives possible |
| O09 | data processing | sensory bottleneck, sufficiency | information can enter hidden correlations | no Landauer energy |
| O10 | process tensors | non-Markovian interventions | multitime records distinguish memory | full instrument required |
| O11 | neural operators | maps between function spaces | learn \(\Sigma(z)\), not only a matrix | learning is not ontology |
| O12 | hierarchical reservoirs | fading memory, continued fractions | fixed point may live in Stieltjes space | theorem remains open |
| O13 | RG defect propagation | relevant/irrelevant directions | bridge defect may have a flow exponent | target coding forbidden |
| O14 | optimal transport | minimum conversion cost | orientation is an operational resource | cost does not choose sign |
| O15 | principal bundles | calibration and frames | two independent missing objects | sections remain external |

### Physics

The strongest mechanical analogies are Feshbach projection/self-energy,
Mori–Zwanzig memory kernels, open-system process tensors, linear response to a
twist or flux, and effective actions after field elimination. There is no
equivalence to a concrete field theory because FIN still lacks local
spacetime, fields, a measure, an action scale, and quantization rules.

### Mathematics

The natural new space is

\[
\operatorname{Stieltjes}_+(E)=
\left\{\Sigma:(0,\infty)\to\operatorname{End}(E)\mid
(-1)^n\Sigma^{(n)}(z)\succeq0\right\}.
\]

Nested Schur complements form a diagram over the context poset. Sheaf/global
section language remains a research analogy until covering data and gluing
axioms are defined.

### Computer science

O11 resembles neural operators because the learned object is the map

\[
\text{process record}\longmapsto\Sigma(z),
\]

not a single parameter vector. O12 resembles reservoir computing: the hidden
sector retains a trace of history and the current record is a readout. This
does not establish a biological or computational ontology for the nadsoliton.

### Biology

The defensible analogy concerns fading memory and adaptive weights: hidden
modes represent decay scales, plasticity would modify the memory spectrum, and
readout sees only a projection of the full state. **[Speculative]** No
biological data or cell/synapse-to-FIN-node map exists.

### Cybernetics

O03, O07, O08, and O10 form the scheme

\[
\text{system}\to\text{observer}\to\text{identification}
\to\text{calibration}\to\text{decision}.
\]

The missing observer need not be a new term in \(A\); it may be a family of
projections, instruments, and calibration sections.

## 9. Third generation: combining O-objects

![Construction dependencies between the O and G generations.](FIN_Second_Generation_Atlas_Figures/generation2_generation3_map.png)

| ID | Combination | Higher-order object | Definition/role | Status |
|---|---|---|---|---|
| G01 | O01+O02 | Dynamic Effective Action Bundle | Stieltjes self-energy plus exact source action | [Proven] |
| G02 | O01+O03+O12 | Stieltjes–Schur Context Functor | contexts, compatible reductions, completely monotone self-energies | [Proven] |
| G03 | O04+O05+O14 | Chiral Memory Resource Bundle | odd memory susceptibility plus asymmetry-resource budget | [Strong evidence] |
| G04 | O06+O12 | Memory–Diffusion Geometry Tower | geometry indexed by diffusion time and hiding depth | [Moderate evidence] |
| G05 | O07+O08 | Calibrated Fingerprint Bundle | scale-free test plus explicitly external scale section | [Proven] |
| G06 | O09+O10 | Operational Information-Balance Tensor | multitime information-contraction ledger | [Moderate evidence] |
| G07 | O01+O11 | Adaptive Self-Energy Flow | learned trajectory in Stieltjes-function space | [Speculative] |
| G08 | O12+O13 | Dynamic Completion-Defect Tower | tests whether bridge defect survives dynamic reduction | [Moderate evidence] |
| G09 | O03+O14+O15 | Equivariant Context-Resource Diagram | contexts, scale, orientation as separate sectional structures | [Speculative] |
| G10 | O05+O09 | Signed Information-Loss Susceptibility | derivative of information contraction with respect to chiral twist | [Speculative] |

## 10. Principal “shadow shapes”

### Shadow 1 — the missing object space is Stieltjes space

Rather than ask which next static kernel arises after compression, the
Green+Schur+dynamics combination shows that the proper object is the function
\(z\mapsto\Sigma(z)\). The self-similarity question changes from “is \(K\) a
fixed point?” to “is the family of memory functions closed or attracted to a
fixed class?” **[Strong evidence]** as a direction; no fixed point was found.

### Shadow 2 — observation and compression both select context

Compression selects retained degrees of freedom; an observer selects
accessible observables. Their common shadow is the context \(E\) and its
effective operator \(B_E(z)\). **[Moderate evidence]** Physical observation is
not thereby reduced to a Schur complement: an instrument and record remain
necessary.

### Shadow 3 — “negative information” is often flow into an invisible sector

O09 and O01 suggest

\[
\text{loss of accessible distinguishability}
\longleftrightarrow
\text{correlations/excursions in sector }H.
\]

Information may remain in a full reversible description while disappearing
from the selected record. **[Strong evidence]**

### Shadow 4 — chirality is a tangent, not a point

\[
\Xi_E\in T_{\Sigma_E}\operatorname{Stieltjes}_+(E).
\]

Orientation appears as a direction in memory-operator space. This new derived
characteristic does not choose whether that tangent is positive or negative.
**[Strong evidence]**

### Shadow 5 — physicalization requires two sections and one process

The smallest shadow of missing physics consists of three independent
elements: a scale section, an orientation/sector section, and an operational
process with preparation, instrument, and record. None substitutes for either
of the others. **[Proven]** as a separation of types.

### Shadow 6 — adaptation may act on memory rather than directly on a kernel

Instead of

\[
\dot K=-\nabla_K\mathcal L,
\]

the natural object may be

\[
\partial_t\Sigma(z;t)
=-\operatorname{grad}_{\Sigma}\mathcal L_{\rm process}.
\]

**[Speculative]** Positivity and complete monotonicity could constrain the
learning law.

### Shadow 7 — the bridge defect may be dynamic

Legacy and strict may differ not only by their distance profile but by hidden
memory structure. O13/G08 can test whether the defect decays, remains
constant, grows, or leaves the Stieltjes class. None of these outcomes is a
bridge by assumption; it is a new falsifiable atom.

## 11. Falsification of the strongest interpretations

### 11.1. Is G02 a new category theory?

No. **[Refuted]** Schur associativity and context diagrams are known
mathematics. Potential novelty lies in using them to organize FIN.

### 11.2. Does O05 solve the selector problem?

No. **[Refuted]** The family \(A(\theta)\) requires \(\theta\), and inversion
maps \(\theta\) to \(-\theta\), producing the pair \(\pm\Xi\). O05 is a new
receiver.

### 11.3. Is O09 physical entropy?

Not without additional data. **[Refuted]** Thermodynamics requires at least a
temperature, Hamiltonian, protocol, and an apparatus/reset ledger.

### 11.4. Does O12 prove a fractal fixed point?

No. **[Refuted]** O12 defines the right search space. Existence of an attractor
or self-similar operator measure remains open.

### 11.5. Does O11 prove that FIN is biological?

No. **[Refuted]** Neural operators, reservoir computing, and plasticity provide
mechanistic analogies, but no empirical map.

### 11.6. Does G09 provide a global section?

No on current artifacts. **[Refuted]** The sectional language exposes the
obstruction; it removes neither `QW-2191` nor the scale orbit.

## 12. Research programs 255–266

1. **P255 — Stieltjes theorem for all contexts.** Prove positivity, complete
   monotonicity, operator-measure representation, and the limits at zero and
   infinity. Probability: 0.95.
2. **P256 — Minimal realization of self-energy.** Classify hidden pairs
   \((A_{HH},A_{HE})\) producing the same \(\Sigma_E(z)\). Probability: 0.85.
3. **P257 — Formalize the context functor.** Define the category and verify
   Schur associativity in a proof assistant or independent algebraic core.
   Probability: 0.80.
4. **P258 — Analytic chiral memory susceptibility.** Derive \(\Xi_E(z)\),
   covariance, and norm bounds while separating gauge convention from flux
   observable. Probability: 0.82.
5. **P259 — Dynamic RG in Stieltjes space.** Define
   \(\mathcal R:\Sigma^{(n)}\mapsto\Sigma^{(n+1)}\) and search for fixed points,
   cycles, or no-go results. Probability: 0.45; high reward.
6. **P260 — Multitime information-balance tensor.** Extend O09 to a process
   ledger separating system, apparatus, and hidden-node contraction.
   Probability: 0.75.
7. **P261 — Dynamic completion-defect test.** Freeze one map and test
   \(\Delta_\Sigma(z)=\Sigma_{\rm strict}(z)-
   \mathcal C_\Sigma[\Sigma_{\rm legacy}(z)]\), without fitting after seeing
   the target. Bridge probability: 0.30; valuable no-go: 0.80.
8. **P262 — Calibrated Fingerprint Experiment.** Test spectral ratios first,
   then independently calibrate \(\tau\), with no post-unblinding changes.
   Methodological success: 0.75.
9. **P263 — Tomography of currents and chiral memory.** Measure
   \(C_d\), \(d=1,\ldots,5\), and intermediate interventions; test \(C_\chi\)
   and \(\Xi_E(z)\) separately. Probability: 0.70.
10. **P264 — False-positive atlas O01–O10.** Quantify how often random and
    structured Laplacians pass the FIN tests. Probability: 0.90.
11. **P265 — Identifiability of Adaptive Memory Geometry.** Compare genuine
    generator drift with static operators, apparatus memory, and calibration
    errors using hold-out data. Probability: 0.65.
12. **P266 — Biological/cybernetic benchmark.** Test fading-memory capacity,
    robustness–plasticity, recovery, and readout dimension against matched
    standard reservoirs. Probability: 0.72.

## 13. Strategic priority

| Rank | Program | Reason |
|---:|---|---|
| 1 | P255 | near-certain formal closure of the new object space |
| 2 | P256 | determines what can be inferred about the hidden sector |
| 3 | P258 | closes the new chiral-memory characteristic without selector claim |
| 4 | P257 | protects context language against metaphorical use |
| 5 | P264 | measures specificity of the atlas |
| 6 | P260 | joins information to a complete process |
| 7 | P262 | moves the result toward the laboratory |
| 8 | P263 | separates harmonics from memory experimentally |
| 9 | P259 | highest mathematical reward, high risk |
| 10 | P261 | one honest bridge attack |
| 11 | P265 | tests adaptation and identifiability |
| 12 | P266 | controls biological/cybernetic analogies |

## 14. Guardrail-compatible boundaries

This report does not export a non-premise strict selector, discharge
`QW-2191`, a canonical unit of length/time/mass/energy/action, a complete
legacy-to-strict map, role transfer for \(\sin^2\theta_W\),
\(\alpha_{\rm EM}\), or the gravity hierarchy, a strict-source adaptive law,
\(L_{\rm total}\), the Standard Model, GR, a ToE, a physical observer,
environment or apparatus, or external data validation.

The nadsoliton remains primordial information in a solitonic state. The
contexts \(E/H\) partition accessibility of its degrees of freedom; they are
not a separate informational layer beneath it.

## 15. Final verdict

The most important derived structure is neither another kernel nor another
constant. It is

\[
\boxed{\text{the Stieltjes--Schur Context Functor}},
\]

a contextual family of exact effective operators whose self-energies are
completely monotone, whose reductions compose, and whose time transforms
generate memory.

It organizes the common shadow

\[
\text{Green}\leftrightarrow\text{compression}
\leftrightarrow\text{memory}\leftrightarrow\text{observer}
\leftrightarrow\text{accessible information}.
\]

The most interesting new tangent characteristic is \(\Xi_E(z)\), the chiral
memory susceptibility. The most defensible account of information loss is
\(\mathcal L_t\), contraction of accessible distinguishability. The central
no-go remains the absence of orientation and scale sections.

This does not establish FIN as physics. It does provide a rigorous new object
space in which subsequent questions are better typed and easier to falsify.

## 16. Comparative literature

1. Dörfler, F.; Bullo, F., *Kron Reduction of Graphs with Applications to
   Electrical Networks*, <https://arxiv.org/abs/1102.2950>.
2. Lin, Y. T.; Tian, Y.; Anghel, M.; Livescu, D., *Data-driven learning for
   the Mori–Zwanzig formalism*, <https://arxiv.org/abs/2101.05873>.
3. Pollock, F. A. et al., *Non-Markovian quantum processes: complete framework
   and efficient characterisation*, <https://arxiv.org/abs/1512.00589>.
4. Marvian, I.; Spekkens, R. W., *Modes of asymmetry*,
   <https://arxiv.org/abs/1312.0680>.
5. Coifman, R. R.; Lafon, S., *Diffusion maps*,
   <https://doi.org/10.1016/j.acha.2006.04.006>.
6. Shuman, D. I. et al., *The Emerging Field of Signal Processing on Graphs*,
   <https://arxiv.org/abs/1211.0053>.
7. Kovachki, N. et al., *Neural Operator: Learning Maps Between Function
   Spaces*, <https://arxiv.org/abs/2108.08481>.
8. Abramsky, S.; Brandenburger, A., *The Sheaf-Theoretic Structure of
   Non-Locality and Contextuality*, <https://arxiv.org/abs/1102.0264>.
9. Chamseddine, A. H.; Connes, A., *The Spectral Action Principle*,
   <https://arxiv.org/abs/hep-th/9606001>.

## 17. Reproducible artifacts

- `fin_nadsoliton_second_generation_atlas.py`
- `FIN_Second_Generation_Combination_Search.csv`
- `FIN_Second_Generation_Derived_Objects.csv`
- `FIN_Second_Generation_Generation3.csv`
- `FIN_Second_Generation_Atlas_Results.json`
- `FIN_Second_Generation_Atlas_Figures/combination_search_coverage.png`
- `FIN_Second_Generation_Atlas_Figures/derived_object_score_matrix.png`
- `FIN_Second_Generation_Atlas_Figures/generation2_generation3_map.png`
- `FIN_Second_Generation_Atlas_Derived_Nadsoliton_Objects_EN.md`
