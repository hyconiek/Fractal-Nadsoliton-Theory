# Can FIN Become Predictive Physics?

## Twenty executed research programs, the wave–diffusion observer problem, and ten next recommendations

**Date:** 18 July 2026  
**Language:** English  
**Scope:** mathematical physics only; no cosmological, philosophical, or “Theory of Everything” claims  
**Repository state:** local HEAD 27b28fae, including the state-map through P3170  
**Reproducibility:** [core suite](./fin_ten_programs_experiments.py) and [observer-dynamics suite](./fin_observer_dynamics_experiments.py)

---

## Confidence convention

Every substantive conclusion carries exactly one required label:

- **[Proven]** — established by a proof, a standard theorem under stated hypotheses, or an exact finite computation.
- **[Strong evidence]** — supported by independent arguments and stable computations, but not a theorem covering all extensions.
- **[Moderate evidence]** — credible with important unresolved hypotheses.
- **[Weak evidence]** — suggestive but fragile.
- **[Speculative]** — a research possibility without adequate supporting evidence.
- **[Refuted]** — contradicted by proof or explicit counterexample in the stated scope.

Probabilities below are planning estimates, not mathematical confidence levels.

---

# Executive summary

**[Proven]** The strict and legacy kernels are exactly related on the fixed twelve-cycle: each is a degree-at-most-six polynomial in the cycle Laplacian and, because their seven radial spectral values are separately distinct, each is a polynomial in the other.

**[Proven]** This finite bridge is not a continuum bridge. Its raw interpolation matrix has condition number about \(8.0\times10^6\), and applying the same polynomial to natural \(C_N\) extensions gives relative errors \(9.17\times10^2\), \(1.32\times10^4\), and \(6.22\times10^5\) at \(N=16,24,48\). It is seven-point interpolation, not a robust completion law.

**[Refuted]** Different inertia does not forbid a real finite functional relation between the kernels. A real polynomial may cross zero; the exact bridge is a counterexample. A positivity-preserving or Stieltjes bridge may still be obstructed, but that stronger requirement must be stated.

**[Proven]** The reflection-even operator algebra cannot produce an absolute orientation. Yet a system–apparatus pair has a complete relative invariant

\[
(\Delta,c)=
\bigl(\lambda_a(r_s-r_a)\bmod12,\lambda_s\lambda_a\bigr).
\]

Exhaustive enumeration of all \(576\) frame pairs gives exactly \(24\) diagonal dihedral orbits, each of size \(24\). Relational physics can therefore avoid an absolute selector without claiming that the core generates one.

**[Proven]** Dimensionless invariant data cannot naturally generate a nonzero dimensionful quantity. The full action–length–time span requires three independent calibration scales; \(15\,624\) monomials in six weight-zero receivers contain no weight-one source.

**[Proven]** A single \(C_{12}\) operator does not determine its continuum. Three scaling families agreeing at \(N=12\) have first gaps tending to zero, a nonzero constant, and infinity. The cycle Laplacian has quadratic low-momentum dispersion; its square root has linear dispersion.

**[Proven]** The same \(L\) supports distinguishable unitary and diffusive dynamics. At optimized dimensionless time \(6.77\), their position distributions have Jensen–Shannon divergence \(0.363\) nats and \(D_{\rm KL}=1.392\) nats; ideally eight observations give expected log Bayes factor at least ten.

**[Proven]** More strongly, the actual strict QW2118 matrix has positive off-diagonal weights and constant row sum \(s=1.660307278766099\). The Perron shift \(A=sI-K_s\) is simultaneously a positive Dirichlet generator for the unitary group \(e^{-itA}\) and the reversible Markov semigroup \(e^{-tA}\).

**[Proven]** These are inequivalent temporal laws, not contradictory descriptions. Coherent escape is quadratic at short time, diffusion is linear, and the Born population matrices violate Chapman–Kolmogorov while the Markov matrices obey it. The operator does not select the physical temporal ray, state cone, instrument, environment, or clock.

**[Proven]** An exact infinite-dimensional unitary dilation reproduces FIN diffusion for the restricted observer, but no finite closed bath can reproduce nonconstant exponential relaxation for all time. Frequent projective observation freezes the coherent walk instead of automatically creating the FIN heat flow.

**[Proven]** Shannon information does not determine energy or temperature. A two-level SWAP process verifies the finite Landauer equality to \(8.3\times10^{-17}\), but only after a Hamiltonian, Gibbs state, temperature, and process are supplied.

**[Refuted]** The naive tuple

\[
(\mathbb C^{12},\mathbb C^{12},D=K_s,J=\text{complex conjugation})
\]

is not a standard real/even spectral triple: the maximum point-projection first-order residual is \(1.0373\), and a compatible even grading is obstructed.

**[Refuted]** \(K_s\) is not directly a stable Euclidean Gaussian action: it has seven negative modes. A shift greater than \(0.68187476238\) is necessary and is ill-conditioned near threshold.

**[Strong evidence]** The realistic route to physics is a conjunction of a refinement/causal limit, a calibrated operational realization, scale-free held-out predictions, and formal proof certificates. Another dimensionless scalar cannot supply these structures.

---

# Part I — Independent structural analysis

## 1. Method

The program ranking was not inherited from prior suggested future work. The procedure was:

1. reconstruct the frozen mathematical data;
2. reduce repository blockers to logically distinct open-problem classes;
3. score candidate programs by expected scientific value;
4. select ten without a diversity quota;
5. execute proofs, counterexamples, computations, and literature comparisons;
6. retain negative results at equal status with positive ones;
7. require any physical proposal to name preparation, dynamics, measurement, calibration, and falsification.

The independent calculations used

\[
K_s(d)=\frac{\cos(0.18575d+0.16250)}{1+d^{1.8}},
\qquad
K_l(d)=4\ln2\,
\frac{\cos(\pi d/4+\pi/6)}{1+0.01d},
\]

with zero diagonal and cyclic distance on \(C_{12}\).

**[Proven]** The deterministic NumPy suite completed without importing generated FIN verdicts.

## 2. Structural open-problem register

“Complete” means complete at the level of dependency classes; historical microtasks are instances of these classes.

| ID | Open problem | Importance | Physical relevance | Status |
|---:|---|---:|---:|---|
| O1 | canonical and stable strict–legacy bridge | 8 | 5 | finite polynomial proved; stability open |
| O2 | parameter-running provenance beyond interpolation | 8 | 7 | open |
| O3 | legacy physical-role transfer | 5 | 8 | blocked on a physical bridge |
| O4 | uniqueness of the selected filters inside \(C^*(L)\) | 9 | 8 | open |
| O5 | out-of-sample parameter identifiability | 7 | 8 | open |
| O6 | absolute origin/orientation selector | 8 | 7 | no-go proved |
| O7 | relational selector for all observables | 9 | 10 | theorem proved; audit open |
| O8 | odd-sector adaptive bifurcation | 9 | 8 | open |
| O9 | internal dimensionful source | 9 | 10 | invariant-input no-go |
| O10 | minimal calibration package | 6 | 10 | rank-three package proved |
| O11 | physical clock and time flow | 9 | 10 | open |
| O12 | noncircular adaptive potential \(V\) | 9 | 7 | open |
| O13 | self-referential fixed points and stability | 9 | 6 | open |
| O14 | FIN-native refinement \(K_N\) | 10 | 10 | open |
| O15 | norm/strong-resolvent continuum limit | 10 | 10 | open |
| O16 | causal order or propagation cone | 10 | 10 | exact finite cone absent |
| O17 | Lorentzian limiting principal symbol | 10 | 10 | open |
| O18 | physical state preparation | 8 | 10 | absent |
| O19 | instruments and apparatus outcomes | 8 | 10 | absent |
| O20 | classical/quantum operational probability law | 9 | 10 | absent |
| O21 | bath, KMS state, temperature, allowed operations | 8 | 9 | conditional only |
| O22 | bounded variational action selecting \(K_s\) | 9 | 8 | shift/potential required |
| O23 | Osterwalder–Schrader positivity | 10 | 10 | undefined before temporal refinement |
| O24 | metric family and stress–energy variation | 9 | 10 | absent |
| O25 | repaired real/even spectral triple | 9 | 8 | naive version refuted |
| O26 | nonabelian gauge group and representation | 9 | 10 | absent |
| O27 | fermions, chirality, KO dimension | 10 | 10 | absent |
| O28 | held-out dimensionless predictions | 7 | 10 | candidate fingerprint computed |
| O29 | physical twelve-mode implementation | 6 | 10 | open |
| O30 | robustness to noise and graph changes | 8 | 9 | open |
| O31 | machine-certified no-go/source proofs | 8 | 5 | feasible |
| O32 | justified ensemble or universality limit | 7 | 5 | absent |
| O33 | selected subsystem/tensor factorization | 6 | 5 | choice underdetermined |
| O34 | theorem-level originality | 7 | 4 | core largely standard |

**[Strong evidence]** O7, O14–O20, O23, O28, and O29 have the largest capacity to change FIN’s epistemic status.

## 3. Candidate ranking and selected programs

Expected scientific value combines depth, originality, robustness, physical relevance, testability, feasibility, and redundancy. It is a decision score, not a theorem.

| Rank | Direction | Score /100 | Decision |
|---:|---|---:|---|
| 1 | calibrated prediction and Bayesian discrimination | 82 | Program 9 |
| 2 | refinement, continuum, causality, Lorentz | 79 | Program 4 |
| 3 | relational selector/reference frames | 77 | Program 2 |
| 4 | exact bridge plus generalization falsification | 74 | Program 1 |
| 5 | adaptive bootstrap and dynamical identifiability | 70 | Program 5 |
| 6 | variational stability and OS gate | 67 | Program 8 |
| 7 | units, RG, modular time | 65 | Program 3 |
| 8 | information thermodynamics/open systems | 62 | Program 6 |
| 9 | repaired noncommutative geometry | 58 | Program 7 |
| 10 | machine-assisted invariant/source proofs | 56 | Program 10 |
| 11 | tensor-network compression | 43 | auxiliary only |
| 12 | information geometry/optimal transport | 41 | not selected |
| 13 | topological phases/\(K\)-theory | 36 | not selected |
| 14 | causal-set reinterpretation | 33 | not selected |
| 15 | random-matrix robustness | 31 | not selected |
| 16 | sheaves/derived geometry/HoTT | 28 | language, not generator |
| 17 | reverse-engineered SM/GR receivers | 18 | too circular |
| 18 | free probability of the present pair | 0 | refuted |

**[Proven]** The present commuting non-scalar kernel pair cannot be freely independent in a faithful trace: freeness gives \(\tau(abab)=0\), while commutation and freeness give \(\tau(abab)=\tau(a^2)\tau(b^2)>0\).

---


# Part II — Ten executed research programs

# Program 1 — Exact spectral bridge and its generalization test

## 1. Motivation

The first question is whether the strict and legacy kernels are genuinely different mathematical objects or merely different coordinates on the same finite spectral datum. A robust answer would clarify what can be transported between them and what cannot.

**Research hypothesis.** On the fixed cycle they are related by functional calculus, but no size-independent completion law follows from that fact. **[Strong evidence]**

The expected breakthrough was either a canonical bridge or a precise theorem exposing why the apparent bridge is only interpolation. Estimated success probability was \(0.85\), impact \(8/10\), difficulty \(6/10\), and research time one to three months for a formal stability paper. The program is logically prior to role transfer, but not to the operational projects.

## 2. Background

Let \(L_N\) be the combinatorial Laplacian of \(C_N\). Every real radial circulant is diagonal in the Fourier basis and belongs to \(C^*(L_N)\). On \(C_{12}\), radial symmetry leaves seven distinct Laplacian eigenvalues,

\[
\ell_k=4\sin^2(\pi k/12),\qquad k=0,\ldots,6.
\]

The spectral theorem therefore already implies that any such kernel is \(p(L_{12})\) for a real polynomial \(p\) of degree at most six. **[Proven]** This is standard finite functional calculus, not a FIN-specific new theorem.

## 3. Mathematical formulation

Write

\[
K_s=p_s(L_{12}),\qquad K_l=p_l(L_{12}).
\]

The stronger bridge problem asks for functions \(q,r\) such that

\[
K_s=q(K_l),\qquad K_l=r(K_s).
\]

Such functions exist exactly when each target eigenvalue is constant on every degeneracy class of the source. A canonical physical bridge would additionally need stability, a declared admissible function class, compatibility with refinement, and preservation of any role-bearing structure.

## 4. Research methodology

The calculation proceeded in four stages:

1. build both \(12\times12\) circulants from their declared radial profiles;
2. diagonalize them in the common Fourier basis;
3. solve the seven-point Vandermonde systems for \(p_s,p_l,q,r\);
4. freeze the polynomial \(q\), extend the two profile formulas naturally to \(C_N\), and test \(q(K_{l,N})\approx K_{s,N}\) at unseen sizes.

Residuals were evaluated spectrally and in Frobenius norm. Conditioning and minimum eigenvalue separation were recorded to expose interpolation fragility.

## 5. Detailed derivations

The interpolation lemma is elementary but decisive.

**Theorem 1.** If a self-adjoint matrix \(A\) has \(m\) distinct eigenvalues and a commuting self-adjoint matrix \(B\) is scalar on every eigenspace of \(A\), then \(B=q(A)\) for a real polynomial of degree at most \(m-1\). **[Proven]**

**Proof.** Write \(A=\sum_{j=1}^m a_jP_j\) and \(B=\sum_{j=1}^m b_jP_j\). Lagrange interpolation gives a real polynomial \(q\) with \(q(a_j)=b_j\), hence \(q(A)=\sum_jq(a_j)P_j=B\). \(\square\)

Both kernels have seven distinct radial eigenvalues. Therefore

\[
K_s=q(K_l),\qquad K_l=r(K_s)
\]

for polynomials of degree at most six. The computed residuals were

| identity | relative residual |
|---|---:|
| \(K_s-p_s(L)\) | \(2.53\times10^{-14}\) |
| \(K_l-p_l(L)\) | \(1.28\times10^{-13}\) |
| \(K_s-q(K_l)\) | \(7.80\times10^{-14}\) |
| \(K_l-r(K_s)\) | \(2.91\times10^{-13}\) |

These are floating-point representations of exact finite interpolation identities. **[Proven]**

The minimum separations of the seven source nodes are \(0.04358\) for the strict spectrum and \(0.43466\) for the legacy spectrum; the raw Laplacian Vandermonde condition number is \(8.04\times10^6\). Thus coefficient values depend strongly on basis and perturbation even though the evaluated finite identity is exact. **[Strong evidence]**

## 6. Proof attempts

An attempted obstruction based on different inertia fails. A real polynomial is not sign preserving, so \(A\) and \(q(A)\) may have arbitrary sign patterns on a finite spectrum. The exact \(q\) above is a counterexample to the unrestricted inertia claim. **[Refuted]**

A stronger theorem remains possible: a positive, completely monotone, operator-monotone, Stieltjes, local, or parameter-sparse bridge may be incompatible with the two spectral orderings. That problem must name its function class before inertia or monotonicity can be used. **[Moderate evidence]**

No proof was found that the interpolation polynomial is selected by the kernel formulas, by an action, or by a functor across graph sizes. **[Strong evidence]**

## 7. Counterexamples

The fixed-size theorem is maximally nonunique outside seven nodes: for any polynomial \(h\),

\[
q_h(x)=q(x)+h(x)\prod_{j=1}^{7}(x-\lambda_j(K_l))
\]

defines the same \(K_s\) on \(C_{12}\). **[Proven]** Thus the matrix identity alone selects neither extrapolation nor physical semantics.

Role transfer also fails logically: if two matrices are related by a polynomial, labels such as energy, coupling, or force do not transfer unless the preparation, action, observables, and units are intertwined. **[Proven]**

## 8. Numerical investigations

The degree-six bridge fitted at \(N=12\) was evaluated on natural profile extensions. Relative Frobenius errors were

| unseen size | \(\|q(K_{l,N})-K_{s,N}\|_F/\|K_{s,N}\|_F\) |
|---:|---:|
| 16 | \(9.17\times10^2\) |
| 24 | \(1.32\times10^4\) |
| 48 | \(6.22\times10^5\) |

The extrapolation is not merely inaccurate; it diverges rapidly for this natural continuation. **[Proven]** This does not prove that every possible refinement bridge fails, but it falsifies the claim that finite functional equivalence supplies one automatically.

## 9. Literature comparison

Finite circulant functional calculus is a routine consequence of the spectral theorem. Graph spectral filters are likewise functions of graph Laplacians; the modern signal-processing formulation is exemplified by [Hammond, Vandergheynst and Gribonval](https://arxiv.org/abs/0912.3848). The FIN-specific content is the exact two-way interpolation and its failed out-of-size generalization, not a new spectral theorem.

## 10. Final conclusions

The strict and legacy kernels are two coordinates on one commutative seven-point spectral algebra at \(N=12\). **[Proven]**

They are not thereby one refinement law, one physical operator, or one role-bearing object. **[Proven]**

The most useful continuation is to search for a low-complexity bridge stable under a separately justified \(K_N\) family, not to repeat fixed-size interpolation. **[Strong evidence]**

## 11. Confidence assessment

- Exact finite polynomial equivalence: **[Proven]**
- Inertia-only obstruction to real functional equivalence: **[Refuted]**
- Existing polynomial as a canonical completion map: **[Refuted]**
- Existence of some sparse, positive, refinement-stable bridge: **[Speculative]**
- Estimated probability of such a bridge: \(0.20\)
- Scientific impact if found: \(8/10\)
- Difficulty: \(8/10\)
- Time: 6–18 months after a refinement family exists
- Dependencies: Program 4 and a precise admissible bridge category

---

# Program 2 — Relational selector and reference-frame completion

## 1. Motivation

The current core is reflection symmetric, so an absolute clockwise/counterclockwise choice is obstructed. Physics, however, normally measures a system relative to an apparatus. The program asks whether the selector obligation should be replaced by a relational theorem instead of solved by an impossible absolute construction.

**Research hypothesis.** The unoriented kernel cannot select a frame, but a system–apparatus pair admits complete invariant relative coordinates. **[Strong evidence]**

Estimated probability of mathematical success was \(0.95\), probability of resolving the physical selector problem \(0.75\), impact \(10/10\), difficulty \(5/10\), and time one to three months for a complete downstream audit.

## 2. Background

Let reflection act by \((R\psi)(x)=\psi(-x)\). Radial \(L,K_s,K_l\) satisfy \(RKR=K\). Functional calculus, spectral measures, traces, entropy, and any natural automorphism-preserving construction remain even.

The oriented-frame torsor is

\[
F=\mathbb Z_{12}\times\{\pm1\}.
\]

The dihedral group acts freely and transitively by

\[
(a,\varepsilon)\cdot(r,\lambda)
=(a+\varepsilon r,\varepsilon\lambda).
\]

## 3. Mathematical formulation

An absolute selector is an equivariant map from the current invariant input to \(F\) or to the orientation torsor \(\{+,-\}\). A relational selector is instead an invariant of a pair \((f_s,f_a)\in F\times F\) under the simultaneous diagonal action.

For \(f_s=(r_s,\lambda_s)\) and \(f_a=(r_a,\lambda_a)\), define

\[
\Delta=\lambda_a(r_s-r_a)\pmod {12},
\qquad c=\lambda_s\lambda_a.
\]

## 4. Research methodology

The work combined:

1. a naturality/no-fixed-point proof;
2. representation decomposition into reflection-even and reflection-odd operators;
3. an analytic orbit-classification theorem for \(F\times F\);
4. exhaustive enumeration of all \(576\) pairs and all \(24\) dihedral transformations;
5. a minimal Gibbs interaction and an odd-current operational probe.

## 5. Detailed derivations

**Theorem 2.** No automorphism-natural absolute orientation can be constructed from the current even FIN core. **[Proven]**

**Proof.** If input \(D\) is fixed by reflection and \(S\) is equivariant, then \(S(D)=S(RD)=R S(D)\). Neither orientation is reflection fixed, a contradiction. \(\square\)

Equivalently, if

\[
J=\frac{T-T^*}{2i},\qquad RJR=-J,
\]

then

\[
C^*(I,L,K_s,K_l)\cap\{A:RAR=-A\}=\{0\}.
\]

**Theorem 3.** The pair \((\Delta,c)\) is a complete invariant of the diagonal \(D_{12}\)-action on \(F\times F\); equivalently,

\[
(F\times F)/D_{12}\cong D_{12}.
\]

**[Proven]**

**Proof.** In torsor notation the invariant is \(f_a^{-1}f_s\). It is unchanged by simultaneous left multiplication. Conversely, equality \(f_a^{-1}f_s=f_a'^{-1}f_s'\) gives the unique group element \(g=f_a'f_a^{-1}\) carrying both members of the first pair to the second. The coordinate expression is exactly \((\Delta,c)\). \(\square\)

Every diagonal-invariant observable therefore factors through these 24 relative values. **[Proven]**

## 6. Proof attempts

Index theory, eta invariants, cohomology, categorical fixed points, and spectral flow were tested as potential hidden selectors. They cannot help while their inputs and paths remain reflection even: natural constructions inherit the automorphism. **[Proven]**

The cycle has \(H^1(C_{12};\mathbb R)\cong\mathbb R\), but reflection acts by \(-1\), so the dihedral-invariant part is zero. Orientability supplies two sections, not a preferred one. **[Proven]**

An internal absolute selector could exist only after adding a nonzero odd state, source, boundary condition, current, chiral density, or thermodynamic sector. **[Proven]**

## 7. Counterexamples

Choosing vertex zero, the first stored eigenvector, or a lexicographic sign returns a unique frame but changes under relabelling. Such routines are conventions, not natural selectors. **[Refuted]**

A symmetric pitchfork produces a pair of branches \(\pm m\), not a canonical branch. At every finite size and zero external field the invariant equilibrium mean remains zero. **[Proven]**

The relational solution is insufficient if a downstream claim explicitly depends on an apparatus-independent absolute sign. In that case the missing odd resource remains a genuine extra axiom. **[Proven]**

## 8. Numerical investigations

Enumeration produced 24 relative invariants, each fibre containing 24 of the 576 frame pairs, and passed invariance under all 24 group elements. **[Proven]**

For orientation signs \(s,a\in\{\pm1\}\), the invariant interaction

\[
H_{\rm int}=-gsa
\]

at \(\beta g=2\) gives

\[
\langle s\rangle=\langle a\rangle=0,
\quad \langle sa\rangle=\tanh2=0.96402758,
\quad P(s=a)=0.98201379.
\]

Reliable relative alignment coexists with no absolute orientation. **[Proven]**

For a calibrated odd probe \(H(h)=K_s-hJ\),

\[
E_k(h)-E_{12-k}(h)=-2h\sin(2\pi k/12),
\]

and the computed Gibbs current was \(-0.105343,0,+0.105343\) at \(h=-0.2,0,+0.2\) for the declared numerical setting. **[Strong evidence]**

## 9. Literature comparison

This is the finite torsor form of operational quantum reference frames, where invariant relations replace inaccessible absolute quantities. See [Gour and Spekkens](https://arxiv.org/abs/0711.0043), [Loveridge, Miyadera and Busch](https://arxiv.org/abs/1703.10434), and [Carette, Głowacki and Loveridge](https://arxiv.org/abs/2303.14002). The group-quotient theorem is standard; the FIN-specific result is its exact application as a replacement test for the selector.

## 10. Final conclusions

An absolute selector is impossible from the present even core. **[Proven]**

A complete relational frame exists once an apparatus is included. **[Proven]**

The scientifically decisive audit is whether every proposed observable can be expressed through \((\Delta,c)\). If yes, the absolute selector is gauge; if not, a physical odd resource must be supplied and measured. **[Strong evidence]**

## 11. Confidence assessment

- Absolute natural-selector no-go: **[Proven]**
- Complete relational-frame theorem: **[Proven]**
- Replacement of every FIN absolute selector claim: **[Moderate evidence]**
- Internal absolute direction from present scalar data: **[Refuted]**
- Estimated probability that relational reformulation closes the physical selector problem: \(0.75\)
- Impact: \(10/10\)
- Difficulty: \(6/10\)
- Time: 1–3 months for observable-by-observable rewriting
- Dependencies: Program 9's apparatus model and Program 5's dynamics

---

# Program 3 — Dimensional covariance, RG, and modular time

## 1. Motivation

FIN is dimensionless. The program tests whether heat asymptotics, zeta regularization, renormalization, or modular flow can create action, length, time, mass, and energy without inserting units by hand.

**Research hypothesis.** These constructions can organize relative scales but cannot select a nonzero dimensionful standard from invariant input. **[Strong evidence]**

The expected negative theorem has high foundational value: it prevents hidden unit insertion. Estimated success probability \(0.95\), impact \(9/10\), difficulty \(7/10\), time two to six months, with Program 4 providing any future refinement needed to reopen the question.

## 2. Background

A physical quantity transforms under the positive unit-rescaling group

\[
G=(\mathbb R_{>0})^r.
\]

Dimensionless spectra, ratios, entropies, graph labels, and phases lie in the trivial representation. A nonzero length, time, or action lies in a nontrivial one-dimensional representation.

Heat kernels and zeta functions extract asymptotic information only from a family with a scale regime. Renormalization requires a coarse-graining map and reference condition. Tomita–Takesaki flow requires a von Neumann algebra and faithful state.

## 3. Mathematical formulation

Let \(X\) be the total dimensionless FIN input and \(V_\chi\) the one-dimensional representation with character \(\chi\ne1\). The source question is whether a natural equivariant map

\[
F:X\longrightarrow V_\chi\setminus\{0\}
\]

exists.

For a calibrated completion, use independent dimensions action–length–time and ask for the minimum number of charged scales spanning

\[
[\hbar_*]=(1,0,0),\quad
[\ell_*]=(0,1,0),\quad
[\tau_*]=(0,0,1).
\]

## 4. Research methodology

The program used:

1. a representation-theoretic equivariance proof;
2. exact rank tests for one, two, and three calibration sources;
3. exhaustive integer-weight enumeration of \(15\,624\) receiver monomials;
4. heat-dimension and zeta-determinant calculations on the shifted strict spectrum;
5. two exact Gaussian coarse-graining maps;
6. an algebraic audit of modular time on \(\mathbb C^{12}\) and \(M_{12}(\mathbb C)\).

## 5. Detailed derivations

**Theorem 4 — dimensional-source obstruction.** If \(G\) acts trivially on \(X\) and through nontrivial \(\chi\) on \(V_\chi\), every equivariant map \(F:X\to V_\chi\) is zero. **[Proven]**

**Proof.** For every \(g\),

\[
F(x)=F(gx)=gF(x)=\chi(g)F(x).
\]

Choose \(g\) with \(\chi(g)\ne1\). Then \(F(x)=0\). \(\square\)

Thus Shannon numbers, eigenvalue ratios, phases, determinants, and dimensionless couplings cannot naturally output a nonzero SI-valued quantity. **[Proven]**

Three independent calibrations define

\[
E_*=\frac{\hbar_*}{\tau_*},\qquad
p_*=\frac{\hbar_*}{\ell_*},\qquad
c_*=\frac{\ell_*}{\tau_*},\qquad
m_*=\frac{\hbar_*\tau_*}{\ell_*^2}.
\]

Their weight matrix has rank three. Removing any one leaves at least one of action, length, time, energy, or mass underdetermined. **[Proven]**

This is a minimal generic calibration package, not a proof that nature supplies those values from FIN.

## 6. Proof attempts

For \(A=K_s-\lambda_{\min}I\ge0\),

\[
Z(t)=\operatorname{Tr}e^{-tA},\qquad
d_s(t)=\frac{2t\,\operatorname{Tr}(Ae^{-tA})}{Z(t)}.
\]

Every finite matrix satisfies \(d_s(t)\to0\) as \(t\downarrow0\), and with a zero mode also as \(t\to\infty\). Its finite zeta function is entire and has no Weyl pole. **[Proven]** A transient heat dimension is therefore not an emergent physical dimension.

Dimensional transmutation has the form

\[
\Lambda_{\rm RG}=\mu
\exp\!\left[-\int^g\frac{dg'}{\beta(g')}\right].
\]

It converts a coupling specified at reference scale \(\mu\) into another scale; it does not eliminate \(\mu\). **[Proven]**

On the commutative algebra \(\mathbb C^{12}\), every faithful state is tracial and modular flow is trivial. On \(M_{12}\), a supplied state \(\rho_\beta\propto e^{-\beta K}\) gives

\[
\sigma_t^\rho(A)=e^{-it\beta K}Ae^{it\beta K},
\]

but arbitrary \(\beta\) rescales the clock. **[Proven]**

## 7. Counterexamples

Three graph families can agree at \(N=12\) while using \(A_N\), \(N^2A_N\), or \(N^4A_N\); the first nonzero gaps then tend to zero, a constant, or infinity. **[Proven]** The finite spectrum cannot choose a scaling orbit representative.

A code that returns “one unit” does not solve the problem: under a change from metres to centimetres the physical numerical value must transform, while an invariant algorithm receives unchanged input. **[Refuted]**

An RG blocking rule is not encoded by \(K\). Decimation and pair averaging are equally definable and yield different effective operators. **[Proven]**

## 8. Numerical investigations

The exhaustive weight search found zero weight-one monomials among \(15\,624\) products of six weight-zero receivers. Rank for one, two, and three independent calibration scales was \(1,2,3\). **[Proven]**

The shifted strict spectrum has a maximum transient heat dimension about \(1.09178\) near \(t=20.69\), but zero ultraviolet and infrared limits. **[Proven]**

For \(A=K_s+0.8I\), exact Gaussian decimation and normalized pair averaging produced effective precision matrices separated by Frobenius distance \(2.47996\). Scheme independence is absent at this finite stage. **[Proven]**

## 9. Literature comparison

Zeta regularization defines determinants after an operator scale is given; see [Hawking](https://doi.org/10.1007/BF01626516). The [spectral action](https://arxiv.org/abs/hep-th/9606001) explicitly contains its cutoff. Wilsonian RG explicitly begins with a coarse-graining and scale convention; see [Wilson and Kogut](https://doi.org/10.1016/0370-1573(74)90023-4). The thermal-time proposal of [Connes and Rovelli](https://arxiv.org/abs/gr-qc/9406019) presupposes an algebra and state. These frameworks offer conditional physicalization tools, not an exception to the equivariance theorem.

## 10. Final conclusions

No internal nonzero physical unit can be generated from the current dimensionless invariant core. **[Proven]**

A calibrated theory can be built after supplying three independent action–length–time scales or an equivalent rank-three package. **[Proven]**

Heat, zeta, RG, and modular constructions may become useful only after a refinement, state, blocking prescription, and calibration are supplied. **[Strong evidence]**

## 11. Confidence assessment

- Equivariant dimensional-source no-go: **[Proven]**
- Generic rank-three calibration package: **[Proven]**
- Absolute modular clock from the current core: **[Refuted]**
- FIN-native dimensional transmutation after a new refinement: **[Speculative]**
- Probability that current data alone determine a unit: \(0\)
- Probability that a calibrated extension is mathematically coherent: \(0.95\)
- Impact: \(9/10\)
- Difficulty: \(7/10\)
- Time: 2–6 months for a complete unit-covariant formalization
- Dependencies: Programs 4, 6, and 9

---

# Program 4 — Refinement, continuum, causality, and Lorentz structure

## 1. Motivation

A physical theory needs a controlled large-system or continuum limit. The isolated twelve-point operator could be a discretization of many inequivalent objects, so the program asks whether FIN itself selects a refinement, propagation law, causal order, or Lorentzian principal symbol.

**Research hypothesis.** A continuum can be constructed after adding a refinement package, but no unique limit follows from the \(C_{12}\) snapshot. **[Strong evidence]**

The target theorem would have the highest physical impact, \(10/10\), difficulty \(10/10\), a 6–24 month horizon, and an estimated \(0.25\) probability for a natural refinement and \(0.10\) for a causal Lorentzian limit.

## 2. Background

Discrete Laplacians can converge to continuum Laplace operators when meshes, embeddings, measures, and normalizations are specified. Strong-resolvent or Mosco convergence then transports spectra and semigroups. None of these data is contained in one finite matrix.

Lorentzian physics additionally needs a hyperbolic local principal symbol, a time orientation, and a causal interpretation. The inertia of a global matrix is not a spacetime metric signature.

## 3. Mathematical formulation

A sufficient target object is a sequence

\[
(\mathcal A_N,H_N,D_N,\iota_N,a_N,\tau_N,\Gamma_N)
\]

with refinement maps \(\iota_N\), spatial and temporal scales \(a_N,\tau_N\), measures, an orientation datum \(\Gamma_N\), controlled locality, and a declared convergence topology.

A FIN-native theorem would need to derive parameter running for at least \(\omega_N,\beta_N,\eta_N\), the diagonal/normalization convention, and the measure.

## 4. Research methodology

The program used:

1. a nonuniqueness construction for finite prefixes;
2. the cycle Laplacian as a positive convergence control;
3. mutually incompatible extensions agreeing exactly at \(N=12\);
4. low-momentum regression for \(L_N\) and \(\sqrt{L_N}\);
5. short-time propagation tests for the strict long-range matrix;
6. algebraic separation of matrix inertia from Lorentz signature.

## 5. Detailed derivations

**Theorem 5 — single-snapshot obstruction.** A fixed finite operator is compatible with infinitely many refinement sequences having mutually inequivalent limits or no limit. **[Proven]**

**Proof.** Keep \(D_{12}\) fixed. For \(N>12\), choose cycle Laplacians converging to \(S^1\), grid Laplacians converging to higher-dimensional tori, direct sums converging to disconnected spaces, or an unbounded oscillating sequence. Their common finite member imposes no asymptotic relation. \(\square\)

For the cycle Laplacian,

\[
\lambda_k(L_N)=4\sin^2(\pi k/N),
\]

and \(N^2L_N\) converges spectrally at fixed \(k\) to the unit-circumference circle Laplacian with eigenvalues \((2\pi k)^2\). **[Proven]**

| \(N\) | relative error, \(k=1\) | \(k=2\) | \(k=3\) |
|---:|---:|---:|---:|
| 12 | \(-2.264\%\) | \(-8.811\%\) | \(-18.943\%\) |
| 24 | \(-0.570\%\) | \(-2.264\%\) | \(-5.036\%\) |
| 48 | \(-0.143\%\) | \(-0.570\%\) | \(-1.279\%\) |
| 96 | \(-0.0357\%\) | \(-0.143\%\) | \(-0.321\%\) |
| 192 | \(-0.00892\%\) | \(-0.0357\%\) | \(-0.0803\%\) |

This proves that a continuum is possible after the missing data are supplied, not that FIN selects them.

## 6. Proof attempts

Two equally simple strict-profile continuations agree at \(N=12\): fixed lattice parameters, and running parameters preserving physical wavelength and damping length,

\[
\omega_N=\omega_{12}\frac{12}{N},\qquad
\beta_N=\beta_{12}\left(\frac{12}{N}\right)^\eta.
\]

At fixed physical separation \(x=1/4\), the first tends toward zero with oscillations while the second stays at \(0.0914286\). **[Proven]** No canonicality argument selected one.

Low-momentum fits gave power \(1.999999\) for \(L_N\) and \(0.9999996\) for \(\sqrt{L_N}\). The former suggests diffusive/Schrödinger \(z=2\); the latter wave-like \(z=1\). The spatial operator alone does not select temporal order or dynamic exponent. **[Proven]**

## 7. Counterexamples

For strict \(K_{ij}\ne0\) at every cyclic separation on \(C_{12}\),

\[
e^{-itK}_{ij}=-itK_{ij}+O(t^2).
\]

Every vertex therefore acquires nonzero amplitude at arbitrarily small time; there is no exact finite causal cone. **[Proven]**

Orienting all cycle edges produces a directed closed loop, not an acyclic causal order. **[Proven]**

The seven-negative/five-positive inertia of \(K_s\) lives in Fourier-mode space. A Lorentz metric is a local field whose principal symbol has one timelike direction at every point. Equating them is a category error. **[Refuted]**

## 8. Numerical investigations

The first gap of \(L_N\), \(N^2L_N\), and \(N^4L_N\) illustrates three incompatible scalings. At \(N=192\) the respective gaps are

\[
0.00107083,\qquad 0.274131,\qquad 70.1776.
\]

All three constructions contain exactly the same unscaled \(L_{12}\) at the reference size when normalization is declared only afterward. **[Proven]**

The farthest strict coupling is small but nonzero:

\[
K_s(6)/K_s(1)\approx0.0236.
\]

This may permit approximate long-range propagation bounds in a justified family, but it cannot yield exact microcausality on the present graph. **[Moderate evidence]**

## 9. Literature comparison

Discrete-to-continuum Laplacian convergence requires explicit triangulations and normalization; see [Dodziuk](https://doi.org/10.2307/2373615) and modern graph-like convergence work by [Post and Simmer](https://arxiv.org/abs/1704.00064). Approximate cones for local quantum systems descend from [Lieb and Robinson](https://doi.org/10.1007/BF01645779). Causal order determines at most conformal Lorentz geometry under strong manifold hypotheses; see [Malament](https://doi.org/10.1063/1.523436). FIN currently supplies none of those hypotheses uniquely.

## 10. Final conclusions

A positive continuum construction exists as a conditional mathematical route. **[Proven]**

The present finite object neither selects the refinement nor supplies exact causality, Lorentz signature, time orientation, or physical scale. **[Proven]**

The next admissible high-value theorem is a refinement functor with running parameters and convergence guarantees, not another interpretation of the \(N=12\) inertia. **[Strong evidence]**

## 11. Confidence assessment

- Single-snapshot continuum obstruction: **[Proven]**
- Conditional convergence of normalized cycle Laplacians: **[Proven]**
- Exact causal cone on current strict matrix: **[Refuted]**
- FIN-native natural refinement: **[Speculative]**
- Lorentzian scaling limit: **[Speculative]**
- Probability of a useful refinement theorem: \(0.25\)
- Probability of a causal Lorentzian limit: \(0.10\)
- Impact: \(10/10\)
- Difficulty: \(10/10\)
- Time: 6–24 months
- Dependencies: Programs 3, 7, 8, 11–14, and explicit \(N\)-dependent data

---

# Program 5 — Adaptive bootstrap and dynamical identifiability

## 1. Motivation

FIN includes adaptive or self-referential language. The program asks whether the proposed learning law is genuinely variational, whether its fixed points select the strict kernel, and whether observed dynamics identify a unique update law rather than merely fit one.

**Research hypothesis.** A clean Lyapunov theorem exists for an externally fixed covariance, but feedback covariance introduces an omitted derivative and destroys automatic gradient structure. **[Strong evidence]**

Expected impact was \(8/10\), difficulty \(8/10\), time 3–12 months, with estimated success \(0.65\) for a conditional theorem and \(0.25\) for a noncircular FIN-specific potential.

## 2. Background

For self-adjoint \(K\), a common matrix learning law is

\[
\dot K=-\Pi\bigl(V'(K)-C\bigr),
\]

where \(\Pi\) projects onto an admissible tangent space, \(V\) is a scalar spectral potential, and \(C\) is a target covariance. If \(C\) depends on \(K\), the system is a feedback law rather than the gradient of the fixed-target objective unless an integrability condition holds.

## 3. Mathematical formulation

For fixed \(C\), define

\[
\mathcal F_C(K)=\operatorname{Tr}V(K)-\operatorname{Tr}(KC).
\]

For feedback \(C=C(K)\), define the one-form

\[
\alpha_K(H)=\operatorname{Tr}\!\bigl(C(K)H\bigr).
\]

The central questions are whether \(\alpha\) is exact, whether a Lyapunov function is determined before seeing the target kernel, and whether odd-sector instability can break reflection symmetry.

## 4. Research methodology

The program:

1. differentiated both fixed-target and feedback objectives;
2. derived the integrability condition for the covariance one-form;
3. constructed a polynomial interpolation counterexample to potential uniqueness;
4. decomposed the Jacobian into even and odd reflection sectors;
5. compared unitary and dissipative trajectories generated from the same spatial operator;
6. tested the standard pitchfork normal form as a conditional selector mechanism.

## 5. Detailed derivations

For fixed \(C\),

\[
\nabla\mathcal F_C(K)=V'(K)-C,
\]

and projected gradient flow gives

\[
\frac{d\mathcal F_C}{dt}
=-\|\Pi\nabla\mathcal F_C\|_F^2\le0.
\]

This is a valid Lyapunov theorem under a fixed covariance and orthogonal tangent projection. **[Proven]**

If instead

\[
\widetilde{\mathcal F}(K)
=\operatorname{Tr}V(K)-\operatorname{Tr}(K\,C(K)),
\]

then for a variation \(H\),

\[
D\widetilde{\mathcal F}_K[H]
=\operatorname{Tr}\!\left[(V'(K)-C(K))H\right]
-\operatorname{Tr}\!\left[K\,DC(K)[H]\right].
\]

Hence

\[
\nabla\widetilde{\mathcal F}
=V'(K)-C(K)-DC(K)^*[K].
\]

The frequently written feedback law omits the final term and is not the gradient of this naive objective unless it vanishes. **[Proven]**

More generally, a potential for \(C(K)\) exists locally only if its Jacobian is symmetric in the Frobenius pairing,

\[
\langle DC(K)[H_1],H_2\rangle
=\langle H_1,DC(K)[H_2]\rangle.
\]

**[Proven]**

## 6. Proof attempts

Suppose a desired fixed point \(K_*\) and commuting \(C_*\) are already known. On the finite spectrum of \(K_*\), polynomial interpolation can always choose \(V'\) so that

\[
V'(K_*)=C_*.
\]

Therefore existence of a fitted potential is not evidence that the potential predicts \(K_*\). **[Proven]** Predictivity requires \(V\) to be specified independently and validated on held-out operators or sizes.

A reflection-equivariant law can undergo an odd-sector pitchfork only if its Jacobian has an odd eigenvalue crossing zero transversely and a stabilizing nonlinear coefficient. Center-manifold theory then yields

\[
\dot m=a(\mu-\mu_c)m-bm^3+O(m^5).
\]

The generic conditional theorem is valid, but the required FIN Jacobian crossing has not been derived. **[Speculative]**

## 7. Counterexamples

If the initial condition is exactly reflection even and the vector field is equivariant, uniqueness of ordinary differential equations keeps the trajectory in the even fixed subspace. It cannot generate orientation spontaneously at finite deterministic size. **[Proven]**

At finite equilibrium and zero field, symmetry pairs \(+m\) and \(-m\), so \(\langle m\rangle=0\). Noise selects random realized branches, not a canonical mathematical sign. **[Proven]**

A potential reverse-engineered after observing \(K_s\) is circular and has no out-of-sample content. **[Refuted]**

## 8. Numerical investigations

For the normal form \(\dot m=m-m^3+h\), exact symmetric initialization remains at \(m=0\) when \(h=0\); initial values \(\pm10^{-9}\) converge to \(\pm1\), and fields \(\pm10^{-4}\) select corresponding branches. **[Strong evidence]**

The same graph operator also permits unitary and diffusive flows. In the executed comparison, the optimized position distributions were strongly distinguishable, yet both were compatible with the same frozen spectrum. This confirms that the update/evolution law is independent data, a question developed fully in Programs 11–20. **[Strong evidence]**

## 9. Literature comparison

Matrix learning flows resemble Oja-type adaptive dynamics; the original principal-component rule is [Oja 1982](https://doi.org/10.1007/BF00275687). Equivariant bifurcation theory rigorously classifies branch orbits, but does not choose one without a state or perturbation. The FIN-specific unresolved point is provenance of \(V\), \(C(K)\), and the odd tangent source.

## 10. Final conclusions

The fixed-covariance law has a valid variational interpretation. **[Proven]**

The feedback law is not automatically the gradient of the advertised objective; an adjoint derivative or an independently exact covariance one-form is required. **[Proven]**

Neither a fitted potential nor a symmetric pitchfork uniquely selects the strict kernel or orientation. **[Proven]**

## 11. Confidence assessment

- Fixed-target Lyapunov theorem: **[Proven]**
- Automatic feedback-gradient interpretation: **[Refuted]**
- Post-hoc potential as evidence of prediction: **[Refuted]**
- FIN-specific odd bifurcation: **[Speculative]**
- Probability of an independently sourced FIN potential: \(0.25\)
- Probability of a useful conditional adaptive model: \(0.65\)
- Impact: \(8/10\)
- Difficulty: \(8/10\)
- Time: 3–12 months
- Dependencies: Programs 2, 8, 9, and a fully typed covariance law

---

# Program 6 — Information thermodynamics and the open-system bridge

## 1. Motivation

The repository uses information language, but Shannon entropy is dimensionless while thermodynamic entropy, heat, and work require a physical model. This program asks for the smallest rigorous information-to-physics bridge and tests whether it can avoid arbitrary constants.

**Research hypothesis.** Relative entropy becomes thermodynamic only after adding a Hamiltonian, temperature, bath, and allowed process; no unique energetic bridge follows from a probability distribution alone. **[Strong evidence]**

Estimated success probability for a conditional bridge was \(0.90\), for a FIN-specific prediction \(0.10\), impact \(7/10\), difficulty \(6/10\), and time 2–6 months.

## 2. Background

Shannon entropy

\[
H(p)=-\sum_i p_i\log p_i
\]

measures uncertainty in nats. Thermodynamic entropy is \(S_{\rm th}=k_BH\) only for a declared physical ensemble and convention. A Gibbs model adds energies \(E_i\) and inverse temperature \(\beta\),

\[
p_i=Z^{-1}e^{-\beta E_i}.
\]

Landauer's principle concerns a process coupling a memory to a thermal reservoir, not an algebraic property of \(H(p)\).

## 3. Mathematical formulation

A minimum thermodynamic completion is

\[
(\mathcal H,H,\rho,\mathcal R,H_R,\beta,U,\mathcal O),
\]

consisting of system Hilbert space and Hamiltonian, state, reservoir, reservoir Hamiltonian and temperature, an allowed global process, and observables.

The identification question is whether \(p\) or \(K\) uniquely determines any member of this tuple.

## 4. Research methodology

The program used:

1. inverse Gibbs construction to test uniqueness;
2. scale-orbit invariance under \(H\mapsto cH,\beta\mapsto\beta/c\);
3. an exact two-level erasure-by-SWAP calculation;
4. relative-entropy decomposition of heat;
5. comparison with Jaynes inference, KMS equilibrium, quantum resource theories, and modular entropy;
6. open-system counterexamples.

## 5. Detailed derivations

**Theorem 6 — Gibbs nonidentifiability.** Every strictly positive finite probability vector is a Gibbs distribution for infinitely many energy scales. **[Proven]**

**Proof.** For any \(\theta>0\) and constant \(C\), set

\[
E_i=-\theta\log p_i+C,\qquad \beta=\theta^{-1}.
\]

Then \(e^{-\beta E_i}=e^{-\beta C}p_i\), and normalization recovers \(p_i\). Varying \(\theta\) changes all energy gaps. \(\square\)

Thus neither temperature nor energy is determined by Shannon data.

For an initially uncorrelated system and Gibbs reservoir, global unitary evolution gives the exact finite Landauer identity

\[
\beta Q
=\Delta S
+I(S':R')
+D(\rho_R'\|\rho_R),
\]

where \(\Delta S=S(\rho_S)-S(\rho_S')\) is the system entropy reduction and \(Q\) is heat deposited in the reservoir. **[Proven]** The inequality \(\beta Q\ge\Delta S\) follows from nonnegativity of mutual and relative entropy.

## 6. Proof attempts

Jaynes' maximum-entropy principle derives a Gibbs family after the measured constraints and their physical observables are specified. It does not determine which operator is energy or the numerical energy unit. **[Proven]**

KMS equilibrium similarly requires a \(C^*\)-dynamical system \((\mathcal A,\alpha_t)\); entropy alone supplies neither \(\alpha_t\) nor calibrated time. **[Proven]**

Modular Hamiltonians \(-\log\rho\) are dimensionless. Calling them energy requires a temperature or time normalization. **[Proven]**

Thermodynamic resource theory can make state conversion operationally precise, but its free states and allowed operations are additional axioms. **[Strong evidence]**

## 7. Counterexamples

Under

\[
H\mapsto cH,\qquad \beta\mapsto\beta/c,
\]

the Gibbs state is unchanged. Every equilibrium probability and Shannon entropy is therefore compatible with a continuous orbit of energy/temperature calibrations. **[Proven]**

A pure closed state evolves unitarily with constant von Neumann entropy while the Shannon entropy of one chosen measurement basis can vary. Measurement uncertainty is not automatically thermodynamic entropy production. **[Proven]**

The statement “one bit has energy \(k_BT\ln2\)” is incomplete: the bound concerns erasure at a specified temperature and ideal limit, not possession of a bit in isolation. **[Refuted]**

## 8. Numerical investigations

A two-level system/reservoir SWAP at \(\beta=1.3\) yielded

\[
\Delta S=0.1737242,\qquad
Q=0.2858350,\qquad
D(\rho_R'\|\rho_R)=0.1978613,
\]

with zero final mutual information. The residual in

\[
\beta Q-\Delta S-D
\]

was \(8.33\times10^{-17}\). **[Proven]**

The Gibbs scale-orbit test gave zero numerical difference after simultaneous \(H\) and \(\beta\) rescaling. **[Proven]**

These calculations validate the conditional thermodynamic identity while demonstrating that its physical scale came from supplied \(\beta\) and energy gaps.

## 9. Literature comparison

The information measure is due to [Shannon](https://doi.org/10.1002/j.1538-7305.1948.tb01338.x); the maximum-entropy statistical construction to [Jaynes](https://doi.org/10.1103/PhysRev.106.620). A rigorous finite-size Landauer equality and corrections are given by [Reeb and Wolf](https://arxiv.org/abs/1306.4352). Quantum thermodynamic resource theories, such as [Brandão et al.](https://arxiv.org/abs/1305.5278), begin from an explicit Hamiltonian and bath. FIN currently offers an information representation, not an equivalent thermodynamic model.

## 10. Final conclusions

Shannon information can participate in physics only through an operational thermodynamic tuple. **[Proven]**

No arbitrary constant is needed inside the dimensionless relative-entropy identity, but converting it to joules requires a supplied energy/temperature calibration. **[Proven]**

The best bridge is conditional and falsifiable; it is not an internal emergence of temperature or energy. **[Strong evidence]**

## 11. Confidence assessment

- Gibbs nonidentifiability theorem: **[Proven]**
- Finite Landauer identity under stated assumptions: **[Proven]**
- Shannon entropy as physical entropy without a model: **[Refuted]**
- Coherent calibrated FIN thermodynamics: **[Moderate evidence]**
- Unique FIN temperature or energy scale: **[Refuted]**
- Probability of a useful conditional thermodynamic model: \(0.65\)
- Probability of a new FIN-specific thermodynamic prediction: \(0.10\)
- Impact: \(7/10\)
- Difficulty: \(6/10\)
- Time: 2–6 months
- Dependencies: Programs 3, 11–17, and an explicit bath

---

# Program 7 — Repaired noncommutative geometry, gauge fields, and fermions

## 1. Motivation

The repository invokes a spectral triple, while genuine noncommutative geometry can produce metric, gauge, and fermionic structures only when its axioms and representations are satisfied. This program tests the minimal tuple and designs the smallest honest repair.

**Research hypothesis.** The naive vertex triple fails standard real/even conditions; a doubled graph triple is feasible but does not uniquely generate nonabelian gauge or fermion content. **[Strong evidence]**

Estimated success \(0.70\) for a repaired graph triple, below \(0.05\) for unique physical gauge content, impact \(8/10\), difficulty \(8/10\), time 4–12 months.

## 2. Background

A finite spectral triple \((\mathcal A,H,D)\) needs a represented algebra, Hilbert space, and self-adjoint Dirac operator. A real/even triple additionally needs a real structure \(J\), grading \(\gamma\), KO signs, first-order condition, and usually orientability and duality constraints.

The tuple

\[
(\mathbb C^{12},\mathbb C^{12},K_s)
\]

is a finite spectral triple in the weak sense: commutators are bounded and the resolvent is compact. **[Proven]** Those properties are automatic in finite dimension and do not establish the stronger geometry.

## 3. Mathematical formulation

Let \(\mathcal A=\mathbb C^{12}\) act diagonally on \(H=\mathbb C^{12}\), \(D=K_s\), and let \(J\) be componentwise conjugation. Test

\[
[[D,a],Jb^*J^{-1}]=0
\]

and seek \(\gamma^2=1\) with \([\gamma,a]=0\) and \(\{D,\gamma\}=0\).

Then compare with the doubled incidence candidate

\[
H_{\rm graph}=\ell^2(V)\oplus\ell^2(E),\qquad
D_{\rm graph}=
\begin{pmatrix}0&d^*\\ d&0\end{pmatrix},
\qquad
\gamma=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

## 4. Research methodology

The program used:

1. symbolic component formulas for the first-order condition;
2. exhaustive point-projection residuals;
3. a graph-coloring contradiction for the grading;
4. unitary-group analysis of \(\mathbb C^{12}\);
5. a constructive doubled vertex/edge proposal;
6. a provenance test separating valid added axioms from consequences of \(K_s\).

## 5. Detailed derivations

For diagonal \(a,b\),

\[
[[D,a],Jb^*J^{-1}]_{ij}
=D_{ij}(a_j-a_i)(b_j-b_i)
\]

up to the immaterial conjugation convention. Choosing point projections exposes a nonzero entry whenever \(D_{ij}\ne0\). **[Proven]**

The executed maximum Frobenius residual over point projections was \(1.037269\), so the naive real first-order condition fails. **[Refuted]**

For an even grading, multiplicity-one diagonal representation forces

\[
\gamma=\operatorname{diag}(\gamma_i),\qquad \gamma_i=\pm1.
\]

Then

\[
\{D,\gamma\}_{ij}=D_{ij}(\gamma_i+\gamma_j).
\]

Every off-diagonal strict entry is nonzero, so every pair would need opposite signs. Three vertices already give a contradiction. No compatible even grading exists on this representation. **[Proven]**

## 6. Proof attempts

The doubled incidence construction has a valid formal grading and a Laplace-type square. It creates room for a nontrivial real structure and differential calculus. **[Strong evidence]** However, first-order, orientability, distance, and role of the long-range kernel must be checked on the completed representation.

A precise research target is whether a stable filter \(f\) exists with

\[
K_s=f(D_{\rm graph}^2)
\]

under a refinement compatible with Program 4. Even if it exists, \(f\) and the graph calculus need an independent selection principle. **[Moderate evidence]**

## 7. Counterexamples

The unitary group of the current algebra is

\[
U(\mathbb C^{12})=U(1)^{12},
\]

which is abelian. It cannot contain derived \(SU(2)\) or \(SU(3)\) factors. **[Proven]**

Nonabelian gauge fields require matrix-algebra summands and a bimodule. Chiral fermions require left/right representations, KO dimension, grading, and internal Dirac/Yukawa data. The twelve eigenvalues do not determine those choices. **[Proven]**

Calling \(K_s\) a Dirac operator does not repair failed axioms. **[Refuted]**

## 8. Numerical investigations

The first-order residual scan and grading contradiction were exact finite tests. A diagonal choice of \(D\) makes the first-order residual zero, but then all commutators with \(\mathcal A\) vanish and the induced metric is trivial or infinite between distinct points. **[Proven]**

This establishes a useful repair constraint: removing off-diagonal structure solves one axiom only by destroying the geometry.

## 9. Literature comparison

Finite graph spectral triples are legitimate but depend on the selected calculus; see [Requardt](https://arxiv.org/abs/hep-th/9708010). Finite real triples can generate Yang–Mills/Higgs sectors after algebra and bimodule data are supplied; see [Krajewski](https://arxiv.org/abs/hep-th/9701081). The [spectral action](https://arxiv.org/abs/hep-th/9606001) is a conditional construction, not a proof that the current FIN tuple satisfies real/even axioms.

## 10. Final conclusions

The weak finite spectral triple is mathematically valid but physically underpowered. **[Proven]**

The naive real/even FIN triple is falsified. **[Refuted]**

A repaired doubled graph geometry is a tractable extension, but any nonabelian gauge group or fermion sector will be additional structured data unless a uniqueness theorem is found. **[Strong evidence]**

## 11. Confidence assessment

- Weak finite spectral-triple status: **[Proven]**
- Naive first-order and even-grading claims: **[Refuted]**
- Feasibility of a repaired graph triple: **[Moderate evidence]**
- Unique Standard-Model-like content from \(K_s\): **[Speculative]**
- Probability of a rigorous repaired triple: \(0.70\)
- Probability of uniquely forced nonabelian gauge/fermion content: \(0.05\)
- Impact: \(8/10\)
- Difficulty: \(8/10\)
- Time: 4–12 months
- Dependencies: Programs 4, 8, and a new algebra/bimodule

---

# Program 8 — Variational stability and the Osterwalder–Schrader gate

## 1. Motivation

A genuine field theory requires a bounded action, a probability measure or quantum phase, equations of motion, and eventually a route from Euclidean correlators to unitary observables. This program tests whether the strict matrix can already serve as a quadratic action and what additional structure is minimally necessary.

**Research hypothesis.** The unshifted strict quadratic form is unstable; a positive completion is easy but nonunique, and reflection positivity is undefined until a temporal theory is supplied. **[Strong evidence]**

Estimated success \(0.80\) for a stable finite model, \(0.15\) for an OS-positive continuum theory, impact \(9/10\), difficulty \(9/10\), time 6–18 months.

## 2. Background

For a real field \(\phi\in\mathbb R^{12}\), the Gaussian Euclidean action

\[
S_2[\phi]=\frac12\phi^TA\phi
\]

defines a normalizable measure only when \(A>0\) after zero-mode treatment. A field theory additionally needs interactions, a temporal lattice or manifold, reflection, locality/control, and a continuum limit.

The Osterwalder–Schrader reconstruction theorem converts Euclidean correlation functions to a Hilbert-space quantum theory only after regularity, Euclidean covariance, symmetry, clustering, and reflection positivity.

## 3. Mathematical formulation

Test the family

\[
A_\mu=K_s+\mu I,
\qquad
S_{\mu,\lambda}[\phi]
=\frac12\phi^TA_\mu\phi
+\lambda\sum_i\phi_i^4.
\]

A graph-Laplacian alternative is

\[
A_{\rm M}=sI-K_s,
\qquad s=\sum_j(K_s)_{ij}.
\]

For stress–energy, one would need a metric-dependent family \(S[\phi;g]\) and define a response by functional differentiation with respect to \(g\).

## 4. Research methodology

The program:

1. computed the full inertia and stability threshold;
2. measured conditioning near the positive boundary;
3. compared scalar shifts with the Perron/Laplacian completion;
4. derived the dependence of observables on a metric family rather than one matrix value;
5. formulated an OS-positive product-theory target;
6. searched for counterexamples to uniqueness and direct physical interpretation.

## 5. Detailed derivations

The strict matrix has

\[
\lambda_{\min}=-0.68187476238
\]

and seven negative eigenvalues. Therefore \(S_2=\phi^TK_s\phi/2\) is unbounded below along seven independent directions and its Gaussian integral diverges. **[Proven]**

The shifted action is positive definite exactly when

\[
\mu>0.68187476238.
\]

**[Proven]** A quartic interaction with \(\lambda>0\) can bound the total finite-dimensional integral even when the quadratic term is indefinite, but the chosen interaction and phase structure are additional data.

The Perron completion \(A_{\rm M}=sI-K_s\) is positive semidefinite because \(K_s\) is symmetric and entrywise positive off diagonal with row sum \(s\). It has the Dirichlet form

\[
\phi^TA_{\rm M}\phi
=\frac12\sum_{ij}K_{ij}(\phi_i-\phi_j)^2.
\]

**[Proven]** This is a mathematically distinguished graph-Laplacian completion on \(C_{12}\), but its role as physical action is not selected by the original interpretation.

## 6. Proof attempts

A conditional Euclidean target is

\[
S_E[\phi]
=\frac12\langle\phi,(-\Delta_\tau+A_{\rm M}+m^2)\phi\rangle
+\lambda\sum_{(\tau,i)}\phi_{\tau i}^4
\]

on a reflected temporal lattice. For the free part with \(A_{\rm M}+m^2>0\), standard transfer-matrix methods make reflection positivity plausible; the exact lattice measure and interactions still require proof. **[Moderate evidence]**

No theorem derives the temporal lattice, reflection plane, field representation, \(\lambda\), or continuum scaling from \(K_s\). **[Proven]**

A spectral action does not solve selection because both its cutoff and test function are explicit inputs.

## 7. Counterexamples

A single value \(K(g_0)\) cannot determine stress–energy. The families

\[
K_1(g)=K_s,\qquad K_2(g)=K_s+(g-g_0)B
\]

agree at \(g_0\) but have arbitrary different derivatives \(0\) and \(B\). **[Proven]** Metric response requires a declared family, not only a matrix.

Adding any \(\mu\) above threshold gives a stable Gaussian, so stability does not uniquely choose a mass or action. **[Proven]**

Matrix inertia is not a Wick-rotated spacetime signature, and shifting negative modes is not a derivation of Lorentzian time. **[Refuted]**

## 8. Numerical investigations

There are seven negative and five positive strict modes. Setting

\[
\mu=|\lambda_{\min}|+\varepsilon
\]

gave condition numbers

| \(\varepsilon\) | \(\kappa(A_\mu)\) |
|---:|---:|
| 0.001 | 2343.18 |
| 0.01 | 235.22 |
| 0.1 | 24.42 |
| 1 | 3.34 |

The near-threshold theory is extremely sensitive to perturbations. **[Strong evidence]** A physically credible fit cannot place the shift near threshold without uncertainty propagation.

## 9. Literature comparison

The OS axioms and reconstruction are due to [Osterwalder and Schrader](https://doi.org/10.1007/BF01645738). Their theorem requires a hierarchy of Euclidean correlators with temporal reflection, not one spatial kernel. The spectral action of [Chamseddine and Connes](https://arxiv.org/abs/hep-th/9606001) supplies a powerful conditional variational framework but explicitly introduces its cutoff and algebraic geometry.

## 10. Final conclusions

The unshifted strict kernel is not a stable Gaussian action. **[Refuted]**

Stable finite actions are easy to construct after a shift or quartic term, but they are nonunique. **[Proven]**

Reflection positivity, stress–energy, and field-theoretic observables cannot yet be evaluated because the required temporal, metric, state, and refinement structures are absent. **[Proven]**

## 11. Confidence assessment

- Instability and shift threshold: **[Proven]**
- Perron graph-Laplacian completion: **[Proven]**
- Unique variational action selected by the kernel: **[Refuted]**
- Conditional OS-positive product theory: **[Moderate evidence]**
- FIN-native OS-positive continuum: **[Speculative]**
- Probability of a stable finite model: \(0.80\)
- Probability of a genuine OS-positive continuum: \(0.15\)
- Impact: \(9/10\)
- Difficulty: \(9/10\)
- Time: 6–18 months
- Dependencies: Programs 3, 4, 7, and 11–20

---

# Program 9 — Calibrated predictive fingerprint and Bayesian experiment

## 1. Motivation

The decisive transition from mathematics to physics is a risky, held-out prediction linked to a preparation and measurement protocol. This program seeks a dimensionless signature that survives unit calibration and an experiment that can distinguish competing dynamics generated from the same operator.

**Research hypothesis.** A twelve-mode simulator can test FIN-specific spectral and dynamical fingerprints, but reproducing an engineered operator does not establish that FIN describes an external physical domain. **[Strong evidence]**

Estimated success \(0.70\) for an executable simulator test, \(0.05\) for relevance to fundamental physics without a domain map, impact \(10/10\), difficulty \(7/10\), time 3–12 months.

## 2. Background

An operational theory specifies preparations, transformations, and measurement probabilities. A bare matrix has no empirical likelihood until mapped to this structure.

To avoid arbitrary energy offsets and scales, spectral tests should use affine-invariant ratios. To avoid post-selection, calibration modes must be fixed first and the remaining modes treated as held out.

## 3. Mathematical formulation

Let sorted nondegenerate radial eigenvalues be \(\lambda_0,\ldots,\lambda_6\). Define an affine-normalized fingerprint

\[
x_j=\frac{\lambda_j-\lambda_{\min}}
{\lambda_{\max}-\lambda_{\min}}.
\]

For dynamics, compare simple hypotheses

\[
\mathsf H_U:\ q_x(t)=|\langle x|e^{-itA}|0\rangle|^2,
\qquad
\mathsf H_D:\ p_x(t)=\langle x|e^{-tA}|0\rangle,
\]

when \(A\) is a positive Markov/Laplacian generator. A likelihood-ratio experiment is meaningful after preparation, clock, and site measurement are calibrated.

## 4. Research methodology

The program:

1. computed affine-invariant strict and legacy spectral fingerprints;
2. fitted scale and offset using declared calibration information;
3. evaluated held-out normalized errors;
4. generated unitary and diffusive position distributions;
5. optimized a discrimination time under information criteria;
6. designed positive, negative, and apparatus-reflection controls;
7. separated simulator validation from claims about nature.

## 5. Detailed derivations

If measured energies obey

\[
E_j=a\lambda_j+b,
\qquad a>0,
\]

then all \(x_j\) are unchanged. **[Proven]** They are legitimate unit-free observables once the mode-to-mode correspondence is specified.

The best affine legacy-to-strict spectral fit had scale \(0.109565\) and offset \(0.542413\). Its held-out normalized RMSE was \(0.273796\), with maximum error \(0.354079\). **[Strong evidence]** Thus a mere two-parameter calibration does not make the two kernel spectra empirically interchangeable.

For independent observations generated under \(\mathsf H_U\), the expected log Bayes factor against \(\mathsf H_D\) is

\[
\mathbb E_U[\log B_{UD}]
=nD_{\rm KL}(q\|p).
\]

This converts an operator comparison into a sample-size forecast under simple, fully specified hypotheses. **[Proven]**

## 6. Proof attempts

The operational map can be stated as

\[
\mathcal M:
\{\text{laboratory preparations, controls, outcomes}\}
\longrightarrow
\{\rho,\Phi_t,M_y\}.
\]

No spectral theorem constructs \(\mathcal M\). **[Proven]** It must be supplied and falsified independently.

A candidate FIN-specific test requires at least two held-out spectral ratios or a full time trace not used for parameter fitting. One fitted frequency or one visual profile is never sufficient because scale and offset can absorb it. **[Proven]**

Full process tomography near \(t=0\) would identify the calibrated generator; population-only one-time data generally do not. This motivates the multitime programs below.

## 7. Counterexamples

Any programmable simulator can be engineered to implement a chosen matrix. Agreement then validates the engineering model, not the claim that the same matrix governs an unrelated natural system. **[Proven]**

The directional odd-current response is generic to even circulants once an odd control is added. Observing it alone does not uniquely validate FIN. **[Proven]**

Clock error, mode permutation, disorder, dephasing, and flexible nuisance parameters can imitate or erase a nominal fingerprint. A valid test must preregister their ranges and retain held-out observables. **[Strong evidence]**

## 8. Numerical investigations

The normalized adjacent-gap fingerprints were

\[
g_s=(0.01860,0.04228,0.10168,0.16410,0.35135,0.32197),
\]

\[
g_l=(0.37268,0.02033,0.08844,0.02948,0.10795,0.38110).
\]

Their incompatible order structure supplies a scale-free discriminator. **[Strong evidence]**

In the executed unitary/diffusive comparison, an optimized dimensionless time \(t=6.77\) gave

\[
D_{\rm JS}=0.36293,\qquad
D_{\rm KL}(q\|p)=1.39199.
\]

Under ideal independent sampling, eight observations give expected log Bayes factor exceeding ten. The unitary and diffusive position entropies were \(1.04034\) and \(2.45800\), and return probabilities \(0.60433\) and \(0.11069\). **[Strong evidence]**

These figures are design diagnostics, not laboratory evidence; Programs 11–20 refine the protocol and nuisance analysis.

## 9. Literature comparison

Operational probabilistic theories make preparation and measurement primitives explicit; the purification-based quantum framework of [Chiribella, D'Ariano and Perinotti](https://arxiv.org/abs/0908.1583) illustrates how extra operational axioms select quantum structure. Bayesian Hamiltonian learning and experiment design are developed by [Hincks et al.](https://arxiv.org/abs/1806.02427). Operational state-discrimination limits are reviewed in [Guff et al.](https://arxiv.org/abs/1906.09737). FIN has not yet supplied an empirical domain to which these methods can be applied, but they define the correct gate.

## 10. Final conclusions

A calibrated twelve-mode spectral/dynamical experiment is feasible in principle. **[Moderate evidence]**

The strict and legacy fingerprints remain distinguishable after affine calibration. **[Strong evidence]**

Only preregistered held-out predictions in an independently motivated physical system can change FIN from a mathematical model into physics. **[Proven]**

## 11. Confidence assessment

- Affine-invariant fingerprint construction: **[Proven]**
- Reported numerical separation: **[Strong evidence]**
- Feasible engineered twelve-mode test: **[Moderate evidence]**
- Emulator agreement as evidence for fundamental physics by itself: **[Refuted]**
- Probability of executable simulator experiment: \(0.70\)
- Probability current FIN already identifies a natural empirical domain: \(0.05\)
- Impact: \(10/10\)
- Difficulty: \(7/10\)
- Time: 3–12 months
- Dependencies: Programs 2, 3, 4, and 11–20

---

# Program 10 — Machine-assisted invariant discovery and proof certificates

## 1. Motivation

The repository contains thousands of generated artifacts and repeated searches for selectors, sources, and closure. This program asks where exhaustive computation and formal proof can genuinely close classes of possibilities instead of creating more candidate prose.

**Research hypothesis.** Symmetry and dimensional no-go results are well suited to machine certification; unrestricted automated search cannot manufacture a physically meaningful source from invariant inputs. **[Strong evidence]**

Estimated success \(0.90\) for formal no-go certificates, \(0.05\) for discovery of a new strict physical source from current data, impact \(7/10\), difficulty \(6/10\), time 2–6 months.

## 2. Background

Finite group actions, representation weights, polynomial identities, spectral residuals, and SAT-size combinatorics can be exhaustively checked. A proof assistant can then certify the general theorem separately from floating-point evidence.

The central danger is search-space leakage: a system can appear to find a selector by using array order, file names, target values, or conventions that are not invariant mathematical inputs.

## 3. Mathematical formulation

Let \(G\) act on an input representation \(V\) and target \(W\). Candidate source maps lie in

\[
\operatorname{Hom}_G(V,W).
\]

If \(V\) is trivial and \(W\) is reflection odd or dimension weight one, this Hom-space is zero.

For combinatorial profiles \(x\in\{0,1\}^{12}\), the task is to enumerate orbits under translations and dihedral symmetry and distinguish selection of an orbit from selection of a representative.

## 4. Research methodology

The program used:

1. exact representation decompositions;
2. exhaustive enumeration of invariant Boolean predicates and score functions;
3. dimensional-weight monomial search;
4. orbit canonicalization for all nonzero binary profiles;
5. numerical residual scripts with fixed seeds and tolerances;
6. a formal-certificate architecture separating theorem, finite witness, provenance, and physical interpretation.

## 5. Detailed derivations

**Theorem 7 — invariant-to-odd exclusion.** If the input representation is \(G\)-trivial and target \(W\) has no fixed vector, then

\[
\operatorname{Hom}_G(V,W)=0.
\]

**[Proven]** This simultaneously covers the absolute orientation torsor at the set level and odd linear sources at the representation level.

For dimensional characters, products add weights. The subring generated by weight-zero quantities has weight zero, so it contains no nonzero weight-one element. **[Proven]**

An invariant algorithm on a transitive orbit can return an orbit label but cannot select a representative equivariantly. **[Proven]**

These are short theorems suitable for proof-assistant formalization.

## 6. Proof attempts

A valid automated discovery pipeline should require each candidate to export:

1. a typed formula;
2. its transformation law;
3. provenance from admitted inputs;
4. exact or interval-certified nonzero value;
5. a theorem coupling it to the target;
6. a held-out falsification test.

Without all six, a numerical pattern is not a new source. **[Strong evidence]**

SMT/SAT can exhaust finite selector predicates; invariant theory can classify polynomial covariants; interval arithmetic can certify spectral signs; Lean or Isabelle can encode the no-fixed-point and weight proofs. No fundamental obstacle was found for these bounded tasks. **[Strong evidence]**

## 7. Counterexamples

A neural network trained on labelled “clockwise” arrays learns the labelling convention unless data augmentation and equivariance remove it. High validation accuracy on the same convention is not selector discovery. **[Refuted]**

Symbolic regression can reconstruct the degree-six finite kernel bridge, but Program 1 shows that exact fit may extrapolate catastrophically. Formula simplicity on seven points is not a physical law. **[Proven]**

Automated theorem generation cannot prove a false existential statement if its checker is sound; it should be used to certify no-go boundaries, not to promote uncited interpretations. **[Strong evidence]**

## 8. Numerical investigations

All \(2^{12}-1=4095\) nonzero binary profiles were enumerated. They form 351 translation orbits and 223 dihedral orbits. **[Proven]** Dihedral quotienting removes orientation distinctions but never selects an absolute representative.

The even-to-odd equivariant linear-map dimension was zero. The search of \(15\,624\) weight-zero monomials produced zero weight-one candidates. **[Proven]**

The accompanying reproducibility suite also certifies polynomial residuals, selector symmetry, stability thresholds, continuum scaling, and fingerprint mismatches.

## 9. Literature comparison

Computer algebra, invariant theory, SAT, and interactive theorem proving are established methods; their value here is disciplined scope. Free-probability or random-matrix searches are not currently justified: the strict and legacy matrices commute, and there is no ensemble or large-\(N\) law. Modern universality results such as [Tao and Vu](https://arxiv.org/abs/1202.0068) presuppose precisely those absent hypotheses.

## 10. Final conclusions

Machine assistance can close finite symmetry and source classes with auditable certificates. **[Strong evidence]**

It cannot convert a chart convention, target leakage, or dimensionless receiver into a physical source. **[Proven]**

The immediate deliverable should be a small formal library of no-go and provenance theorems, linked to deterministic witnesses, rather than another unbounded scan. **[Strong evidence]**

## 11. Confidence assessment

- Finite orbit and weight calculations: **[Proven]**
- Formal certification feasibility: **[Strong evidence]**
- Automated discovery of an absolute selector from current even data: **[Refuted]**
- Discovery of a genuinely new typed source after new input: **[Speculative]**
- Probability of formalizing the main no-go results: \(0.90\)
- Probability current artifacts hide an admissible new physical source: \(0.05\)
- Impact: \(7/10\)
- Difficulty: \(6/10\)
- Time: 2–6 months
- Dependencies: all programs supply theorem statements and witnesses

---

# Part III — Maximum-creativity cross-field search

The following directions were tested as generators rather than merely analogies.

| Field | Result of admission test | Status |
|---|---|---|
| spectral geometry | exactly describes the finite functional core; needs a refinement for dimension/continuum | **[Proven]** |
| operator algebras | current algebra is finite and commutative; richer dynamics requires added algebra/state | **[Proven]** |
| modular theory | trivial on \(\mathbb C^{12}\); state-dependent and uncalibrated on \(M_{12}\) | **[Proven]** |
| index theory/spectral flow | no selector without an odd path, grading, or chiral source | **[Proven]** |
| equivariant topology | proves the no-section obstruction; relational torsor quotient is constructive | **[Proven]** |
| \(K\)-theory/topological phases | a single finite filter lacks a gapped parameter family and occupied-state prescription | **[Strong evidence]** |
| category/sheaf/homotopy/derived methods | useful language for naturality and torsors; functoriality preserves the blocking automorphism | **[Proven]** |
| information geometry | Fisher metrics quantify parameter distinguishability but need a likelihood family | **[Proven]** |
| optimal transport | can describe dissipative gradient flows after a cost/metric is chosen; it does not select that choice | **[Strong evidence]** |
| renormalization | potentially decisive only with a refinement, blocking map, and reference scale | **[Proven]** |
| causal sets | require an acyclic order; an oriented cycle is not causal | **[Proven]** |
| tensor networks | strict operator has low Schmidt rank under several factorizations, but factorization is not selected | **[Strong evidence]** |
| random matrices | no ensemble or large-\(N\) law currently exists | **[Proven]** |
| free probability | impossible for the present commuting non-scalar pair in a faithful trace | **[Refuted]** |
| neural operators | can approximate update maps, but training data and loss replace the missing physical law | **[Proven]** |
| complexity theory | can rank inference cost after an input/output task exists; it cannot create units or observables | **[Strong evidence]** |
| Bayesian physics | gives the strongest immediate route to model discrimination and honest nuisance treatment | **[Strong evidence]** |
| machine-discovered mathematics | valuable for invariant classification and certificates, dangerous under target leakage | **[Strong evidence]** |
| holomorphic semigroups | unexpectedly unifies wave/unitary and diffusion as different complex-time restrictions | **[Proven]** |
| operational reference frames | replaces an impossible absolute selector with complete relative observables | **[Proven]** |

The tensor test deserves a quantitative note. Operator-Schmidt ranks of the strict matrix were \(2,3,4,3,3\) for sequential \(2\times6\), \(3\times4\), \(4\times3\), \(6\times2\), and CRT \(3\times4\) factorizations. **[Strong evidence]** This is useful compression but not evidence for a unique physical subsystem structure.

The most creative positive outcome is not a new exotic ontology. It is the conjunction of two standard but previously separated ideas:

\[
\text{one positive spectral generator}
\quad+\quad
\text{an operational temporal/reference choice}.
\]

It yields both the relative-frame theorem and the wave–diffusion ambiguity analyzed next. **[Strong evidence]**

---

# Part IV — Ten additional executed programs: one operator, two temporal worlds

## Observer paradox: precise statement

Let \(W=K_s\), let its constant row sum be

\[
s=1.660307278766099,
\]

and define the positive graph Laplacian

\[
A=sI-W.
\]

Because every off-diagonal QW2118 weight is positive,

\[
\langle f,Af\rangle
=\frac12\sum_{x,y}W_{xy}|f_x-f_y|^2\ge0.
\]

Thus the same \(A\) defines

\[
U_t=e^{-itA}
\quad\text{and}\quad
P_t=e^{-tA}.
\]

The first is a reversible unitary group on amplitudes; the second is an irreducible, uniform-reversible Markov semigroup on probabilities. **[Proven]**

The “observer paradox” is therefore not consciousness-dependent. It is the operational underdetermination of dynamics by a spatial spectral generator when the temporal ray, state cone, preparation, instrument, environment, record, and clock are unspecified. **[Proven]**

# Program 11 — Dual semigroup theorem for the strict FIN operator

## 1. Motivation

The user's central observation is exact for the QW2118 sign pattern: one positive operator supports both coherent and diffusive evolution. The first program determines whether this is merely formal or yields two fully valid dynamical laws.

**Research hypothesis.** Both laws are mathematically canonical once a temporal category is chosen, but neither category is selected by the operator. **[Strong evidence]**

Estimated success \(0.99\), impact \(10/10\), difficulty \(4/10\), time two to four weeks. It is the foundation for Programs 12–20.

## 2. Background

A self-adjoint \(A\) generates a unitary group \(e^{-itA}\). A graph Laplacian with nonpositive off-diagonal entries and zero row sums also generates a stochastic contraction semigroup \(e^{-tA}\).

The strict matrix \(W\) has positive weights at every nonzero cyclic distance and constant row sum \(s\). Its Perron shift therefore has both structures.

## 3. Mathematical formulation

Set

\[
Q=W-sI,\qquad A=-Q=sI-W.
\]

Test:

\[
Q_{ij}\ge0\ (i\ne j),\qquad Q\mathbf1=0,
\qquad Q=Q^T,
\]

and compare the state laws

\[
\psi(t)=U_t\psi_0,\qquad p(t)=P_tp_0.
\]

## 4. Research methodology

All six radial weights, the row sum, and the Fourier eigenvalues were recomputed. Positivity was proved through the Dirichlet form. Markov conservation and reversibility were checked algebraically. Unitary and Markov returns and entropies were then evaluated from one localized preparation.

## 5. Detailed derivations

Since \(Q\) has nonnegative off-diagonal entries and zero row sums, \(e^{tQ}=e^{-tA}\) is entrywise nonnegative and stochastic. Symmetry makes it doubly stochastic and reversible with respect to the uniform distribution. Irreducibility follows from positive connectivity. **[Proven]**

The Markov rates are

\[
0,\ 0.754121154^{(2)},\ 1.577049514^{(2)},\
1.961406862^{(2)},\ 2.199568849^{(2)},\
2.298606272^{(2)},\ 2.342182041.
\]

Self-adjointness gives \(U_t^*U_t=I\). Because \(A=sI-W\), evolution generated by \(A\) and \(W\) differs only by a global phase and sign convention at the level of unitary transition probabilities. **[Proven]**

## 6. Proof attempts

**Theorem 8 — trivial intersection.** A real matrix that is both nonnegative stochastic and unitary is a permutation matrix. Any continuous one-parameter family of such matrices starting at \(I\) is identically \(I\). **[Proven]**

**Proof.** Nonnegative orthogonal columns have disjoint support. Twelve nonzero columns on twelve coordinates therefore each have one entry, whose norm forces it to be one. Permutations are discrete, so a continuous path from \(I\) is constant. \(\square\)

Hence a nontrivial law cannot simultaneously treat the same entries as quantum amplitudes and classical transition probabilities. The state cone and probability rule are indispensable.

## 7. Counterexamples

The uniform state is fixed by both dynamics. An observer restricted to uniform preparation and population measurements sees identical constant data. **[Proven]** Identifiability is therefore preparation-relative.

The Markov interpretation is not automatically stable under naive fixed-formula refinement. The profile gives \(K(8)=-0.00179588\); at \(N=16\) this becomes an off-diagonal weight and violates the generator positivity condition. **[Proven]**

Thus exact Markov duality is established for QW2118 on \(C_{12}\), not for every proposed FIN extension.

## 8. Numerical investigations

From a localized vertex:

| \(t\) | unitary return | unitary position entropy | Markov return | Markov entropy |
|---:|---:|---:|---:|---:|
| 0.1 | 0.994634 | 0.040332 | 0.849359 | 0.716153 |
| 1 | 0.584353 | 1.293108 | 0.262829 | 2.225996 |
| 5 | 0.142549 | 2.259728 | 0.087250 | 2.484376 |

The Markov entropy tends to \(\log12=2.48490665\). The unitary position entropy is nonmonotone; on \(5\le t\le200\) it ranged approximately from \(0.8565\) to \(2.4802\). **[Strong evidence]**

## 9. Literature comparison

Continuous-time quantum walks generated by graph operators were formulated by [Farhi and Gutmann](https://arxiv.org/abs/quant-ph/9706062). Reversible continuous-time Markov chains use the same Laplacian class with real rather than imaginary time. The FIN-specific theorem is that the actual strict QW2118 weights satisfy both admission conditions after the Perron shift.

## 10. Final conclusions

The same strict spectral generator produces two legitimate but inequivalent temporal worlds. **[Proven]**

The spectrum fixes their mode shapes and rates/frequencies but does not select the state ontology, real versus imaginary temporal ray, or observation rule. **[Proven]**

The observer paradox is a precise model-underdetermination theorem. **[Proven]**

## 11. Confidence assessment

- Exact QW2118 unitary/Markov duality: **[Proven]**
- Physical equivalence of the two dynamics: **[Refuted]**
- Stable Markov interpretation under naive refinement: **[Refuted]**
- Probability the finite theorem survives independent recomputation: \(0.999\)
- Impact: \(10/10\)
- Difficulty: \(4/10\)
- Time: 2–4 weeks for formal publication
- Dependencies: none beyond the declared weights

---

# Program 12 — Holomorphic time, Wick continuation, and its obstruction

## 1. Motivation

If wave and diffusion are restrictions of one complex-time family, perhaps analyticity selects or reconstructs the physical branch. This program tests that tempting resolution.

**Research hypothesis.** Analytic continuation is exact in finite dimension but physically nonselective and exponentially unstable; Euclidean reconstruction needs reflection positivity and additional axioms. **[Strong evidence]**

Estimated success \(0.95\), impact \(9/10\), difficulty \(7/10\), time 2–6 months.

## 2. Background

For bounded \(A\ge0\),

\[
F(z)=e^{-zA}
\]

is entire. The positive real axis gives a contraction semigroup and the imaginary axis a unitary group. In field theory, continuation from Euclidean to Lorentzian correlators is justified only for specially constrained data.

## 3. Mathematical formulation

Study

\[
P_\tau=F(\tau),\qquad
U_t=F(it),
\]

the inverse problem \(P_\tau\mapsto A\mapsto U_t\), its conditioning, and the additional conditions required to interpret \(P_\tau\) as Euclidean quantum data.

## 4. Research methodology

The program derived the entire functional calculus, diagonalized the inverse continuation, computed noise amplification from the largest QW2118 rate, and compared a spatial heat kernel with the Osterwalder–Schrader hypotheses.

## 5. Detailed derivations

For finite \(A\),

\[
F(z)=\sum_{n=0}^\infty\frac{(-zA)^n}{n!}
\]

converges everywhere, so

\[
P_\tau=U_{-i\tau}
\]

as a matrix identity. **[Proven]**

If \(A=\sum_ma_mP_m\), backward heat continuation multiplies errors in the \(m\)-th mode by \(e^{\tau a_m}\). Its worst condition number is

\[
\kappa(\tau)=e^{\tau a_{\max}},
\qquad a_{\max}=2.342182041.
\]

**[Proven]**

## 6. Proof attempts

The OS route would require temporally indexed correlation functions, Euclidean covariance, symmetry, regularity, clustering, and reflection positivity. A single spatial \(e^{-\tau A}\) does not prove any such hierarchy. **[Proven]**

Even exact knowledge of \(F\) does not choose which ray \(z=e^{i\theta}t\) is physical. Contractivity selects \(|\theta|<\pi/2\); norm preservation selects the imaginary axis; those are extra dynamical axioms. **[Proven]**

## 7. Counterexamples

Born probabilities

\[
|\langle i|U_t|j\rangle|^2
\]

are quadratic in amplitudes and are not obtained by analytic continuation of Markov probabilities \(\langle i|P_t|j\rangle\). **[Proven]** Common analytic ancestry is not observational equivalence.

Backward diffusion is well defined for exact finite data but unstable under arbitrarily small noise. Therefore “Wick rotate the observed diffusion” is not a robust empirical reconstruction procedure. **[Refuted]**

A positive heat kernel alone need not satisfy reflection positivity for a proposed time reflection and observable algebra. **[Proven]**

## 8. Numerical investigations

Noise amplification was

| \(\tau\) | \(\kappa(\tau)\) |
|---:|---:|
| 1 | \(1.04\times10^1\) |
| 2 | \(1.08\times10^2\) |
| 5 | \(1.22\times10^5\) |
| 10 | \(1.49\times10^{10}\) |
| 20 | \(2.21\times10^{20}\) |

The finite analytic continuation is numerically unusable at moderate long times without strong regularization. **[Strong evidence]**

## 9. Literature comparison

Euclidean reconstruction under reflection positivity is the content of [Osterwalder and Schrader](https://doi.org/10.1007/BF01645738). Holomorphic semigroups and unitary boundary values are standard operator theory. The new FIN insight is their exact application to the QW2118 dual generator and the quantified inverse instability.

## 10. Final conclusions

Wave/unitary and diffusion are exact complex-time shadows of one positive generator. **[Proven]**

Analyticity does not choose the physical temporal axis, probability rule, or clock calibration. **[Proven]**

A genuine Wick bridge remains conditional on an OS-positive temporal field theory. **[Moderate evidence]**

## 11. Confidence assessment

- Entire common family \(e^{-zA}\): **[Proven]**
- Physical equivalence from analytic continuation alone: **[Refuted]**
- Exponential inverse ill-conditioning: **[Proven]**
- Future OS-positive FIN reconstruction: **[Speculative]**
- Probability of an OS-positive completion: \(0.15\)
- Impact: \(9/10\)
- Difficulty: \(9/10\)
- Time: 6–18 months for a full OS test
- Dependencies: Programs 4 and 8

---

# Program 13 — Short-time laws and calibrated model discrimination

## 1. Motivation

If the operator does not select a temporal world, an observer should discriminate the worlds experimentally. This program derives the earliest measurable difference and tests whether unknown clock calibration can hide it.

**Research hypothesis.** Localized short-time data distinguish coherent and Markov evolution through quadratic versus linear escape, but a single uncalibrated time can be misleading. **[Strong evidence]**

Estimated success \(0.95\), impact \(9/10\), difficulty \(5/10\), time 1–3 months.

## 2. Background

For localized initial vertex \(j\), define

\[
q_i(t)=|U_t(i,j)|^2,\qquad
p_i(t)=P_t(i,j).
\]

The total-variation distance determines the minimum simple-hypothesis error with equal priors:

\[
P_{\rm err}=\frac{1-\operatorname{TV}(p,q)}2.
\]

## 3. Mathematical formulation

The experiment must estimate:

1. escape scaling near zero;
2. the full position distribution at multiple times;
3. a common clock or separate nuisance scales;
4. likelihood under both dynamics.

The key nuisance is that \(p(ct)\) with fitted \(c\) can imitate one unitary snapshot more closely than \(p(t)\).

## 4. Research methodology

Taylor expansions were derived entry by entry. The variance coefficient was computed from the QW2118 weights. Total variation was scanned over time, then reoptimized after allowing a free diffusive time scale. Single- and two-time designs were compared.

## 5. Detailed derivations

For \(i\ne j\),

\[
P_t(i,j)=W_{ij}t+O(t^2),
\]

while

\[
|U_t(i,j)|^2=W_{ij}^2t^2+O(t^3).
\]

Survival obeys

\[
p_j(t)=1-st+O(t^2),
\]

\[
q_j(t)
=1-\sigma_W^2t^2+O(t^4),
\qquad
\sigma_W^2=\sum_{i\ne j}W_{ij}^2=0.537963580.
\]

**[Proven]**

A clock reparametrization \(\tau=ct^2\) could match edge \(ij\) only if \(c=W_{ij}\). Because FIN has six distinct radial weights, no single \(c\) matches all channels. **[Proven]**

## 6. Proof attempts

With full calibrated semigroup data near zero,

\[
A=-\left.\frac{dP_t}{dt}\right|_0,\qquad
A=i\left.\frac{dU_t}{dt}\right|_0.
\]

The generator is unique in either category. **[Proven]**

At one sampled unitary time, Hamiltonian eigenvalues are aliased modulo \(2\pi/t\). Multiple incommensurate times plus a spectral bound can remove this branch ambiguity. **[Proven]**

Population-only unitary data remain invariant under diagonal vertex-phase gauge and global energy shifts, so phase-sensitive interferometry is needed for full identification. **[Proven]**

## 7. Counterexamples

Uniform preparation makes both population distributions stationary and completely indistinguishable. **[Proven]**

With a freely fitted diffusion clock, one short-time snapshot becomes nearly degenerate: at unitary \(t=0.10\), total variation fell from \(0.145276\) under a common clock to \(0.001595\). **[Strong evidence]**

A visually similar spreading profile therefore does not establish dynamical equivalence. **[Refuted]**

## 8. Numerical investigations

For common calibrated time:

| \(t\) | total variation | Bayes error |
|---:|---:|---:|
| 0.10 | 0.145276 | 0.427362 |
| 0.50 | 0.403905 | 0.298048 |
| 0.595 | 0.411241 | 0.294380 |
| 2.00 | 0.262472 | 0.368764 |
| 8.305 | 0.400265 | 0.299868 |

The maximum on \(0.005\le t\le10\) occurred near \(0.595\). With a free diffusive clock, optimized TV values at unitary times \(0.10,0.25,0.50,1.00\) became \(0.00160,0.00875,0.03264,0.10888\). **[Strong evidence]**

Using both \(t\) and \(2t\) breaks the leading degeneracy because coherent escape scales by four and diffusive escape by two. **[Proven]**

## 9. Literature comparison

The quadratic survival law underlies the quantum Zeno effect formalized by [Misra and Sudarshan](https://doi.org/10.1063/1.523304). Classical Markov escape is linear at short time. The FIN-specific contribution is the exact coefficient from QW2118 and the all-edge no-reparametrization test.

## 10. Final conclusions

Calibrated multitime localized measurements can distinguish the two temporal worlds. **[Proven]**

One uncalibrated snapshot can make them nearly observationally equivalent. **[Strong evidence]**

The observer needs a clock, preparation, and outcome model; the operator alone supplies none of them. **[Proven]**

## 11. Confidence assessment

- Linear versus quadratic escape theorem: **[Proven]**
- All-edge clock-reparametrization obstruction: **[Proven]**
- Reported discrimination scan: **[Strong evidence]**
- Single-snapshot identification without calibration: **[Refuted]**
- Probability of a successful simulator discrimination: \(0.90\)
- Impact: \(9/10\)
- Difficulty: \(5/10\)
- Time: 1–3 months
- Dependencies: Program 9's operational map

---

# Program 14 — Multitime records, instruments, and the Chapman–Kolmogorov test

## 1. Motivation

An observer is not only a reader of terminal distributions. Intermediate measurements create records and may alter later evolution. This program asks whether multitime data resolve the ambiguity and whether a history distribution exists independently of the instrument.

**Research hypothesis.** Markov populations obey Chapman–Kolmogorov, whereas Born population matrices from uninterrupted unitary evolution do not; the difference is an operational interference witness. **[Strong evidence]**

Estimated success \(0.98\), impact \(10/10\), difficulty \(5/10\), time 1–3 months.

## 2. Background

Let

\[
M_t(i,j)=|U_t(i,j)|^2.
\]

Although \(U_{t+s}=U_tU_s\), generally \(M_{t+s}\ne M_tM_s\) because amplitudes interfere. The product \(M_t^2\) corresponds to a projective position measurement at the intermediate time with the outcome forgotten.

By contrast \(P_{t+s}=P_tP_s\) exactly.

## 3. Mathematical formulation

Compare:

\[
\mathcal C_U(t)=M_{2t}-M_t^2,\qquad
\mathcal C_D(t)=P_{2t}-P_t^2.
\]

For retained two-step records, compare joint laws

\[
q(x_1,x_2)=M_t(x_1,0)M_t(x_2,x_1),
\]

\[
p(x_1,x_2)=P_t(x_1,0)P_t(x_2,x_1).
\]

## 4. Research methodology

The semigroup identities were proved algebraically. Frobenius and total-variation residuals were computed over several times. Joint record discrimination was compared with terminal-only coarse-graining, and Chernoff exponents were evaluated.

## 5. Detailed derivations

For diffusion,

\[
\mathcal C_D(t)=0
\]

identically. **[Proven]**

For unitary populations,

\[
M_{2t}(i,j)
=\left|\sum_kU_t(i,k)U_t(k,j)\right|^2,
\]

which includes cross terms absent from

\[
(M_t^2)(i,j)
=\sum_k|U_t(i,k)|^2|U_t(k,j)|^2.
\]

Thus \(\mathcal C_U\) measures interference between intermediate alternatives. **[Proven]**

For translation-covariant two-step experiments, the two increments are independent under each measured protocol, so the Chernoff information of the full record is twice the single-increment value. **[Proven]**

## 6. Proof attempts

There is no instrument-independent classical joint distribution for sequential quantum outcomes. A joint probability arises only after specifying the intermediate instrument and its state update. **[Proven]**

Data processing gives

\[
\operatorname{TV}(p_{\rm terminal},q_{\rm terminal})
\le
\operatorname{TV}(p_{\rm record},q_{\rm record}),
\]

because forgetting \(x_1\) is a stochastic channel. **[Proven]**

This theorem explains why an observer who preserves records can discriminate models that appear nearly identical at the endpoint.

## 7. Counterexamples

Interpreting \(M_t\) as a classical transition semigroup is false except in special trivial/permutation cases. **[Refuted]**

The failure of Chapman–Kolmogorov does not mean unitarity lacks composition; amplitudes compose perfectly. It falsifies only an autonomous classical population process without an instrument. **[Proven]**

If the apparatus dephases after every step, the effective law is \(M_t^n\), not uninterrupted \(M_{nt}\). The observer's intervention is part of the dynamics.

## 8. Numerical investigations

Frobenius residuals were

| \(t\) | \(\|\mathcal C_D(t)\|_F\) | \(\|\mathcal C_U(t)\|_F\) |
|---:|---:|---:|
| 0.05 | \(2.8\times10^{-15}\) | 0.01077 |
| 0.10 | \(2.6\times10^{-15}\) | 0.04259 |
| 0.50 | \(1.6\times10^{-15}\) | 0.72319 |
| 1 | \(1.1\times10^{-15}\) | 0.70674 |
| 2 | \(5.8\times10^{-16}\) | 0.64683 |

For two recorded increments, maximum TV on the scanned range was \(0.551749\) near \(t=0.425\). At \(t=5\), record TV was \(0.315212\) while terminal-only TV was \(0.046216\). **[Strong evidence]**

## 9. Literature comparison

Quantum instruments originate in the operational formalism of [Davies and Lewis](https://doi.org/10.1007/BF01647093). Modern process tensors make multitime intervention dependence explicit; see [Pollock et al.](https://doi.org/10.1103/PhysRevLett.120.040405). The noninvasive-measurement issue is central to [Leggett–Garg](https://doi.org/10.1103/PhysRevLett.54.857). The FIN calculation supplies an unusually direct finite witness.

## 10. Final conclusions

The Chapman–Kolmogorov residual is a high-signal observer/intervention test. **[Proven]**

Full records are strictly more informative than terminal marginals in the computed regimes. **[Strong evidence]**

Multitime FIN predictions are undefined until instruments are included in the theory. **[Proven]**

## 11. Confidence assessment

- Markov semigroup identity: **[Proven]**
- Unitary-population interference residual: **[Proven]**
- Reported record advantage: **[Strong evidence]**
- Instrument-independent quantum histories: **[Refuted]**
- Probability of a successful two-time simulator test: \(0.95\)
- Impact: \(10/10\)
- Difficulty: \(5/10\)
- Time: 1–3 months
- Dependencies: Programs 9 and 13

---

# Program 15 — Unitary dilation: diffusion as an observer-restricted quantum world

## 1. Motivation

Could the two temporal worlds be compatible rather than rivals—unitary for a global observer and diffusive for an observer who ignores an environment? This program constructs the compatibility exactly and identifies its cost.

**Research hypothesis.** Every FIN diffusion admits an exact infinite-dimensional unitary dilation, but no finite bath can reproduce irreversible relaxation for all times. **[Strong evidence]**

Estimated success \(0.95\), impact \(9/10\), difficulty \(7/10\), time 2–6 months.

## 2. Background

A contraction semigroup can be represented as a compression of unitary dynamics on a larger Hilbert space. Open quantum dynamics can also reproduce classical Markov populations through Lindblad jump operators.

A dilation explains compatibility, not uniqueness: many environments and couplings can produce the same reduced channel.

## 3. Mathematical formulation

Let

\[
A=\sum_m a_mP_m,\qquad a_m\ge0.
\]

For \(a>0\), use the Cauchy probability measure

\[
d\mu_a(\omega)
=\frac{a}{\pi(a^2+\omega^2)}\,d\omega,
\]

whose characteristic function is \(e^{-a|t|}\). Seek an isometry \(J:H\to\mathcal K\) and self-adjoint bath Hamiltonian \(H_B\) with

\[
J^*e^{-itH_B}J=e^{-tA},\qquad t\ge0.
\]

## 4. Research methodology

The program built the dilation mode by mode, proved the characteristic-function identity, derived a finite-environment no-go theorem from almost periodicity, and independently built a GKSL population embedding from strict edge weights.

## 5. Detailed derivations

Take

\[
\mathcal K=
\bigoplus_m P_mH\otimes L^2(\mathbb R,\mu_{a_m}),
\]

let \(H_B\) multiply by \(\omega\), and embed a system eigenvector as the constant bath function. Then

\[
J^*e^{-itH_B}J
=\sum_me^{-a_m|t|}P_m.
\]

For \(t\ge0\), this equals \(e^{-tA}\). **[Proven]**

A second construction uses jump operators

\[
L_{ij}=\sqrt{W_{ij}}\,|i\rangle\langle j|,\qquad i\ne j.
\]

The Lindblad equation induces on diagonal populations

\[
\dot p_i=\sum_{j\ne i}W_{ij}p_j-sp_i=(Qp)_i.
\]

Thus FIN diffusion is exactly the population sector of a valid Markovian open quantum system. **[Proven]**

## 6. Proof attempts

**Theorem 9 — finite-bath obstruction.** A nonconstant exponentially relaxing semigroup cannot be the exact compression of a finite-dimensional unitary group for all \(t\ge0\). **[Proven]**

**Proof.** Finite-dimensional unitary matrix elements are finite sums \(\sum_kc_ke^{-i\omega_kt}\) and hence almost periodic. An almost-periodic function with a limit at infinity is constant. But \(e^{-at}\to0\) for \(a>0\) and is nonconstant. \(\square\)

An exact irreversible realization therefore needs an infinite bath, continual reset, a singular limit, or time-dependent control.

## 7. Counterexamples

The dilation is not a proof that the universe uses the constructed bath. The same reduced semigroup has nonunique Hamiltonian dilations, jump decompositions, environmental records, and unravelings. **[Proven]**

A finite environment may approximate decay for a finite window but must recur. Treating a short simulation without recurrence as an exact irreversible model is false. **[Refuted]**

Diffusion as a compression is not the same map as \(e^{-itA}\) on the original twelve-dimensional system. **[Proven]**

## 8. Numerical investigations

The QW2118 rates are nonnegative, so all Cauchy measures in the explicit dilation are well defined. Direct matrix exponentiation verified the stochastic row sums, positivity, and the GKSL diagonal population equation to machine precision. **[Strong evidence]**

The largest rate \(2.34218\) sets the shortest decay time, while an exact bath requires a continuous frequency support. Truncating that support produces visible long-time approximation error and recurrence. **[Strong evidence]**

## 9. Literature comparison

The construction is a concrete specialization of unitary dilation theory; see [Schäffer](https://doi.org/10.1090/S0002-9939-1955-0068740-7). Markovian quantum generators were classified independently by [Gorini, Kossakowski and Sudarshan](https://doi.org/10.1063/1.522979) and [Lindblad](https://doi.org/10.1007/BF01608499). FIN contributes a particular shared spectral generator, not a new dilation theorem.

## 10. Final conclusions

Global unitarity and observer-level diffusion are exactly compatible after adding an environment and restriction map. **[Proven]**

Exact irreversible relaxation cannot come from a finite closed bath for all time. **[Proven]**

The operator does not select the environment, monitored algebra, or record; those remain physical axioms. **[Proven]**

## 11. Confidence assessment

- Exact infinite-dimensional dilation: **[Proven]**
- Exact GKSL population embedding: **[Proven]**
- Exact finite-bath realization for all time: **[Refuted]**
- Natural FIN-selected bath: **[Speculative]**
- Probability of useful finite-window laboratory dilation: \(0.80\)
- Impact: \(9/10\)
- Difficulty: \(7/10\)
- Time: 2–6 months
- Dependencies: Programs 6 and 14

---

# Program 16 — Measurement, decoherence, and the quantum-Zeno falsification

## 1. Motivation

A common observer narrative says that measurement or decoherence turns the wave into diffusion. This program tests that statement on the exact FIN operator rather than accepting it qualitatively.

**Research hypothesis.** Projective monitoring generally does not yield the FIN heat semigroup; in the frequent-measurement limit it freezes the walk. **[Strong evidence]**

Estimated success \(0.95\), impact \(8/10\), difficulty \(5/10\), time 1–3 months.

## 2. Background

After a unitary step \(\delta\), unread projective position measurement maps populations by

\[
M_\delta(i,j)=|U_\delta(i,j)|^2.
\]

After \(n\) such steps over total time \(T\), populations are \(M_{T/n}^np_0\). The quantum Zeno effect follows from quadratic short-time survival.

## 3. Mathematical formulation

Compare three protocols:

1. uninterrupted unitary survival \(q_0(T)\);
2. \(n\) repeated projective tests,
   \[
   S_n^Q(T)=q_0(T/n)^n;
   \]
3. a classical no-jump event under rate \(s\),
   \[
   S^D(T)=e^{-sT}.
   \]

Also compare \(M_{T/n}^n\delta_0\) with \(e^{-TA}\delta_0\).

## 4. Research methodology

Short-time expansions were proved, Zeno limits were derived, repeated-measurement matrices were exponentiated, and survival and entropy were scanned for \(T=2,5\) over increasing \(n\).

## 5. Detailed derivations

Because

\[
q_0(\delta)=1-\sigma_W^2\delta^2+O(\delta^4),
\]

\[
S_n^Q(T)
=\exp\!\left[-\frac{\sigma_W^2T^2}{n}+O(n^{-3})\right]
\longrightarrow1.
\]

**[Proven]**

Moreover,

\[
M_\delta=I+O(\delta^2),
\]

so

\[
M_{T/n}^n\longrightarrow I,
\]

not \(e^{-TA}\). **[Proven]**

The passive classical probability of no jump tends to \(e^{-sT}<1\). Measurement creates a different quantum process rather than merely updating knowledge.

## 6. Proof attempts

Could a finite dephasing rate yield diffusion? A GKSL model with specified dephasing and hopping may have a diffusive effective limit, but its rates depend on coupling strengths and scaling. It is not determined by \(W\) alone. **[Moderate evidence]**

No universal map “decohere \(e^{-itA}\)” was found that returns exactly \(e^{-tA}\) with the same time parameter and no environmental constants. **[Strong evidence]**

## 7. Counterexamples

At \(T=5\), increasing measurements from \(n=1\) to \(n=2\) decreases survival before the asymptotic increase. Therefore the claim that every additional finite measurement monotonically enhances survival is false; the theorem concerns the large-\(n\) limit. **[Refuted]**

Unread measurement is not noninvasive. It removes interference and changes subsequent probabilities. **[Proven]**

“Decoherence automatically produces FIN diffusion” is falsified by the projective model. **[Refuted]**

## 8. Numerical investigations

| \(T\) | \(n\) | measured quantum survival | uninterrupted survival | classical no-jump |
|---:|---:|---:|---:|---:|
| 2 | 10 | 0.806328 | 0.199719 | 0.040450 |
| 2 | 200 | 0.989298 | 0.199719 | 0.036326 |
| 2 | 1000 | 0.997850 | 0.199719 | 0.036170 |
| 5 | 10 | 0.259926 | 0.142549 | 0.000526 |
| 5 | 1000 | 0.986641 | 0.142549 | 0.000250 |

The analytic classical limits are \(e^{-2s}=0.0361306\) and \(e^{-5s}=0.0002481\). **[Strong evidence]**

At \(T=2,n=1000\), repeatedly dephased population retained \(0.997851\) at the origin and had TV distance \(0.862727\) from FIN diffusion. **[Strong evidence]**

## 9. Literature comparison

The quantum-Zeno theorem is due to [Misra and Sudarshan](https://doi.org/10.1063/1.523304). Different environmental generators yield different quantum-walk decoherence limits; an example is [Bressanini, Benedetti and Paris](https://doi.org/10.1007/s11128-022-03647-x). The executed FIN counterexample makes the dependence concrete.

## 10. Final conclusions

Frequent projective observation freezes the FIN quantum walk instead of producing its heat flow. **[Proven]**

A diffusive reduced law requires a particular open-system mechanism and scaling, not “observation” in the abstract. **[Proven]**

The observer must be represented by instruments and couplings, never by an undefined interpretive label. **[Proven]**

## 11. Confidence assessment

- Zeno limit for QW2118: **[Proven]**
- Projective-decoherence-to-heat claim: **[Refuted]**
- Alternative finite-rate environmental diffusion: **[Moderate evidence]**
- Probability of laboratory Zeno verification on a simulator: \(0.80\)
- Impact: \(8/10\)
- Difficulty: \(5/10\)
- Time: 1–3 months
- Dependencies: Programs 14 and 15

---

# Program 17 — Entropy, records, and the arrow of time

## 1. Motivation

Diffusion appears irreversible and unitary evolution reversible. This program asks whether entropy or records let the operator itself select a time arrow and what exactly an observer contributes.

**Research hypothesis.** The arrow comes from the contractive temporal branch, a nonequilibrium preparation, and coarse-grained records—not from the real spectral operator alone. **[Strong evidence]**

Estimated success \(0.95\), impact \(9/10\), difficulty \(6/10\), time 2–4 months.

## 2. Background

Three entropies must be separated:

1. von Neumann entropy of a quantum state;
2. Shannon entropy of outcomes in a chosen basis;
3. thermodynamic or stochastic entropy production of a process.

They coincide only under additional hypotheses.

## 3. Mathematical formulation

For unitary \(\rho_t=U_t\rho_0U_t^*\),

\[
S_{\rm vN}(\rho_t)=S_{\rm vN}(\rho_0).
\]

For position outcomes \(q(t)\), study \(H(q(t))\). For reversible Markov \(p(t)=P_tp_0\), study \(H(p(t))\), relative entropy to \(\pi=1/12\), and pathwise entropy production.

## 4. Research methodology

Entropy identities were proved from unitary invariance and doubly stochastic majorization. Time traces were scanned at fine resolution. Forward and reversed path probabilities were compared at equilibrium. Records were classified as system, apparatus, or environment variables.

## 5. Detailed derivations

Unitary conjugation preserves the eigenvalues of \(\rho\), hence its von Neumann entropy. **[Proven]**

For \(t_2>t_1\),

\[
p(t_2)=P_{t_2-t_1}p(t_1).
\]

Because \(P_t\) is doubly stochastic, \(p(t_2)\) is majorized by \(p(t_1)\), so Shannon entropy is nondecreasing. **[Proven]**

At uniform equilibrium, symmetry gives detailed balance

\[
\pi_iP_t(j,i)=\pi_jP_t(i,j).
\]

Every path and its time reverse then have equal probability, so stationary entropy production is zero. **[Proven]**

## 6. Proof attempts

A monotone arrow can be obtained by choosing the \(e^{-tA}\) branch and a nonequilibrium localized initial condition. This is a coherent conditional construction. **[Proven]**

No proof can extract the sign of time from the real self-adjoint \(A\): both \(U_t\) and \(U_{-t}\) are allowed, while the contraction choice itself is extra. **[Proven]**

A low-entropy boundary condition or environmental record can make an arrow physically meaningful, but neither follows from the equilibrium operator. **[Strong evidence]**

## 7. Counterexamples

Quantum position entropy can rise and fall while the global state stays pure. Treating it as irreversible thermodynamic entropy is false. **[Refuted]**

Equilibrium reversible diffusion has stationary path law and zero entropy production despite its contractive generator. Therefore diffusion alone does not imply a nonzero thermodynamic arrow. **[Proven]**

An observer's memory can increase entropy through coupling, but an unspecified observer cannot be used as a mathematical source of irreversibility. **[Refuted]**

## 8. Numerical investigations

On \(0\le t\le20\) with step \(0.005\):

- 1,889 steps decreased unitary position entropy;
- zero steps decreased Markov position entropy beyond numerical tolerance;
- the first unitary decrease occurred near \(t=3.185\);
- \(H_U(3.18)=2.338515\), later falling to \(H_U(16.54)=1.670900\);
- \(H_D(20)=\log12\) within \(8\times10^{-14}\).

**[Strong evidence]**

These results sharply separate measurement entropy from state entropy and semigroup mixing.

## 9. Literature comparison

Network entropy production is formalized by [Schnakenberg](https://doi.org/10.1103/RevModPhys.48.571); quantum dynamical entropy production requires an open-system generator and state. FIN contributes a symmetric finite test case, not a new thermodynamic arrow.

## 10. Final conclusions

The strict spectral generator is time-reversal compatible. **[Proven]**

Monotone mixing arises only after selecting the Markov branch and a nonequilibrium preparation. **[Proven]**

Physical irreversibility additionally needs an operational record, bath, boundary condition, or coarse-graining. **[Strong evidence]**

## 11. Confidence assessment

- Unitary von Neumann entropy constancy: **[Proven]**
- Markov Shannon monotonicity here: **[Proven]**
- Intrinsic arrow from \(A\) alone: **[Refuted]**
- FIN-specific thermodynamic arrow after an open-system completion: **[Moderate evidence]**
- Probability of a rigorous record-level entropy model: \(0.70\)
- Impact: \(9/10\)
- Difficulty: \(6/10\)
- Time: 2–4 months
- Dependencies: Programs 6, 15, and 16

---

# Program 18 — The relational observer and apparatus-calibration theorem

## 1. Motivation

The observer paradox can be confused with subjective consciousness. This program removes that ambiguity by representing the observer as a calibrated apparatus with a frame, instrument, and record.

**Research hypothesis.** Both temporal worlds are covariant under simultaneous relabelling, and only system–apparatus relations are observable; apparatus miscalibration quantitatively destroys directional information. **[Strong evidence]**

Estimated success \(0.95\), impact \(9/10\), difficulty \(5/10\), time 1–3 months.

## 2. Background

Program 2 classified relative frames. Here the theorem is extended to dynamical outcome distributions and a noisy orientation calibration.

An apparatus frame is \(f_a=(r_a,\lambda_a)\); the system frame is \(f_s=(r_s,\lambda_s)\).

## 3. Mathematical formulation

For a translation/reflection-covariant base distribution, define the apparatus outcome

\[
p_{\rm obs}(\delta|f_s,f_a)
=p_{\rm base}(r_a+\lambda_a\delta-r_s).
\]

The simultaneous dihedral action changes both frames but should not change relative outcome probabilities.

For an odd probe, let the apparatus report its orientation correctly with probability \(r\).

## 4. Research methodology

All \(576\) system–apparatus frame pairs and 24 dihedral transformations were enumerated for both \(U_t\) and \(P_t\). A binary noisy calibration channel was then applied to a directional probe \(H(h)=W-hJ\).

## 5. Detailed derivations

The complete invariant remains

\[
(\Delta,c)
=\bigl(\lambda_a(r_s-r_a),\lambda_s\lambda_a\bigr).
\]

Because both \(U_t\) and \(P_t\) are functions of the radial \(A\), they commute with every dihedral relabelling. Therefore every apparatus distribution factors through relative coordinates. **[Proven]**

If the ideal directional accuracy is \(a\), calibration reliability \(r\) gives

\[
a_{\rm eff}=ra+(1-r)(1-a).
\]

At \(r=1/2\), \(a_{\rm eff}=1/2\) for every \(a\): the sign information vanishes exactly. **[Proven]**

## 6. Proof attempts

A relational apparatus can identify whether a displacement is clockwise relative to itself and whether the system and apparatus orientations agree. It cannot define an absolute orientation invariant under simultaneous reflection. **[Proven]**

Likewise, an apparatus clock can compare durations but must be calibrated against a physical standard; it does not make a dimensionful second emerge from \(A\). **[Proven]**

No role for consciousness is needed. All observer effects are represented by preparation maps, instruments, frames, and records. **[Proven]**

## 7. Counterexamples

A randomly oriented apparatus with \(r=1/2\) erases all odd information even if the underlying probe is perfect. **[Proven]**

Fixing array index zero as the apparatus origin creates repeatable data in one implementation but is not invariant under relabelling. It is a laboratory convention, not an internal FIN selector. **[Proven]**

An apparatus-relative result cannot discharge a claim that explicitly demands an apparatus-independent absolute direction. **[Refuted]**

## 8. Numerical investigations

Across \(13\,824\) transformed cases, maximum covariance residuals were

| dynamics | maximum residual |
|---|---:|
| unitary | \(2.08\times10^{-16}\) |
| diffusive | \(1.80\times10^{-16}\) |

For \(t=2,h=0.1\),

\[
P(+1)=0.255947607,\qquad P(-1)=0.159296564.
\]

Conditioned ideal directional accuracy was \(a=0.616379\). Calibration reliabilities \(r=0.50,0.60,0.80,0.95,1.00\) yielded mutual information \(0,0.00156,0.01412,0.03189,0.03944\) bits. **[Strong evidence]**

The \(hJ\) term is an external diagnostic, not a strict internal selector source. **[Proven]**

## 9. Literature comparison

The result belongs to resource theories of reference frames; see [Gour and Spekkens](https://arxiv.org/abs/0711.0043) and relational observables in [Loveridge, Miyadera and Busch](https://arxiv.org/abs/1703.10434). The FIN-specific addition is exact covariance for both temporal branches and the calibration-information curve.

## 10. Final conclusions

The observer is mathematically an operational reference system, not a consciousness postulate. **[Proven]**

Relative observables are complete and covariant for both FIN dynamics. **[Proven]**

Calibration determines accessible directional information but cannot create an absolute selector or physical clock scale. **[Proven]**

## 11. Confidence assessment

- Dynamical relational covariance: **[Proven]**
- Calibration information loss: **[Proven]**
- Absolute selector from apparatus relation: **[Refuted]**
- Physical implementation of the odd probe: **[Speculative]**
- Probability of a complete relational reformulation: \(0.75\)
- Impact: \(9/10\)
- Difficulty: \(5/10\)
- Time: 1–3 months
- Dependencies: Programs 2, 9, and 14

---

# Program 19 — Dynamic exponent, continuum scaling, and causal diagnosis

## 1. Motivation

“Wave dynamics” can mean Schrödinger propagation or a second-order hyperbolic wave equation. They have different continuum scaling and causal behavior even when built from the same spatial operator. This program prevents that ambiguity from being mistaken for Lorentz emergence.

**Research hypothesis.** The operator supports several dynamic exponents; \(z=1\) requires a square root/Dirac or second-order law and is not selected by unitarity. **[Strong evidence]**

Estimated success \(0.90\) for classification, \(0.10\) for a FIN-native causal limit, impact \(10/10\), difficulty \(9/10\), time 6–18 months.

## 2. Background

If a continuum spatial symbol satisfies \(A(k)\sim c|k|^2\), then different temporal equations yield different spectral multipliers and scaling between time and length.

A dynamic exponent \(z\) means characteristic time grows as length to the power \(z\).

## 3. Mathematical formulation

Compare four laws:

| law | multiplier | expected \(z\) |
|---|---|---:|
| heat | \(e^{-tA}\) | 2 |
| Schrödinger | \(e^{-itA}\) | 2 |
| wave | \(\cos(t\sqrt A)\), \(\sin(t\sqrt A)/\sqrt A\) | 1 |
| Poisson | \(e^{-t\sqrt A}\) | 1 |

The first and fourth are dissipative; the second and wave propagator are reversible in their appropriate state spaces.

## 4. Research methodology

Low-mode power laws were regressed, cycle relaxation/wave times were computed for \(N=12,\ldots,192\), QW2118 effective exponents were compared across mode pairs, and propagation support was analyzed from power-series coefficients.

## 5. Detailed derivations

For \(g_N=4\sin^2(\pi/N)\),

\[
\tau_{\rm diff}=g_N^{-1}\sim\frac{N^2}{4\pi^2},
\qquad
\tau_{\rm wave}=g_N^{-1/2}\sim\frac{N}{2\pi}.
\]

**[Proven]**

Doubling spatial resolution multiplies diffusive/Schrödinger time by approximately four and wave time by two. The fitted low-momentum powers \(2.0000\) for \(L\) and \(1.0000\) for \(\sqrt L\) confirm the distinction. **[Strong evidence]**

Unitarity constrains norm preservation, not dispersion. \(e^{-itL}\) is unitary with \(z=2\). **[Proven]**

## 6. Proof attempts

A local \(z=1\) continuum target may be built from a doubled Dirac/incidence operator or a second-order equation

\[
\partial_t^2\phi+A_N\phi=0.
\]

It still requires a refinement, clock/length relation, local limit, orientation, and convergence of commutators or Green functions. **[Moderate evidence]**

QW2118 mode-rate estimates gave \(z_{\rm eff}\approx1.064\) from modes \(1,2\) and \(0.870\) from \(1,3\). Their disagreement prevents a stable exponent claim from twelve nonlocal modes. **[Strong evidence]**

## 7. Counterexamples

For irreducible Markov \(P_t\), every entry is positive for every \(t>0\): diffusion has instantaneous support. **[Proven]**

Finite graph unitary and wave matrices are analytic; a distant entry appears at a finite order in \(t\), yielding suppression rather than an exact cone. **[Proven]**

Therefore neither “unitary” nor “wave-like plot” establishes Lorentz causality. **[Refuted]**

The Poisson semigroup has \(z=1\) scaling yet is dissipative and generally nonlocal, showing that \(z=1\) alone is insufficient. **[Proven]**

## 8. Numerical investigations

| \(N\) | \(\tau_{\rm diff}\) | \(\tau_{\rm wave}\) |
|---:|---:|---:|
| 12 | 3.7321 | 1.9319 |
| 24 | 14.6739 | 3.8306 |
| 48 | 58.4444 | 7.6449 |
| 96 | 233.5274 | 15.2816 |
| 192 | 933.8594 | 30.5591 |

The asymptotic ratios are the clearest computational discriminator of temporal universality class. **[Strong evidence]**

## 9. Literature comparison

Dynamic universality classes are systematically classified in [Hohenberg and Halperin](https://doi.org/10.1103/RevModPhys.49.435). Quantum walks on graphs are not automatically relativistic; see [Farhi and Gutmann](https://arxiv.org/abs/quant-ph/9706062). FIN needs an independently justified local \(z=1\) refinement before Lorentz questions become mathematically well posed.

## 10. Final conclusions

One spatial operator supports \(z=1\) and \(z=2\), reversible and dissipative temporal laws. **[Proven]**

The current finite spectrum does not choose among them or determine a causal cone. **[Proven]**

A causal continuum remains a high-risk, high-impact new construction. **[Speculative]**

## 11. Confidence assessment

- Dynamic-law classification: **[Proven]**
- Low-mode cycle scaling: **[Proven]**
- Unitarity as sufficient for Lorentz symmetry: **[Refuted]**
- FIN-native \(z=1\) causal refinement: **[Speculative]**
- Probability of a useful \(z=1\) construction: \(0.35\)
- Probability of a Lorentzian causal limit: \(0.10\)
- Impact: \(10/10\)
- Difficulty: \(9/10\)
- Time: 6–18 months
- Dependencies: Programs 4, 7, and 12

---

# Program 20 — Minimal operational completion and nonselection theorem

## 1. Motivation

The final observer program asks for the smallest object that turns the common generator into empirical predictions and proves which components cannot be removed.

**Research hypothesis.** The minimal physical object is not a richer spectral operator but an operational process tuple; removing any of its typed components restores concrete nonidentifiability. **[Strong evidence]**

Estimated success \(0.90\), impact \(10/10\), difficulty \(7/10\), time 3–6 months for formalization.

## 2. Background

A physical probability is conditional on preparation, dynamics, instrument, apparatus reference, clock, and record rule. Spectral data constrain one member—the generator—but do not define the conditional experiment.

## 3. Mathematical formulation

Define the operational completion

\[
\mathfrak P=
\bigl(
A,\mathcal S,\rho_0,\gamma,
\{\Phi_\tau\},
\{\mathcal I_{a|x}\},
F_{\rm app},
\mathcal R
\bigr),
\]

where:

- \(A\) is the spectral generator;
- \(\mathcal S\) is the state space and probability rule;
- \(\rho_0\) is a preparation;
- \(\gamma\) calibrates dimensionless to physical time;
- \(\Phi_\tau\) is the dynamical channel;
- \(\mathcal I_{a|x}\) are instruments;
- \(F_{\rm app}\) is the apparatus frame;
- \(\mathcal R\) turns outcomes into persistent records.

## 4. Research methodology

Each component was removed in turn. Programs 11–19 supply explicit witness pairs that become observationally indistinguishable or undefined after removal. The two temporal models were then treated as models of a common axiom set to prove logical nonselection.

## 5. Detailed derivations

**Theorem 10 — dynamical nonselection.** Let \(A\ne0\) be finite, positive, with zero mode and graph-Laplacian signs. The axioms concerning \(A\) alone admit at least

\[
\mathcal M_U:\ \psi(t)=e^{-itA}\psi_0
\]

and

\[
\mathcal M_D:\ p(t)=e^{-tA}p_0.
\]

Both satisfy the same spectral premises, but their nontrivial propagators cannot be simultaneously unitary and stochastic. Therefore neither temporal law is a logical consequence of those premises. **[Proven]**

This is a model-theoretic obstruction: another theorem about the same unaugmented \(A\) cannot choose the physics.

In physical units one must write at least

\[
i\hbar\frac{d\psi}{d\tau}=E_0A\psi
\]

or

\[
\frac{dp}{d\tau}=-\kappa Ap.
\]

Only the products \(E_0\tau/\hbar\) or \(\kappa\tau\) are visible to the dimensionless operator. **[Proven]**

## 6. Proof attempts

Necessity witnesses are:

| removed component | witness |
|---|---|
| state space/dynamics | unitary–stochastic intersection theorem |
| clock calibration | \((\gamma,\tau)\mapsto(\gamma/c,c\tau)\) invariance |
| preparation | uniform state hides all population differences |
| instrument | \(M_{2t}\ne M_t^2\) and Zeno effect |
| apparatus frame | random orientation erases odd information |
| record/environment | distinct decoherence laws share terminal marginals |

Each obstruction is independent of cosmological or philosophical interpretation. **[Proven]**

Additional operational axioms can select a category: transition-probability-preserving reversible pure-state dynamics leads toward unitary/antiunitary structure; positive simplex-preserving continuous semigroups lead to Markov generators. The choice precedes the theorem. **[Strong evidence]**

## 7. Counterexamples

Removing \(\rho_0\) permits the uniform state, which makes both worlds stationary. Removing instruments makes sequential quantum probabilities undefined. Removing the clock allows rate re-fitting. Removing the apparatus frame destroys orientation data. **[Proven]**

Adding only “information” does not repair any component because both models admit information measures. **[Refuted]**

Adding only an observer label without a mathematical instrument changes no prediction. **[Refuted]**

## 8. Numerical investigations

The necessity witnesses were executed on QW2118:

- common-generator residuals at machine precision;
- short-time \(O(t)\) versus \(O(t^2)\);
- Chapman–Kolmogorov unitary residual up to \(0.723\);
- exact Markov residual below \(3\times10^{-15}\);
- relational covariance below \(2.1\times10^{-16}\);
- random calibration yielding zero odd mutual information.

**[Strong evidence]**

Together they form a finite acceptance suite for any claimed operational completion.

## 9. Literature comparison

General probabilistic theories explicitly separate states, transformations, and effects. Quantum reconstruction results, such as [Chiribella, D'Ariano and Perinotti](https://arxiv.org/abs/0908.1583), demonstrate that quantum theory follows only after strong operational axioms. Process tensors add multitime instruments. The FIN result is a concrete independence proof showing why the spectral theorem alone cannot replace those axioms.

## 10. Final conclusions

The minimal bridge to physics is the operational tuple \(\mathfrak P\), not a single missing scalar or another function of \(A\). **[Proven]**

The same operator's wave and diffusion shadows make every omitted component experimentally consequential. **[Proven]**

A future FIN theory must choose or derive a temporal category and then survive preregistered multitime tests. **[Strong evidence]**

## 11. Confidence assessment

- Dynamical nonselection theorem: **[Proven]**
- Necessity of operational typing: **[Proven]**
- Recovery of the full tuple from spectral data alone: **[Refuted]**
- Existence of a FIN-native principle that supplies the full tuple: **[Speculative]**
- Probability of formalizing the necessity theorem: \(0.95\)
- Probability current data already determine the physical tuple: \(0.02\)
- Impact: \(10/10\)
- Difficulty: \(7/10\)
- Time: 3–6 months
- Dependencies: synthesis of Programs 2–19

---

# Part V — Ten recommended next studies (Programs 21–30)

These are deliberately **not** reported as completed research. They are the next highest-value programs after the twenty executions above. Their probabilities are planning estimates.

| Rank | Program | Success probability | Impact | Difficulty | Time |
|---:|---|---:|---:|---:|---|
| 21 | intervention-aware process-tensor tomography | 0.85 | 10/10 | 6/10 | 2–6 months |
| 22 | positivity-preserving \(K_N\) refinement | 0.30 | 10/10 | 9/10 | 6–15 months |
| 23 | \(z=1\) Dirac/hyperbolic refinement | 0.35 | 10/10 | 9/10 | 6–18 months |
| 24 | explicit OS-reflection-positivity program | 0.20 | 10/10 | 10/10 | 9–24 months |
| 25 | SPAM- and clock-robust Bayesian identification | 0.80 | 9/10 | 7/10 | 3–8 months |
| 26 | minimal-environment/dilation identifiability | 0.65 | 8/10 | 8/10 | 4–10 months |
| 27 | full relational-observable rewrite audit | 0.75 | 9/10 | 6/10 | 2–5 months |
| 28 | adaptive odd-sector and integrability audit | 0.35 | 9/10 | 9/10 | 6–15 months |
| 29 | machine-checked operational no-go library | 0.90 | 7/10 | 6/10 | 2–6 months |
| 30 | blinded twelve-mode held-out challenge | 0.70 | 10/10 | 8/10 | 6–18 months |

## Program 21 — Intervention-aware process-tensor tomography

**Why now.** Program 14 produced the sharpest finite discriminator, \(M_{2t}\ne M_t^2\). The next step is to reconstruct all two- and three-time intervention probabilities rather than only position projections. **[Strong evidence]**

**Hypothesis.** A modest informationally complete instrument set will separate coherent, Markov, dephasing, and non-Markovian bath models even with realistic noise. **[Moderate evidence]**

- Mathematics: process tensors, quantum combs, identifiability, Chernoff design.
- Computation: optimize times and instruments under finite shots; bootstrap confidence regions.
- Physical experiment: twelve coupled optical, microwave, trapped-ion, or programmable quantum-walk modes; compare runs with and without intermediate intervention.
- Falsification: all candidate models remain likelihood-equivalent after informationally complete data and calibrated nuisance bounds.
- Dependencies: Programs 9, 13–16.
- Estimated existence/success probability: \(0.85\).

## Program 22 — Positivity-preserving FIN refinement

**Why now.** The exact \(C_{12}\) Markov generator loses off-diagonal positivity in the naive \(N=16\) extension. A stochastic continuum is impossible without repairing that failure. **[Proven]**

**Hypothesis.** There may exist a running parameter family \(W_N\ge0\) that matches \(N=12\), retains a nontrivial strict profile, and converges after normalization to a positivity-preserving operator. **[Speculative]**

- Mathematics: completely monotone graph filters, Bernstein functions, Dirichlet forms, Mosco convergence.
- Computation: constrained continuation search with held-out sizes and interval-certified signs.
- Physical experiment: compare relaxation modes on several emulator sizes.
- Falsification: no feasible family satisfies positivity, \(N=12\) fit, parameter simplicity, and stable low-mode limit.
- Dependencies: Programs 1, 4, 11, 19.
- Estimated probability: \(0.30\).

## Program 23 — A \(z=1\) doubled Dirac or hyperbolic refinement

**Why now.** Neither heat nor Schrödinger evolution of a Laplacian establishes Lorentz behavior. A first-order or second-order hyperbolic structure is the shortest credible route. **[Strong evidence]**

**Hypothesis.** A vertex-edge doubled operator can square to a graph Laplacian while admitting a controlled \(z=1\) limit. **[Moderate evidence]**

- Mathematics: discrete exterior calculus, graph Dirac operators, strong-resolvent convergence, Lieb–Robinson bounds.
- Computation: packet velocity, dispersion linearity, commutator tails, robustness under refinement.
- Physical experiment: measure arrival-time scaling \(N\) versus \(N^2\).
- Falsification: persistent nonlinear dispersion, nonlocal tails, or absence of a consistent time orientation.
- Dependencies: Programs 4, 7, 19, 22.
- Estimated probability: \(0.35\).

## Program 24 — Osterwalder–Schrader positivity for explicit FIN correlators

**Why now.** Analytic continuation is physically empty without reflection positivity. This is the hardest rigorous quantum-field gate. **[Proven]**

**Hypothesis.** A temporally extended free theory based on \(A_{\rm M}+m^2\) may be reflection positive; interacting or FIN-selected uniqueness is much less likely. **[Moderate evidence]**

- Mathematics: transfer matrices, reflection-positive kernels, Schwinger functions, constructive field theory.
- Computation: exact finite-lattice reflection matrices and smallest-eigenvalue certificates.
- Physical experiment: none until a reconstructed Hamiltonian yields held-out correlators; later compare Euclidean and real-time simulator data.
- Falsification: a negative reflection form for any declared admissible test function.
- Dependencies: Programs 8, 12, 19, 23.
- Estimated probability: \(0.20\).

## Program 25 — SPAM- and clock-robust Bayesian identification

**Why now.** One uncalibrated snapshot nearly erased wave/diffusion separation. Preparation-and-measurement errors and clock drift must be treated as parameters, not ignored. **[Proven]**

**Hypothesis.** Three or more optimized times plus calibration controls identify the temporal category despite bounded SPAM error. **[Moderate evidence]**

- Mathematics: Bayesian experimental design, structural identifiability, nuisance marginalization.
- Computation: posterior predictive checks, model misspecification, sequential design.
- Physical experiment: randomized order, blind labels, calibration pulses, common-clock and independent-clock controls.
- Falsification: broad nuisance priors absorb every proposed held-out distinction.
- Dependencies: Programs 9, 13, 14, 18.
- Estimated probability: \(0.80\).

## Program 26 — Minimal environment and dilation identifiability

**Why now.** Program 15 proves existence but extreme nonuniqueness of baths. Physical content requires determining which environmental features are inferable from system records. **[Strong evidence]**

**Hypothesis.** Finite-window multitime data can bound bath dimension, memory time, or spectral density, though never identify a unique infinite dilation. **[Moderate evidence]**

- Mathematics: minimal Stinespring dimension, realization theory, pseudomodes, non-Markovianity.
- Computation: fit Cauchy, finite-mode, and reset-bath models under recurrence tests.
- Physical experiment: monitor selected environmental outputs and look for information backflow.
- Falsification: observational equivalence of all bath classes across accessible controls.
- Dependencies: Programs 14–17, 21.
- Estimated probability: \(0.65\).

## Program 27 — Full relational-observable rewrite audit

**Why now.** The relative-frame theorem may remove the absolute selector from physics, but only if every claimed observable and coupling can be rewritten relationally. **[Strong evidence]**

**Hypothesis.** Most operational predictions factor through system–apparatus invariants; any exception will precisely identify the needed odd physical resource. **[Moderate evidence]**

- Mathematics: groupoids, invariant theory, resource theories of asymmetry.
- Computation: type every origin/sign-sensitive artifact and mechanically test covariance.
- Physical experiment: apparatus reflection and random-orientation controls.
- Falsification: a necessary observable remains absolute and no measurable odd resource supplies it.
- Dependencies: Programs 2, 18, 20.
- Estimated probability: \(0.75\).

## Program 28 — Adaptive integrability and odd-sector bifurcation

**Why now.** Feedback gradient claims and spontaneous selector claims remain unproved. The actual adaptive Jacobian must replace generic pitchfork prose. **[Strong evidence]**

**Hypothesis.** Either the covariance one-form fails integrability or an explicit odd eigenvalue crossing can be isolated in a justified \(K_N\) family. **[Speculative]**

- Mathematics: Fréchet derivatives, Helmholtz integrability conditions, equivariant center manifolds, stochastic limits.
- Computation: automatic differentiation, Jacobian representation blocks, interval-certified bifurcation coefficients.
- Physical experiment: perturbation/reversal hysteresis if a real adaptive platform is specified.
- Falsification: no odd tangent source or no transverse eigenvalue crossing.
- Dependencies: Programs 2, 5, 22.
- Estimated probability of a positive FIN-specific bifurcation: \(0.35\).

## Program 29 — Machine-checked operational no-go library

**Why now.** The strongest results are short necessity theorems that should not depend on prose interpretation. **[Strong evidence]**

**Hypothesis.** The selector, scale, finite-bath, unitary–stochastic intersection, and dynamical-nonselection proofs can all be formally certified. **[Strong evidence]**

- Mathematics: finite group actions, linear algebra, representation weights, semigroup axioms.
- Computation: Lean/Isabelle proofs linked to rational or interval-certified QW2118 witnesses.
- Physical experiment: none; the deliverable is auditability.
- Falsification: a hidden premise is found or a checker rejects the proof.
- Dependencies: Programs 1–20.
- Estimated probability: \(0.90\).

## Program 30 — Blinded twelve-mode held-out challenge

**Why now.** FIN needs a prediction test, not another reconstruction from known target values. **[Proven]**

**Hypothesis.** A preregistered subset of strict spectral ratios and multitime responses can predict blinded emulator observations beyond affine calibration and nuisance fitting. **[Moderate evidence]**

- Mathematics: scale-free invariants, likelihood ratios, multiple-testing control.
- Computation: generate blinded synthetic and hardware datasets; lock fitting modes before revealing tests.
- Physical experiment: independently fabricated twelve-mode weighted ring, with both coherent and controlled dissipative regimes.
- Falsification: held-out errors match or exceed generic circulant alternatives.
- Dependencies: Programs 9, 13, 14, 21, 25.
- Estimated probability of a technically decisive emulator result: \(0.70\).
- Estimated probability that the result alone establishes fundamental physics: \(0.02\).

---

# Part VI — Independent comparison, dependency graph, and roadmap

## 1. Objective comparison of the twenty executed programs

Scores are \(0\)–\(10\) judgments. They rank research value, not truth status.

| Overall rank | Program | Depth | Originality | Robustness | Physical relevance | Testability | Transformative potential |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | P20 operational completion/nonselection | 9 | 8 | 10 | 10 | 8 | 10 |
| 2 | P14 multitime instrument test | 8 | 8 | 10 | 10 | 10 | 9 |
| 3 | P9 calibrated predictive challenge | 7 | 7 | 8 | 10 | 10 | 10 |
| 4 | P4 refinement and continuum | 10 | 6 | 9 | 10 | 5 | 10 |
| 5 | P11 exact dual dynamics | 8 | 7 | 10 | 9 | 9 | 9 |
| 6 | P2 relational selector | 8 | 8 | 10 | 9 | 8 | 9 |
| 7 | P13 short-time discrimination | 7 | 7 | 10 | 9 | 10 | 8 |
| 8 | P19 dynamic exponent/causality | 9 | 6 | 9 | 10 | 7 | 9 |
| 9 | P15 unitary dilation | 8 | 6 | 10 | 9 | 8 | 8 |
| 10 | P18 relational apparatus | 7 | 7 | 10 | 9 | 8 | 8 |
| 11 | P3 units/RG/modular time | 9 | 6 | 10 | 10 | 4 | 9 |
| 12 | P17 entropy and time arrow | 7 | 6 | 10 | 9 | 8 | 8 |
| 13 | P8 action and OS gate | 9 | 6 | 9 | 10 | 4 | 9 |
| 14 | P16 decoherence/Zeno | 7 | 6 | 10 | 8 | 9 | 7 |
| 15 | P12 holomorphic time/Wick | 8 | 6 | 10 | 9 | 5 | 8 |
| 16 | P6 information thermodynamics | 8 | 5 | 10 | 9 | 5 | 7 |
| 17 | P7 repaired NCG | 9 | 6 | 9 | 8 | 3 | 7 |
| 18 | P5 adaptive bootstrap | 8 | 7 | 8 | 7 | 4 | 7 |
| 19 | P1 spectral bridge | 7 | 8 | 10 | 5 | 3 | 6 |
| 20 | P10 machine certificates | 6 | 5 | 9 | 5 | 2 | 6 |

**[Strong evidence]** Programs 20, 14, 9, and 4 form the shortest path from mathematical underdetermination to an empirically risky theory. Program 11 supplies the new organizing theorem; Program 2 supplies the best selector reformulation.

## 2. Dependency graph

    Finite strict spectral core
       ├── P1 finite bridge
       ├── P2 relational selector ── P18 relational apparatus ──┐
       ├── P3 units and clock ────────────────┐                 │
       ├── P4 refinement ── P7 geometry ── P8 OS/action        │
       │          └────────────── P19 dynamic exponent          │
       └── P11 dual generator                                  │
              ├── P12 holomorphic time                          │
              ├── P13 short-time test ──────────────────────────┤
              ├── P14 multitime instruments ────────────────────┤
              ├── P15 dilation/bath ── P17 records              │
              └── P16 decoherence/Zeno                          │
                                                               ▼
                                      P20 operational completion
                                                │
                      ┌─────────────────────────┼────────────────────┐
                      ▼                         ▼                    ▼
             P21/P25 tomography       P22/P23/P24 refinement       P30 blind test

The graph has three bottleneck cuts:

1. **operational cut:** no preparation–instrument–record map, no empirical probability;
2. **refinement cut:** no \(K_N\), no continuum, universality, causal limit, or asymptotic dimension;
3. **calibration cut:** no dimensionful clock/length/action source, no numerical physical prediction.

All three must be crossed for fundamental physics. **[Proven]**

## 3. Critical bottlenecks

### B1. Independently motivated empirical domain

A programmable emulator is easy to fit; a natural system selected independently of the kernel is not. Without it, FIN remains mathematics. **[Proven]**

### B2. Refinement and locality

The \(N=12\) Markov property already fails under the naive \(N=16\) profile extension. A valid continuum needs new running laws and positivity/locality control. **[Proven]**

### B3. Temporal category and clock

The common generator supports unitary, heat, wave, and Poisson laws. Choosing one and calibrating its time is separate structure. **[Proven]**

### B4. State, instrument, and record

Uniform preparation hides model differences; intermediate observation changes coherent evolution; records control inferential power. These are not optional semantics. **[Proven]**

### B5. Reflection positivity or an alternative quantum reconstruction

Complex-time analyticity alone does not produce a physical Hilbert space from Euclidean data. **[Proven]**

## 4. High-risk/high-reward directions

- A positivity-preserving, local, \(z=1\), OS-positive refinement would transform the framework, but its joint probability is at most about \(0.10\). **[Speculative]**
- A successful blinded physical-domain test has the largest immediate empirical value and about \(0.20\) probability once a nontrivial domain is identified. **[Moderate evidence]**
- A FIN-native principle deriving the operational tuple has enormous value but no current positive mechanism; probability about \(0.15\). **[Speculative]**
- A repaired graph spectral triple is mathematically feasible but unlikely to force Standard-Model gauge content. **[Moderate evidence]**

## 5. Dead ends and bounded no-go lanes

- extracting an absolute selector from reflection-even scalars: **[Refuted]**
- extracting a nonzero physical unit from dimensionless invariants: **[Refuted]**
- identifying matrix inertia with Lorentz signature: **[Refuted]**
- treating the fixed-size polynomial bridge as a continuum/role bridge: **[Refuted]**
- treating Shannon entropy as heat without Hamiltonian and bath: **[Refuted]**
- using the naive \(D=K_s,J=\) conjugation tuple as a real/even triple: **[Refuted]**
- using \(K_s\) unshifted as a Gaussian Euclidean action: **[Refuted]**
- using free probability for the commuting strict/legacy pair: **[Refuted]**
- assuming projective decoherence produces \(e^{-tA}\): **[Refuted]**
- claiming one terminal spreading profile identifies the dynamics: **[Refuted]**
- continuing scalar or selector scans without a new odd/charged source: **[Refuted]**

## 6. Unexpected discoveries

1. The strict and legacy matrices are exact two-way polynomial functions on \(C_{12}\), invalidating an inertia-only finite obstruction. **[Proven]**
2. The strict QW2118 profile itself yields a positive reversible Markov generator after its Perron shift, while the same shift generates an equivalent unitary walk up to phase. **[Proven]**
3. Both dynamics are restrictions of \(e^{-zA}\), but their physical equivalence is false and inverse continuation is exponentially unstable. **[Proven]**
4. A complete relative selector exists for system–apparatus pairs even though an absolute selector does not. **[Proven]**
5. FIN diffusion has an exact infinite unitary dilation but no exact finite-bath dilation for all time. **[Proven]**
6. Frequent unread position measurement freezes the unitary walk instead of producing FIN diffusion. **[Proven]**
7. Multitime records can retain large discriminating information when the terminal marginal is almost useless. **[Strong evidence]**

## 7. New theorem inventory

The monograph establishes or specializes the following:

- finite two-way kernel interpolation theorem;
- finite bridge nonuniqueness and generalization counterexample;
- natural absolute-selector obstruction;
- complete relative-frame invariant theorem;
- dimensionful-source equivariance theorem;
- single-snapshot continuum obstruction;
- feedback-gradient correction theorem;
- Gibbs energy-scale nonidentifiability theorem;
- naive spectral-triple first-order and grading obstructions;
- strict action stability threshold;
- invariant-to-odd and weight-zero source exclusions;
- exact QW2118 dual unitary/Markov theorem;
- stochastic–unitary trivial-intersection theorem;
- holomorphic temporal-family and inverse-conditioning result;
- all-edge short-time clock-reparametrization obstruction;
- multitime Chapman–Kolmogorov interference witness;
- exact infinite dilation and finite-bath no-go;
- projective-decoherence/Zeno limit theorem;
- entropy/arrow separation;
- dynamical nonselection theorem and operational necessity witnesses.

Each theorem is **[Proven]** in its stated finite or conditional scope. None proves that FIN is a physical theory.

## 8. New conjectures

### Conjecture C1 — relational completeness

Every empirically admissible FIN observable can be written as a diagonal system–apparatus invariant; any nonrelational residue corresponds to a measurable odd resource. **[Speculative]** Estimated probability \(0.75\).

### Conjecture C2 — positive refinement

A nontrivial parameter-running \(W_N\) exists that matches QW2118, preserves nonnegative Dirichlet weights, and has a stable continuum low-mode limit. **[Speculative]** Estimated probability \(0.30\).

### Conjecture C3 — no natural temporal ray

No automorphism-natural rule on an unadorned positive self-adjoint operator selects the imaginary rather than positive real time ray while remaining compatible with both unitary and Markov functional calculi. **[Strong evidence]** A categorical formulation is still required.

### Conjecture C4 — hyperbolic reconstruction

A doubled graph-Dirac refinement can produce \(z=1\) dispersion, but satisfying locality, orientation, OS positivity, and FIN-specific parameter selection simultaneously is rare. **[Speculative]** Estimated probability \(0.15\).

### Conjecture C5 — adaptive odd bifurcation obstruction

The current feedback law lacks either an exact covariance one-form or a strict odd eigenvalue crossing, so it cannot internally select a branch without new data. **[Moderate evidence]**

## 9. Phased roadmap

### Phase 0 — formal freeze, 0–3 months

- reproduce every finite number with two independent implementations;
- formalize the selector, unit, dual-dynamics, and nonselection theorems;
- preregister admissible inputs and nuisance parameters.

Success criterion: machine-checkable statements and no provenance leakage.

### Phase 1 — operational observer tests, 2–6 months

- execute two- and three-time intervention tomography;
- fit coherent, Markov, and open-system alternatives;
- include SPAM and clock uncertainty;
- retain complete records.

Success criterion: Bayes factors and predictive residuals on held-out runs.

### Phase 2 — relational and dimensional completion, 3–9 months

- rewrite every observable relative to an apparatus;
- declare the rank-three calibration package;
- determine whether any absolute odd claim remains.

Success criterion: every predicted number has units, preparation, measurement, and symmetry behavior.

### Phase 3 — refinement competition, 6–18 months

- construct several \(K_N\) families, including positivity-preserving and \(z=1\) candidates;
- compare them at unseen \(N\);
- reject any without stable low-mode and locality behavior.

Success criterion: a convergent family selected by a non-post-hoc criterion.

### Phase 4 — field reconstruction, 9–24 months

- test reflection positivity and transfer matrices;
- audit repaired spectral triples and gauge groups;
- derive correlators and effective actions.

Success criterion: a rigorous reconstructed dynamical theory with held-out observables.

### Phase 5 — blinded physical challenge, 12–36 months

- select an independent physical platform;
- lock calibration modes and predictions;
- reveal held-out spectral and multitime data;
- compare against generic circulant, Markov, and quantum-walk alternatives.

Success criterion: repeatable out-of-sample superiority, not mere implementability.

---

# Part VII — Final scientific judgment

## 1. Deepest interpretation that survives falsification

The smallest robust object is the positive Dirichlet generator

\[
A=sI-K_s.
\]

Its spectral measure generates a holomorphic family \(e^{-zA}\). The unitary, heat, wave, Poisson, Green, and spectral-filter constructions are different functional or temporal images of this generator. **[Proven]**

This does **not** mean they are the same physics. Their state spaces, composition laws, probability rules, interventions, entropy behavior, and continuum scaling differ. **[Proven]**

The deepest surviving interpretation is therefore:

> FIN is presently a finite spectral kinematics with multiple admissible temporal semantics. Its missing physical principle is an operational-relational completion that chooses a state cone and temporal category, calibrates time and scale, specifies preparation/instruments/records, and embeds the finite generator into a convergent family.

**[Strong evidence]**

The observer paradox disappears in this formulation. A phase-sensitive observer with coherent control sees unitary interference; a population observer coupled to or ignoring an environment may see diffusion. Neither observer creates the operator. Their accessible algebra and intervention channel determine which predictions are operationally available. **[Proven]**

## 2. Current verdict

- Mathematical coherence of the finite spectral core: **[Strong evidence]**
- Exact strict finite unitary/Markov duality after Perron shift: **[Proven]**
- Automatic emergence of physics from the spectral theorem: **[Refuted]**
- Current status as a predictive fundamental physical theory: **[Refuted]**
- Possibility of a calibrated, testable extension: **[Moderate evidence]**
- Claim that consciousness is needed to resolve the two dynamics: **[Refuted]**
- Claim that preparation, clock, instrument, environment, and records are needed for empirical probabilities: **[Proven]**

The appropriate scientific response is neither confirmation nor dismissal. It is to stop extending receiver formulas and attack the operational, refinement, and calibration cuts with blinded tests.

---

# What are the three discoveries that would most dramatically increase the probability that FIN becomes a genuine theory of physics?

## Discovery 1 — A FIN-native causal and reflection-positive refinement theorem

**Statement sought.** Construct a non-post-hoc family \(A_N\) matching QW2118 at \(N=12\), preserving positive Dirichlet structure, possessing a local \(z=1\) or otherwise physically justified limit, and yielding reflection-positive Euclidean correlators or a direct unitary continuum reconstruction.

**Why it is fundamental.** It would replace the isolated matrix with a continuum-capable dynamical theory and simultaneously address universality, locality, causality, and the wave/diffusion relation. Without it, twelve-point behavior cannot define spacetime or field physics. **[Proven as a necessity assessment]**

**Does it appear achievable?** Individual components are standard; their simultaneous FIN-specific derivation is difficult and currently unsupported. **[Speculative]**

**Suggested mathematical techniques.**

- Dirichlet-form and Mosco convergence;
- discrete exterior calculus and graph Dirac operators;
- strong-resolvent convergence;
- Lieb–Robinson estimates and hyperbolic principal symbols;
- transfer matrices and Osterwalder–Schrader positivity;
- equivariant orientation and refinement functors.

**Suggested computational experiments.**

- constrained \(K_N\) parameter-running search with interval-certified positivity;
- multi-size dispersion and dynamic-exponent fits;
- commutator-tail and arrival-time scaling;
- reflection-form eigenvalue tests;
- out-of-size bridge validation.

**Suggested physical experiments.**

- implement several weighted-ring sizes;
- compare \(N\) versus \(N^2\) propagation times;
- test real-time and Euclidean correlators under the same preregistered parameters;
- search for recurrence, locality violations, and bath dependence.

**Estimated probability that the full discovery exists:** \(0.10\).

## Discovery 2 — A blinded, independently motivated physical-domain match

**Statement sought.** Identify a natural physical system for reasons independent of the FIN fit, calibrate only a preregistered subset of parameters, and correctly predict at least two held-out spectral ratios plus multitime intervention data better than generic circulant, Markov, quantum-walk, and flexible nuisance alternatives.

**Why it is fundamental.** A theory becomes physics by surviving risky observation. This discovery would supply the absent referent of vertices, clock, preparations, and outcomes and would sharply increase empirical credibility even before a full continuum theory. **[Proven as a methodological necessity]**

**Does it appear achievable?** A simulator challenge is feasible; finding a non-engineered domain and beating alternatives is substantially harder. **[Moderate evidence]**

**Suggested mathematical techniques.**

- operational model mapping;
- structural identifiability;
- affine-invariant spectral statistics;
- process-tensor likelihoods;
- Bayesian model comparison and posterior predictive checking;
- multiple-testing and preregistration discipline.

**Suggested computational experiments.**

- blind synthetic-data challenges with clock and SPAM uncertainty;
- sequential Chernoff-optimal design;
- robustness to disorder, mode relabelling, dephasing, and missing data;
- strict held-out comparison of bridge and refinement families.

**Suggested physical experiments.**

- twelve-mode coherent/dissipative platform with switchable intermediate measurements;
- short-time escape-law test;
- Chapman–Kolmogorov versus interference test;
- apparatus-reflection controls;
- blinded held-out mode spectroscopy.

**Estimated probability that a technically decisive simulator discovery exists:** \(0.70\).

**Estimated probability that an independently motivated fundamental-domain match exists:** \(0.20\).

## Discovery 3 — A source-bearing operational selection principle

**Statement sought.** Find a new axiom or theorem—independent of the current spectral data—that uniquely or nearly uniquely supplies the admissible state space, temporal ray, apparatus relation, environment/record structure, and dimensional calibration, while yielding testable restrictions not fitted afterward.

**Why it is fundamental.** The dynamical nonselection theorem proves that \(A\) alone admits inequivalent unitary and diffusive worlds. A source-bearing operational principle is the smallest kind of result capable of selecting physics without pretending that another scalar function of \(A\) can do so. **[Proven as a logical necessity]**

**Does it appear achievable?** Relational reference frames and open-system dilations provide partial ingredients, but no known principle fixes the whole tuple or its scales. **[Speculative]**

**Suggested mathematical techniques.**

- convex operational reconstructions and purification axioms;
- categorical process theories and groupoids of reference frames;
- \(C^*\)- and von Neumann algebraic state selection;
- modular/local covariance with explicit clock calibration;
- resource theories of asymmetry and thermodynamics;
- independence proofs and categorical naturality.

**Suggested computational experiments.**

- automated model search over minimal axiom subsets;
- countermodel generation for every proposed selection principle;
- formal verification of necessity and independence;
- prediction comparison on unseen graphs and instrument sets.

**Suggested physical experiments.**

- compare coherent, diffusive, and open-system predictions under identical spatial couplings;
- reverse the apparatus frame and environmental monitoring;
- measure full multitime records and entropy flow;
- test whether the proposed principle forbids any actually realizable channel.

**Estimated probability that such a noncircular FIN-specific principle exists:** \(0.15\).

---

## Final answer in one sentence

The same FIN operator genuinely generates wave/unitary and diffusive mathematical shadows, but the object that could turn either shadow into physics is not another spectral formula: it is a calibrated operational-relational continuum theory whose temporal law and held-out observables survive experiment. **[Strong evidence]**

