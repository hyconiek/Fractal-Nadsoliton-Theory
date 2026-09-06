# ST8547 — Entropy rigidity under binary refinement and hidden correlation

Krzysztof Żuchowski — Independent Researcher — Fractal Information Theory Project  
ORCID 0009-0002-0909-3613  
5 September 2026

## Scope and evidence

This study follows ST8397–ST8546. It supplies explicit proofs, an adversarial
counterfamily and independently recomputed numerical tests. It does not count
JSON status records as separate studies or as proof certificates.

The authoritative local inputs consulted were the project rules in AGENTS.md,
the current-state compendium, the kernel-transformation diagram, the
ST8397–ST8546 source and report, and the earlier product-refinement results
indexed in AGENTS.md (ST30, ST245, ST1137). The broader project history enters
through those source summaries. This is not a fresh verification of every
historical report or all numerical certificates in the repository.

The previous 150 programme records include 144 files with an empty JSON payload.
Their hashes certify file integrity. They do not independently certify the
associated mathematical claims. The logarithmic-mean identity itself is valid
and is established below; several interpretations need narrower hypotheses.

Every theorem below is proved in its stated finite-graph or scalar class.
Numerical values are floating-point checks, not interval certificates.
Literature priority is discussed at the end. No claim of globally new
mathematics or experimentally established FIN physics is made.

## 1. Inputs and conventions

Let W be a real symmetric conductance matrix with nonnegative off-diagonal
entries, zero diagonal, and connected support. Write

\[
A=\operatorname{diag}(W\mathbf1)-W.
\]

For the frozen strict FIN example there are twelve vertices and

\[
W_{ij}=\frac{\cos(0.18575d+0.16250)}{1+d^{9/5}},
\quad d=\min(|i-j|,12-|i-j|),\quad i\ne j.
\]

All these weights are positive since the cosine argument at d=1,...,6 lies
between zero and pi/2. The strict parameters are given inputs; this study
does not derive them. The historical raw legacy kernel has signed weights
and cannot be substituted into the positive-conductance theorems below.

Use a positive mass vector p, and a separable convex functional

\[
F_\phi(p)=\sum_i\phi(p_i),\qquad
\phi\in C^2(0,\infty),\quad\phi''>0.
\]

Universal statements use the positive mass cone. For probability states we
also impose sum p_i=1. We use the same scalar phi at each refinement level,
in mass coordinates, without changing a reference measure or clock. These
conventions are part of the theorem, not consequences of the spectral theorem.

For a continuous positive symmetric scalar mobility M define

\[
(K^W_M(p)\xi)_i=\sum_jW_{ij}M(p_i,p_j)(\xi_i-\xi_j).
\]

We require heat matching edge by edge, for all positive a,b:

\[
M(a,b)[\phi'(a)-\phi'(b)]=a-b. \tag{1}
\]

This is stronger than matching one vector field at one state. It follows, for
example, if the same rule must work for every single-edge graph. It should not
be inferred solely by inspecting a fixed complete graph at a single p.

## 2. Established gradient identity and its scope

Equation (1) uniquely determines

\[
M_\phi(a,b)=\frac{a-b}{\phi'(a)-\phi'(b)},\qquad
M_\phi(a,a)=\frac1{\phi''(a)}. \tag{2}
\]

The diagonal is the continuous limit. Summing edge fluxes proves

\[
K^W_{M_\phi}(p)\nabla F_\phi(p)=Ap. \tag{3}
\]

For phi(x)=x log x the mobility is the logarithmic mean

\[
\Lambda(a,b)=\int_0^1a^t b^{1-t}\,dt
=\frac{a-b}{\log a-\log b}. \tag{4}
\]

Consequently negative Shannon entropy F=sum p log p decreases along
p_dot=-Ap; Shannon entropy S=-F increases. Its dissipation is

\[
\mathcal I_W(p)=\sum_{i<j}W_{ij}(p_i-p_j)\log(p_i/p_j)\ge0. \tag{5}
\]

This gradient representation is established mathematics for reversible Markov
chains, not a discovery originating in FIN. Maas 2011 treats the logarithmic
mean and the general divided-difference entropy mobility; see (2.7) and the
heat-gradient theorem in [Maas's original paper](https://arxiv.org/abs/1102.5238).

The construction recovers an Onsager representation of an already supplied
generator A. It does not infer unknown microscopic jump rates or a physical
time unit. Nor does (3) determine an unrestricted matrix K: at a given p,
adding a positive operator annihilating both 1 and log p preserves (3).
Such modifications generally violate the universal edge-local ansatz.

Boundary states need a limiting interpretation: Lambda(a,0)=0, while log(0)
is not finite. The formal product K_p log p at such a state is undefined
term by term. The limiting edge flux a-b and the linear heat equation remain
well defined.

## 3. What refinement compatibility must mean

On an n-by-m product use horizontal weights W and vertical weights B:

\[
\widetilde A=A\otimes I_m+I_n\otimes L_B,\qquad
C=I_n\otimes\mathbf1_m^\top.
\]

C sums masses. It is distinct from the Hilbert isometry
R=I_n tensor 1_m/sqrt(m). A product probability state is
P_{ib}=p_i r_b, where r_b>0 and sum r_b=1.

Linear heat satisfies C Atilde=A C for every fine state and any B. A stronger,
separate requirement is preservation of the full cotangent quadratic form:

\[
C K^{\widetilde W}_M(P)C^\top=K^W_M(p). \tag{6}
\]

Lifted potentials C^T xi are constant within each fiber, so vertical edges
make zero contribution to (6). For a horizontal edge it becomes

\[
\sum_b M(r_b a,r_b b)=M(a,b). \tag{7}
\]

In particular, equal binary splitting requires

\[
2M(a/2,b/2)=M(a,b). \tag{8}
\]

The heat intertwiner alone does not imply (6). Formula (3) holds for every
admissible phi, including examples that fail (6). This is the exact additional
compatibility premise needed for entropy selection.

## 4. Binary refinement alone does not select Shannon

**Theorem 1 — complete dyadic classification in the stated class.**
Under (1), equation (8) holds for all positive masses if and only if

\[
\phi''(x)=\frac{g(\log_2 x)}{x}, \tag{9}
\]

where g is continuous, strictly positive and one-periodic.

**Proof.** Taking a=b=x in (8) and using (2) gives
phi''(x/2)=2 phi''(x). Thus h(x)=x phi''(x) is invariant under doubling,
which is precisely (9). Conversely integrate (9) from b/2 to a/2.
Substitution x=u/2 and periodicity give
phi'(a/2)-phi'(b/2)=phi'(a)-phi'(b). Substitution in (2) proves (8). QED.

An explicit adversary is

\[
\begin{aligned}
\omega&=2\pi/\log2,\quad 0<|\varepsilon|<1,\\
\phi_\varepsilon(x)&=x\log x-
\frac{\varepsilon x}{\omega(1+\omega^2)}
[\cos(\omega\log x)+\omega\sin(\omega\log x)],\\
\phi_\varepsilon''(x)&=
\frac{1+\varepsilon\sin(\omega\log x)}{x}>0. \tag{10}
\end{aligned}
\]

It has a continuous extension at zero, is smooth and strictly convex on the
positive cone, generates exactly the same heat equation via (2), and passes
every equal binary refinement. It is not Shannon plus an affine function.

For epsilon=1/4 and (a,b)=(0.07,0.19), the binary relative defect is about
2.22e-16, while the equal ternary defect is 0.09275176 and the (0.3,0.7) split
defect is 0.08029200. These are numerical evaluations of an analytic
counterexample. Special mass ratios can accidentally hide the defect:
a/b=1/4 is one such example. A single successful pair test is insufficient.

## 5. A minimal conditional rigidity mechanism

**Theorem 2 — binary refinement plus convex transport selects Shannon.**
Assume (1), (8), phi''>0, and concavity on (0,infinity) of the diagonal
mobility f(x)=M(x,x). Then

\[
\boxed{\phi(x)=\kappa x\log x+\alpha x+\beta,\quad
M(a,b)=\Lambda(a,b)/\kappa,\quad \kappa>0.} \tag{11}
\]

**Proof.** Binary compatibility implies f(x/2)=f(x)/2. From (9), continuity
and periodicity bound f(x)/x above and below by positive constants, hence
f(0+)=0. For any positive concave f with f(0)=0, f(x)/x is nonincreasing.
But f(2x)/(2x)=f(x)/x. Monotonicity forces the ratio to be constant on every
interval [x,2x], and hence on all positive x. Write f(x)=x/kappa.
Since f=1/phi'', integration proves (11), and (2) fixes M. QED.

The diagonal-concavity assumption follows from the familiar requirement that
the local flux cost

\[
\Psi(a,b,j)=\frac{j^2}{2M(a,b)} \tag{12}
\]

be jointly convex. Indeed, restrict to a=b=x. For endpoints x1,x2 choose
j1=f(x1), j2=f(x2). Convexity at weight t gives
s^2/[2f(tx1+(1-t)x2)] <= s/2 with
s=t f(x1)+(1-t)f(x2)>0; therefore f is concave.

Conversely, the logarithmic mean is jointly concave: each integrand
a^t b^(1-t) in (4) is concave. The perspective j^2/(2z), convex and decreasing
in z>0, proves joint convexity of (12) for M=Lambda/kappa.

This produces a concrete link between binary refinement, information and
transport. Convexity of flux cost and full mobility compatibility are additional
axioms in this FIN application. They are not supplied by the finite
intertwining equation itself. This theorem does not select kappa (information
normalization), the graph W, a clock, or physical units.

**Axiom-removal checks.**

- Remove diagonal concavity: family (10) survives. At epsilon=1/4 its
  diagonal mobility satisfies f''(1)=2 epsilon^2 omega^2-epsilon omega>8.
- Remove binary compatibility: quadratic entropy phi=x^2/2 with M=1 has
  convex flux cost and exact heat matching, but 2M(a/2,b/2)=2 instead of 1.
- Remove heat matching: the positive homogeneous concave geometric mean
  M=sqrt(ab) has convex flux cost and binary compatibility without forcing
  a given entropy to generate heat.

These are necessity statements for this sufficient package, not a proof of
absolute minimality among all possible entropy axioms.

## 6. Other sufficient conditions, with explicit limits

**Theorem 3 — two multiplicatively independent subdivisions.**
Under (1), continuous compatibility with equal binary and equal ternary splits
forces (11).

**Proof.** The function h(x)=x phi''(x) is invariant under multiplication by
2 and 3. On log coordinates its periods are log 2 and log 3. Their ratio is
irrational, since an integer relation would imply 2^m=3^n. The subgroup they
generate is dense. Continuity makes h constant. QED.

Allowing only 2,4,8,... is insufficient, by Theorem 1. Earlier FIN work
considers fiber sizes 2,3,4,5, but it does not require one common entropy
geometry to satisfy (6) for all those sizes.

Arbitrary binary split fractions are another sufficient condition. Fix a,b
and set F(s)=M(sa,sb). Applying (7) to masses (s+t)(a,b), with split
s/(s+t), gives F(s+t)=F(s)+F(t). Continuity gives degree-one homogeneity.
Its diagonal again implies phi''(x)=kappa/x.

**Reference-measure caveat.** If states are represented as densities p_i/pi_i
and the reference measure is simultaneously refined as pi tensor r, every
convex relative entropy sum pi_i phi(p_i/pi_i) has a compatible divided-
difference construction on product states with matching fiber reference r.
Quadratic relative entropy is an explicit counterexample to unconditional
Shannon selection in that alternative convention. The mass-coordinate and
common-clock premises above must not be silently replaced.

## 7. A correlation defect for an observer inside a layer

Return to Shannon with kappa=1. Let P be an arbitrary strictly positive
joint state, p_i=sum_b P_ib and r_i(b)=P_ib/p_i. Define

\[
\mathscr D(P)=K^W_\Lambda(p)-
C K^{\widetilde W}_\Lambda(P)C^\top. \tag{13}
\]

**Theorem 4 — exact coarse Onsager inequality and equality case.**

\[
\mathscr D(P)\succeq0,
\qquad
\mathscr D(P)=0\ \Longleftrightarrow\ P=p\otimes r
\quad\text{when the base graph is connected.} \tag{14}
\]

**Proof.** The conductance of (13) on edge i,j equals W_ij delta_ij with

\[
\delta_{ij}=\Lambda(p_i,p_j)-\sum_b\Lambda(P_{ib},P_{jb})\ge0. \tag{15}
\]

Using (4), Holder's inequality gives, at every t in (0,1),
sum_b P_ib^t P_jb^(1-t) <= p_i^t p_j^(1-t).
Equality holds precisely when r_i=r_j. Integration preserves this strict
criterion. Thus (13) is a positive weighted graph Laplacian, and it is zero
iff r_i=r_j on every occupied base edge. Connectedness propagates equality
to all vertices. QED.

No choice of the vertical rate in B changes this horizontal statement.
Without connectedness the condition is independence on each connected
component; strict positivity avoids additional support cases.

This is a conditional diagnostic for hidden layer correlations. It is
computed from fine-state data; a coarse-only experiment cannot measure it.

## 8. Exact entropy-production accounting

**Theorem 5 — the hidden horizontal dissipation is a KL sum.**

\[
\begin{aligned}
\sum_b\mathcal I_W(P_{\cdot b})-\mathcal I_W(p)
 =\sum_{i<j}W_{ij}\big[
p_i D(r_i\Vert r_j)+p_j D(r_j\Vert r_i)\big]\ge0. \tag{16}
\end{aligned}
\]

**Proof.** On one edge substitute P_ib=p_i r_i(b) and split
log(P_ib/P_jb)=log(p_i/p_j)+log(r_i(b)/r_j(b)).
The first part sums to (p_i-p_j)log(p_i/p_j); the remaining two terms are
the indicated KL divergences. Summing edges proves the identity. QED.

This quantity vanishes under the same connected-positive independence
condition as (14). It is an entropy dissipation rate in dimensionless graph
time, not heat in joules. Full fine entropy production also includes the
nonnegative vertical-edge contribution.

At epsilon=0.6 in the explicit correlated 24-state test,

- the coarse heat identity error is 1.39e-17;
- the horizontal entropy-production excess is 0.1286853406;
- the independently computed KL expression agrees within 3.06e-16;
- the Onsager defect has eleven positive eigenvalues and one mass zero mode.

Thus an exactly closed coarse heat equation can coexist with strictly missing
fine dissipation and correlations. A uniform coarse distribution can remain
stationary throughout while this hidden quantity is positive.

## 9. What the results change in FIN

The research identifies two mathematically explicit missing requirements:

1. preservation of the full Onsager transport form under splitting, which is
   stronger than semigroup intertwining;
2. convexity of the local flux action, or an alternative such as equal binary
   plus ternary compatibility.

Together, in the stated separable mass-coordinate class, they select the
Shannon/logarithmic-mean pair up to normalization. Binary self-similarity alone
does not do so. This identifies a precise hypothesis to seek in a future
derivation of the FIN information geometry.

The correlation defect and KL identity make the internal-observer idea
quantitative: coarse dynamical closure says nothing by itself about absence
of hidden transport dissipation. They offer a fine/coarse diagnostic that
can be falsified on finite graph models without a laboratory.

Both conclusions apply to every connected positive weighted graph; they
are not a unique fingerprint of the frozen FIN parameters. The strict graph
is an explicit, nontrivial example. The historical signed legacy graph
fails the needed positivity premise; no legacy role transfer is performed.

## 10. Dual dynamics and what is still open

The self-adjoint operator still defines U_t=exp(-itA), P_t=exp(-tA) and
cos(t sqrt(A)). Unitary evolution preserves the density-matrix spectrum;
the heat channel decreases negative Shannon entropy on the vertex simplex.
The extra refinement axioms above constrain an Onsager description of the
heat channel. They neither derive the unitary measurement rule nor equate
all three evolutions.

The graph, its state preparation, the reference measure convention, dimensional
clock, continuum continuation, physical energy and instrument remain inputs.
No selector/QW-2191, Standard Model, gravity or L_total closure follows.

The general no-equivariant-section lemma in prior reports remains correct
when its vertical group is actually a symmetry of the admissible objects.
It is not a universal ban on new selection axioms. In particular arbitrary
tail translations need not preserve positivity, locality or an action; that
group must be justified before using its no-section lemma as physical evidence.

## 11. Literature comparison and originality limits

Maas's [Gradient flows of the entropy for finite Markov chains (2011)](https://arxiv.org/abs/1102.5238)
already gives the logarithmic-mean heat-gradient construction, its
divided-difference generalization, and structural assumptions such as
homogeneity and concavity. The previous FIN wording suggesting that this
construction was a newly discovered microscopic kinetic law is too strong.

Our proofs of dyadic classification, the sufficient entropy-selection
conditions and the explicit FIN refinement-defect calculation are written
out here independently. The periodic counterexample comes from a dilation
functional equation; the inequalities use Holder and the KL decomposition.
We have not established novelty against the entire literature. The claimed
contribution is a checked connection and correction within the FIN programme.

## 12. Recommended next work, in order

1. Determine whether a genuinely independent FIN action yields joint convexity
   of (12). Do not assume convexity merely to obtain Shannon.
2. Test whether its actual coarse instrument requires the full operator
   equality (6), or only the already-known heat intertwiner.
3. Use (13) and (16) to classify hidden correlations on several refinement
   levels, preserving the equality conditions and reference measure.
4. Analyze approximate binary compatibility and approximate convexity with
   explicit error bounds; exact rigidity alone does not give a robust
   empirical inference theorem.
5. Only after those steps attempt a physical preparation/measurement mapping.

The main result is a conditional entropy-rigidity theorem and an exact
correlation/dissipation diagnostic. It is not a proof that FIN generates
the observed universe.
