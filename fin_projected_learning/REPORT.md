# FIN: projected learning, self-consistency and persistent state motion

Krzysztof Żuchowski — final research report, 8–9 September 2026.  
New campaign after ST8620; **30 rounds**, ST8621–ST8650.
No PDF is generated. The prior thirty-round report remains unchanged.

## Scope and main outcome

This campaign begins with an implemented rule, not with a hypothetical
event measure. The complete `nadsoliton_neural_analysis.py`, the adaptive
sections of `The FIN Kernel as an Unknown Mat.md`, and the scoped
ST2208/ST2209 code were inspected. Their different rules and state spaces
must not be merged.

The implemented rule is a projected, noisy, leaky covariance tracker.
Its expected geometry is inherited from its supplied random-phase teacher.
The script's unconditional verification string and unattainable entropy
target are not scientific evidence of kernel emergence.

The more interesting mathematical result is not uniformly negative:
after diagonal deletion, strict can be a self-consistent **mixed-state**
kernel and can match a **pure-state time-averaged** covariance. Those facts
do not imply instantaneous pure-state self-consistency or intrinsic
selection. For the explicit closed learning extension studied below, a
dense positive nonuniform strict kernel cannot be a stationary pure-state
kernel. The proof uses an exact maximum principle, not a failed search.

All new statements are scoped to the specified rules. No Oja/PCA, classical
rate, physical entropy, selector or legacy-role claim is silently imported.

## 1. ST8621 — identify the implemented law before analysing it

Let Π delete the diagonal of a real symmetric matrix. The actual script
implements, with η=0.01 and γ=0.01,

\[
 K_{t+1}=(1-\eta\gamma)K_t+\eta\Pi(x_tx_t^T).
\]

Its explicit symmetrization is redundant in exact arithmetic because both
the previous K and the outer product are symmetric. Diagonal deletion is
not redundant. This is a linear leaky covariance filter conditional on
the input sequence; it is not the normalized single-neuron Oja/PCA update.
The latter's averaged vector equation contains Cw−(w^TCw)w, not C−γK.

In particular, the fixed-covariance unprojected flow has fixed point C/γ,
whereas the implemented projected flow has fixed point Π(C)/γ. A PSD
covariance can have an indefinite diagonal-deleted part. The ST2209
unprojected PSD no-go remains correct in its stated scope, but does not
by itself apply to the projected rule actually implemented here.

Next: calculate the input covariance instead of inferring a mechanism from
a final correlation coefficient.

## 2. ST8622 — exact mean geometry and its teacher dependence

The source draws

\[
 x_i=\cos(\Omega i+\theta)+\epsilon_i,
 \quad \Omega=\pi/4,\quad \theta\sim\mathrm{Unif}[0,2\pi],
 \quad \epsilon_i\sim N(0,1/4)
\]

independently between iterations, with independent noise coordinates.
Therefore

\[
 C_{ij}=\tfrac12\cos[\Omega(i-j)]+\tfrac14\delta_{ij}.
\]

For deterministic K0 and a=1−ηγ,

\[
 E K_T=a^T K_0+(1-a^T)\Pi(C)/\gamma.
\]

This follows directly by taking expectations of the linear recurrence.
Neither the legacy phase offset nor its hyperbolic damping appears in the
off-diagonal teacher covariance. At distance two the mean teacher entry
is zero, whereas the canonical legacy target and strict target both have
nonzero entries. No scalar Frobenius normalization fixes that discrepancy.
The phase π/4 is itself supplied to the teacher. The resulting line-indexed
cosine matrix is not a C12-circulant kernel: period eight does not become
period twelve by selecting twelve sample locations.

For twelve samples, the cosine and sine vectors have squared norm six and
are orthogonal. Thus Π(C) has eigenvalues 5/2 twice and −1/2 ten times.
The projected leaky mean is not a rank-one PCA projector.

The actual archived 30,000-step function was executed with seed 8621,
without executing its file-writing main block. The output correlation was
approximately 0.999892 with the teacher shape, 0.838287 with its legacy
line target, and 0.414173 with strict. These are reproducible numerical
checks of the exact teacher-covariance explanation, not physical evidence.

## 3. ST8623 — fixed-step fluctuations are not convergence or chaos

For i≠j, direct trigonometric and Gaussian moment calculation gives

\[
 \operatorname{Var}(x_ix_j)=\tfrac18+\sigma^2+\sigma^4=7/16,
 \qquad \sigma^2=1/4.
\]

Indeed E[cos²u cos²(u+δ)]=1/8+(cos²δ)/4; subtracting the squared
mean cancels the phase-dependent term. Independent noise contributes
σ²+σ⁴. Conditional on deterministic initialization,

\[
 \operatorname{Var}(K_{T,ij})
 =\eta^2(7/16)\frac{1-a^{2T}}{1-a^2}.
\]

For fixed nonzero η this tends to a positive stationary variance, not zero.
The source's 30,000 iterations also retain a^30000≈0.04978 of the initial
condition in the exact linear recurrence. A stationary stochastic filter
must not be called an exact deterministic convergent kernel merely because
one shape correlation is high.

Conversely, persistent noise does not demonstrate chaos. Two kernels driven
by the same input stream obey ΔK_T=a^T ΔK0 exactly, with conditional
Lyapunov exponent log|a|<0. The physical source of the input and the
iteration clock remain supplied parameters.

## 4. ST8624 — an impossible entropy target and a false verification label

The script normalizes twelve softmax eigenweights to a probability vector,
then compares its Shannon entropy with 4 log2=log16. Every probability
vector on twelve labels satisfies

\[
 H\le\log12<\log16.
\]

The minimum possible relative discrepancy is at least
1−log12/log16≈0.1037594. This is a carrier-size bound, not a numerical
failure to optimize. The four-bit identity belongs to a different count
of alternatives unless another state space or measure is explicitly supplied.
The target-dependent Frobenius normalization also changes the softmax
temperature scale; entropy is not a scale-free geometry statistic.

An adversarial call to the actual analysis function with the negative of
its own target returns correlation −1 while still emitting the string
“SIMULATION VERIFIED: Nadsoliton geometry emerges from Hebbian learning.”
The string is unconditional. This refutes its use as a validation gate.
The historical source is preserved; its output is reclassified here.

## 5. ST8625 — dephasing retains degenerate blocks

For fixed Hermitian K and a state ρ0, the exact infinite-time average is

\[
 \overline\rho=\sum_\lambda P_\lambda\rho_0P_\lambda.
\]

For real K, taking real parts commutes with this expression. The proof is
the finite integral of exp[it(λ−μ)]: unequal frequencies average to zero,
equal ones do not. A degenerate block is not generally a scalar multiple
of its full projector. The two-dimensional example K=0, ρ0=diag(1,0)
stays diag(1,0); replacing it by its total occupation times I gives I,
which even has the wrong trace.

The later ST2205 implementation groups degenerate eigenspaces correctly.
This audit preserves that correction and rejects only the overly broad
scalarized reading of the older equation. Seven strict Fourier-sector
eigenvalue intervals were recomputed with exact rational enclosures and
shown pairwise disjoint. Their multiplicities are 1,2,2,2,2,2,1.

## 6. ST8626 — projected mixed-state self-consistency exists, but does not select strict

Consider the explicitly declared density-state extension

\[
 \dot\rho=-i[K,\rho],\qquad
 \dot K=\eta[\Pi(\operatorname{Re}\rho)-\gamma K],
\quad K\in\mathrm{Sym}_0,\quad \rho\succeq0,\quad\operatorname{Tr}\rho=1.
\]

Throughout this closed extension η>0 and γ>0.

This is a closed candidate extension of the repository's state/operator
rule, not the externally driven process in rounds 1–4. The conventional
Hamiltonian sign is immaterial to the stationary and Lyapunov results;
the older example uses the opposite sign. Allowing a mixed state is an
explicit enlargement of its pure-vector state class.

**Fixed-family theorem [Proven].** For any real symmetric zero-diagonal K
such that ρ_K=I/n+γK is PSD, (K,ρ_K) is an exact stationary pair.
Its trace is one, it commutes with K, and Π(ρ_K)=γK.
For nonzero K, trace zero ensures λ_min(K)<0, and the positive-γ range is
0<γ≤−1/[n λ_min(K)].

For strict and γ=1/20 the minimum density eigenvalue is approximately
0.04923960, with a positive exact rational lower enclosure in the saved
certificate. Both derivatives vanish analytically and in the numerical
replay. Thus “strict can never be a projected self-consistent covariance”
is false. The unprojected ST2209 theorem is not refuted.

This does not constitute emergence of strict: the same construction works
for an open continuum of K matrices and encodes K in the initial density.
Existence of a stationary pair is not selection of that pair.

## 7. ST8627 — a pure time-average witness does not give instantaneous stationarity

Strict has one-dimensional constant/Nyquist sectors and five real
two-dimensional cosine/sine sectors. Let c_λ=1/n+γλ>0. On a two-dimensional
sector with orthonormal vectors u_λ,v_λ, choose the component
sqrt(c_λ)(u_λ+i v_λ); on each one-dimensional sector choose
sqrt(c_λ) times its unit vector. The resulting pure vector ψ has norm one.

Unequal-sector cross terms dephase, while within each two-dimensional block

\[
 \operatorname{Re}[(u+iv)(u+iv)^*]=uu^T+vv^T.
\]

Therefore its real time-averaged covariance is exactly I/n+γK. Its
off-diagonal average matches γK. The strict test has error below 10^-13.
This establishes a pure-state *averaged* witness without changing the
unprojected positivity theorem.

At an instant the real pure covariance has rank at most two. It need not
equal the full-rank average. The constructed strict witness has an
instantaneous off-diagonal mismatch about 0.360062, so K does not stay fixed
under the instantaneous pure-state update. Averaging over a fixed K and
evolving a self-changing K are not interchangeable operations. An actual
slow-learning averaging theorem would require additional analysis.

## 8. ST8628 — a Lyapunov law and conserved state spectra

For the declared coupled density model set

\[
 F(K,\rho)=\frac\gamma2\|K\|_F^2-\operatorname{Tr}(K\rho).
\]

The state term contributes zero along the commutator flow, by cyclicity of
trace. Projection is orthogonal on Sym0. Hence

\[
 \dot F=-\eta\|\Pi(\operatorname{Re}\rho)-\gamma K\|_F^2\le0.
\]

The squared-completion bound gives F≥−||ρ||_F²/(2γ)≥−1/(2γ).
The state evolves by unitary conjugation even when K varies, so every
Tr(ρ^m), its spectrum and von Neumann entropy remain constant. Kernel
learning is not automatically learning or erasing the state's eigenvalue
distribution.

The kernel's variation-of-constants formula bounds
||K(t)||_F by exp(−ηγt)||K0||_F+(1−exp(−ηγt))/γ.
Together with compactness of the density-state orbit and local smoothness,
this excludes finite-time blowup. Limit-set reasoning can restrict motion
to the invariant part of the zero-learning set. It does not yet prove
convergence to one stationary state or selection of strict.

## 9. ST8629 — zero learning can coexist with persistent state motion

For n=2 take K=kσ_x and any positive density with real off-diagonal entry
γk. Rotation about the x axis leaves that real coherence invariant while
rotating its y/z Bloch components. Thus Π(Reρ(t))=γK and Kdot=0 for all
t, although ρdot need not vanish.

The explicit test k=0.7, γ=0.2 and
ρ0=[[0.7,0.14],[0.14,0.3]] has positive eigenvalues. Kernel learning is
exactly zero on its orbit, while the population oscillates and
||ρdot||_F≈0.395980. This falsifies the inference “a decreasing learning
functional forces the whole state to become stationary.” It is a finite
matrix oscillation, not a physical soliton or a derivation of time units.

Next: determine whether the same persistent zero-learning mechanism can
support the actual dense positive strict matrix.

## 10. ST8630 — dense-positive rigidity and a pure-state obstruction

**Maximum principle [Proven].** If n≥3, K has zero diagonal and all
off-diagonal entries are positive, and D=diag(d_i) satisfies
Π([K,[K,D]])=0, then D is scalar.

For i≠j the equation is

\[
 \sum_{l\ne i,j}K_{il}K_{lj}(d_i+d_j-2d_l)=0.
\]

Choose indices of the largest and second-largest d. If the largest is
strictly larger, every summand is positive, a contradiction. If the two
are equal, equality in the positive weighted average forces every other
d to have that same value. This proves scalarity without a numerical rank
test. The strict finite matrix of these constraints has rank 11 as a check.

**Invariant-set consequence.** On a trajectory with Kdot=0 and constant
dense-positive K, write ρ=X+iY, with X real symmetric and Y real skew.
Then X=γK+D(t), Xdot=[K,Y], Ydot=−[K,D]. Differentiating X shows that
[K,[K,D]] must be diagonal. The maximum principle gives D=I/n, fixed by
trace. Thus Xdot=0, Ydot=0 and [K,ρ]=0: the state is stationary. The
two-dimensional counterexample lies outside the n≥3 dense-positive
premise, as required.

**Pure-state corollary [Proven].** Such a stationary pair with rank-one ρ
has the uniform complete-graph kernel K_ij=1/(γn) for i≠j. In particular,
the nonuniform strict W is impossible as an instantaneous pure-state
stationary kernel in this declared model.

To prove this, ρ=ψψ* commuting with K makes ψ a K eigenvector. Because K
is real, Reψ and Imψ lie in the same eigenspace. The range of
Reρ=γK+I/n is therefore contained in that eigenspace. The positive Perron
vector of K has positive eigenvalue under γK+I/n, so it belongs to the
range. Perron simplicity forces ψ to be that real vector up to a global
phase. The diagonal Reρ=1/n then forces its entries all equal, and the
off-diagonal stationary equation gives the stated constant K. A stationary
density permits a physically irrelevant global phase rotation of ψ.

This is a new, correctly scoped obstruction for the **projected pure** law;
it coexists with the mixed-state fixed family and the pure averaged witness.
It does not exclude other nonlinear, delayed, stochastic or constrained
learning laws. It does not show that an arbitrary evolving K remains in
the positive cone, nor that every orbit approaches a dense-positive limit.

## Outcome after round 10

The actual source learns its teacher covariance; its advertised entropy
target and unconditional pass label do not establish FIN emergence. A
closed projected extension has a real Lyapunov structure, but its
self-consistency depends crucially on pure versus mixed state and on
instantaneous versus averaged updating. The next high-information questions
are minimum covariance-completion rank, stability of these fixed families,
and whether the positive strict cone is preserved by a physically meaningful
learning law. Those questions lead to the following investigations.

## 11. ST8631 — a diagonal-independent lower bound on covariance rank

For a pure complex vector ψ=a+ib,
Re(ψψ*)=aa^T+bb^T has real rank at most two. More generally a density
matrix of complex rank r has real-part rank at most 2r. Diagonal deletion
can greatly increase a matrix's rank, so rank(K) alone is not a valid
rank test for its covariance source.

Use instead disjoint row and column sets I,J. If Π(Reρ)=γK with γ≠0,
then (Reρ)_(I,J)=γK_(I,J), entirely independent of all missing diagonal
entries. A nonzero 6×6 disjoint minor therefore forces

\[
 \operatorname{rank}(\operatorname{Re}\rho)\ge6,
 \qquad\operatorname{rank}_{\mathbb C}\rho\ge3.
\]

**Exact finite witnesses [Proven].** Such minors exist for strict, canonical
cyclic legacy and canonical line-indexed legacy. The selected index sets
and rational determinant intervals are retained in `geometry_results.json`.
The strict partition is I=(0,1,4,5,8,9), J=(2,3,6,7,10,11); its
determinant is enclosed strictly below zero. Legacy has independently
certified nonzero minors as well.

Numerical singular values selected the candidate minors, but did not certify
them. Certification uses outward rational weight intervals and the complete
720-term determinant expansion with integer interval products. Consequently
no choice of diagonal or amplitude rescales either declared kernel into
the instantaneous real covariance of one pure vector. This does not exclude
time averaging or an ensemble, and does not identify the minimum feasible
rank as exactly three.

## 12. ST8632 — exact state-rank classification at a strict stationary kernel

At a stationary pair in the closed model, write ρ=X+iY with real symmetric
X and real skew Y. The equations imply ΠX=γK and [K,X]=[K,Y]=0.
Thus X=γK+D with D diagonal and [K,D]=0. On a connected support D must
be scalar; trace fixes D=I/n. Positivity of the edge weights is not needed
for this stationary diagonal argument, only connected support.

For strict, let λ0,…,λ6 denote the seven Fourier-sector eigenvalues,
with λ6 the unique lowest one. The earlier exact spectral intervals certify
this ordering and the multiplicities. Set c_j=1/12+γλ_j. Every stationary
density has singleton blocks c0,c6, and on the kth real doublet a block

\[
 \rho_k=c_k I_2+i b_k\begin{pmatrix}0&1\\-1&0\end{pmatrix},
 \qquad |b_k|\le c_k,\qquad k=1,\ldots,5.
\]

These formulas are exhaustive: commutation forbids cross-eigenvalue blocks,
the required real part is scalar on each doublet, and a real skew 2×2
matrix has one parameter. Positivity gives eigenvalues c_k±b_k.

Let γ*=−1/(12λ6)≈0.1222120805. Stationary densities are possible for
0<γ≤γ*. For 0<γ<γ*, the minimum complex rank is **seven**: one in each
positive doublet and two singleton blocks. At γ=γ*, c6=0 while all
other c_j>0, and the minimum rank is exactly **six**. Both bounds are
attained by choosing |b_k|=c_k in every doublet.

Thus six is the least rank over this allowed γ range; it is not the
minimum at an arbitrary fixed decay parameter such as γ=0.05. A density
of rank r needs at least r pure components in any ensemble decomposition.
Equivalently a purification needs auxiliary Hilbert dimension at least r,
as follows from the rank of its coefficient matrix; spectral purification
attains that bound. These are state-representation resource counts, not
quark flavours, spatial dimensions or proved physical constituents.

## 13. ST8633 — minimum rank or entropy leaves 32 phase choices

At fixed admissible γ, every rank-minimizing stationary state has
b_k=±c_k independently, giving exactly 2^5=32 states in the specified
vertex/eigenprojector setting. Reflection sends all b_k to −b_k, pairing
the states. Rotations leave each circulant stationary state invariant.
Their real covariance and the learned K are identical.

The entropy contribution of a doublet is
−(c+b)log(c+b)−(c−b)log(c−b), and its derivative for b>0 is
−log[(c+b)/(c−b)]<0. Hence minimum von Neumann entropy at fixed K,γ
also occurs at |b_k|=c_k; maximum entropy occurs at b=0. Purity is
likewise maximized at the vertices. These criteria do not select the signs.

They are not mechanisms of the current dynamics: the state spectrum,
entropy and purity are conserved. A generic full-rank initial state cannot
evolve into a lower-rank state under its unitary state equation. Nor are
all 32 states automatically distinct physical worlds: distinguishing them
requires an admitted phase-sensitive preparation or measurement. No
selector or physical chirality source is exported by this static criterion.

## 14. ST8634 — exact global minimum of the learning functional

For fixed state spectrum, let P=Trρ². Complete the square in K and use
ρ=X+iY to obtain the exact identity

\[
 F+\frac{P-1/n}{2\gamma}
 =\frac\gamma2\left\|K-\frac{\Pi X}{\gamma}\right\|_F^2
 +\frac{\|Y\|_F^2+\|\operatorname{diag}X-\mathbf1/n\|_2^2}{2\gamma}.
\]

Therefore the global minimum is at least −(P−1/n)/(2γ), with equality
precisely when ρ is real, its diagonal is uniform, and K=(ρ−I/n)/γ.
The bound is attainable on every density-matrix isospectral orbit.

**Attainability proof.** Start from a real diagonal matrix with the specified
eigenvalues, and subtract I/n. Any real symmetric trace-zero matrix admits
an orthonormal basis of zero-Rayleigh vectors: its smallest and largest
eigenvalues bracket zero, so a vector with Rayleigh value zero exists.
The compression to its orthogonal complement is again trace zero. Induct.
In that basis the density is real with diagonal 1/n, as required.

Thus the strict mixed fixed pair from round 6 is a global minimum on its
own encoded state-spectrum orbit, not merely an accidental stationary
point. This is a conditional positive result. It does not select the input
spectrum or the eigenbasis. The set of all minima is Lyapunov stable on
the fixed-spectrum state space: completing the square controls departure
from its compact zero set in any bounded kernel region. This does not
prove asymptotic stability of an individual strict point or convergence of
every trajectory to the minimum set.

## 15. ST8635 — a 50-dimensional local manifold of equally good minima

For strict ρ0=I/12+γK, the real orthogonal isospectral orbit has dimension
66−5=61: O(12) has dimension 66 and its five doublet stabilizers contribute
five continuous dimensions. The uniform-diagonal condition has local
rank eleven. Indeed its derivative on a skew matrix Z is
diag([Z,ρ0]), whose ith entry is 2 sum_j Z_ij(ρ0)_ij.
This is a weighted connected-graph incidence map onto trace-zero diagonal
vectors. All strict off-diagonal weights are positive, so its rank is eleven.

The regular-level-set theorem therefore gives a smooth **50-dimensional**
local minimum manifold. Strict positivity of kernel edges persists in a
sufficiently small neighborhood. Every point has the same state spectrum
and the same minimal F, but need not have the same kernel geometry.
The numerical tangent ranks 61 and 11 check these analytic dimensions;
they do not replace the argument.

An individual strict equilibrium is consequently not asymptotically stable
against every nearby state on that isospectral space: arbitrarily close
distinct equilibrium points never converge to it. This conclusion is about
the supplied unconstrained Sym0 learning class.

## 16. ST8636 — 39 dimensions remain even after preserving regularity and spectrum

Also require that the uniform vector remain the Perron eigenvector, so K
keeps its strict constant row sum. The real orbit is now generated by O(11)
on the uniform-orthogonal subspace and has dimension 55−5=50. The
uniform-diagonal map still has rank eleven, giving a **39-dimensional**
local minimum manifold with the same spectrum and row sum as strict.

The additional rank is proved rather than inferred from singular values.
A diagonal matrix D annihilates that derivative exactly when
P_perp[K,D]P_perp=0. In Fourier coordinates a nonconstant diagonal Fourier
mode q has entries proportional to λ_k−λ_(k−q), excluding indices zero.
There are ten allowed k for each q≠0. Equality of two strict sector
eigenvalues requires their indices to be opposite, so only solutions of
2k=q can cancel such a difference; there are at most two. Thus every
nonconstant Fourier coefficient of D vanishes. The annihilator is just the
constant diagonal, proving rank eleven.

A small numerically integrated constrained isospectral curve supplies a
representative: row, diagonal and spectral errors are below 10^-13,
all edges stay positive, and an edge value moves away from the original
six values. Its learning derivatives vanish. This curve **parametrizes
different equilibria**; it is not a moving solution of the learning law.
The analytic manifold theorem licenses existence; the integration is a
witness check. Since vertex permutations form a finite set, the local
continuum cannot be explained entirely by relabeling.

## 17. ST8637 — imposing circulant symmetry drastically narrows the result

The previous manifolds generally leave the circulant class. If one instead
requires a real nonnegative circulant zero-diagonal K with the exact strict
spectrum, the constant Fourier mode must carry the simple Perron eigenvalue,
and the Nyquist mode must carry the other simple eigenvalue. The five
distinct doublet eigenvalues may only be assigned in 5!=120 ways.

For every assignment the inverse Fourier weights were enclosed with exact
rational intervals. Exactly two assignments have all positive weights;
each of the other 118 has a certified strictly negative edge. No sign was
left unresolved. The two positive kernels are

\[
 K\quad\text{and}\quad P_5 K P_5^T,\qquad i\mapsto5i\pmod{12}.
\]

Their distance-one and distance-five weights are exchanged. They are one
orbit under group automorphisms of C12. Therefore this restricted inverse
problem determines the kernel up to that automorphism. This is a real
conditional identification result, not a general “spectrum never helps”
claim. It still assumes the full strict spectrum and the exact circulant
carrier class. A distinguished nearest-neighbor metric can make the two
labeled assignments consequential; otherwise their relabeling is a gauge.

Equivariance of a learning law is weaker than requiring every state to be
circulant. Symmetric initial data remain in a symmetry-fixed set, but
generic initial data of the actual Sym0 rule need not. Neither the full
spectrum nor that stronger state restriction is sourced by this census.

## 18. ST8638 — the raw learning law can leave the strict positive cone

At a zero edge K_ij=0, a valid pure state with Reρ_ij<0 gives Kdot_ij<0.
There is no general cone-invariance theorem for the raw projected-Hebb rule.
The failure can start from the actual strictly positive strict matrix,
not merely from an artificial boundary state.

Take ψ=(e0−e1)/sqrt2, γ=1/20, η=1000, and fast time θ=ηt. If the
state were frozen, its edge at θ=1 would be

\[
 K_{01}^{\rm frozen}=e^{-1/20}(K_{01}(0)+10)-10
 <-3233/80000,
\]

using K01<0.47 and e^(-1/20)≤1−1/20+(1/20)²/2=761/800.
The exact spectral enclosures give ||K0||_F<6. On 0≤θ≤1,
||K(θ)||_F≤6+θ and the unitary state equation gives
||ρ(θ)−ρ0||_F≤(12θ+θ²)/η. Variation of constants then bounds the
kernel's deviation from the frozen-state solution by
(6θ²+θ³/3)/η. At θ=1,

\[
 K_{01}(t=0.001)<-3233/80000+19/3000=-8179/240000<0.
\]

This is an analytic counterexample with a numerical integration check
(approximately −0.04064154). It does not assert exit at every learning
rate or every state. Negative learned edges do not define a classical
Markov generator merely because the kernel remains Hermitian.

## 19. ST8639 — a positivity repair is a different learning law

One can replace the raw kernel update by its orthogonal projection onto the
tangent cone of nonnegative symmetric zero-diagonal matrices. An interior
edge retains its old update; a zero edge has update
η max(Reρ_ij,0). For this declared constrained law, positive-cone viability
and the projected-gradient energy identity hold:

\[
 \dot F=-\eta\|\operatorname{Proj}_{T_K}
       (\Pi\operatorname{Re}\rho-\gamma K)\|_F^2\le0.
\]

The identity follows entrywise; a clipped outward component contributes
zero, while every retained component contributes minus its square. Kernel
norm bounds remain unchanged because zero edges contribute nothing to the
norm derivative. A construction via the scalar reflecting-map solution
for each exponentially rescaled edge, coupled locally to the smooth state
equation, supplies the projected dynamics; its Lipschitz reflection map
gives local uniqueness and the bounds permit continuation.

This is not a repair derived from the old source. It is an explicit change
of admissible dynamics. At a dense positive equilibrium it agrees with the
old law, so the pure-state strict obstruction still holds.

New boundary equilibria are possible: a pure vector with six components
+1/sqrt12 and six components −1/sqrt12, together with
K=max(ΠReρ,0)/γ, gives two disconnected positive six-cliques and zero
state/kernel derivatives. This proves existence, not their stability or
selection. They are not the connected strict kernel or physical particles.

## 20. ST8640 — adaptive duality requires more than positive weights

The implemented closed state's Hamiltonian is K (up to sign convention),
while the natural classical graph generator would use
A(K)=diag(K1)−K. The commutator [K,A(K)] has entries
(s_j−s_i)K_ij, with s_i=sum_j K_ij. Hence they commute exactly when
degree is constant on each connected support component. For connected
strict-like positive support this means a single constant row sum.

The raw learning rule can break regularity before any edge becomes negative.
For the same pure initial state used in round 18,
s'_0(0)−s'_2(0)=−η/2. Since K02(0)>0, continuity and this nonzero
derivative prove loss of commutation for sufficiently small positive time,
while all edges are still positive.
In the tested early segment all edges remain positive, but the degree
spread is approximately 0.0498752 and ||[K,A(K)]||_F≈0.0558889.
Thus Hermiticity plus entrywise positivity alone does not preserve the
shared instantaneous spectral projectors of the frozen strict duality.

Choosing A(K) as the state Hamiltonian instead would change the supplied
equation and generally invalidate the old F calculation; it cannot be a
silent substitution. Even if regularity is enforced at each time, matrices
at different times may not commute, so time-ordered propagation is not
functional calculus of one fixed A. The frozen identities U_t=e^-itA and
P_t=e^-tA remain correct. Their automatic extension to an arbitrary
adaptive learning trajectory is what is rejected.

## Outcome after round 20

The model admits exact mixed self-consistency and large global-minimum
families, but not unique selection of strict. Imposing the full circulant
spectral class gives a meaningful conditional uniqueness up to an
automorphism; that restriction must not be mistaken for a derived source.
State rank, initial spectrum, positivity, regularity and temporal
composition are distinct obligations. The remaining investigations test
whether a more complete dynamical or variational prescription resolves them.

## 21. ST8641 — strict can be learned exactly from a state that already encodes it

Set ρ0=I/n+γW and K0=0 in the original closed density-learning model.
Then the exact trajectory is

\[
 \rho(t)=\rho_0,\qquad K(t)=(1-e^{-\eta\gamma t})W.
\]

Every K(t) commutes with ρ0, so the state equation vanishes. Substituting
in the kernel equation verifies the remaining derivative. Positivity,
regularity and common spectral projectors are preserved on this particular
trajectory. This is a removal test for an overly broad interpretation of
rounds 18–20: they did not exclude every successful trajectory.

However, the state already contains W in its off-diagonal entries. This
construction demonstrates reproduction, not an intrinsic selection of
strict from less specific data. By contrast, the maximally mixed state I/n
supplies no off-diagonal drive and makes K decay. A uniform pure amplitude
has the same vertex probabilities but supplies a uniform-clique drive.
Population values alone do not specify the coherence resource being learned.

## 22. ST8642 — a dual-compatible Dirichlet gradient has its own obstruction

To use the graph Laplacian A(K) as the state Hamiltonian, consider the
explicit candidate functionals

\[
 F_\pm(K,\rho)=\frac\gamma2\|K\|_F^2
                  \pm\operatorname{Tr}(A(K)\rho),
 \qquad \dot\rho=-i[A(K),\rho].
\]

Let p_i=ρ_ii and

\[
 q_{ij}=\frac{p_i+p_j}{2}-\operatorname{Re}\rho_{ij},\qquad q_{ii}=0.
\]

Positivity of ρ implies q_ij≥0, by the 2×2 minor inequality
Reρ_ij≤sqrt(p_i p_j)≤(p_i+p_j)/2. In the symmetric zero-diagonal
Frobenius metric, the gradient of Tr(A(K)ρ) is q. The Hamiltonian state
part preserves that expectation at fixed K.

Thus descent of F_+ has Kdot=−η(q+γK): it cannot generate positive
edges, and under nonnegative-cone projection every edge is bounded above
by its exponentially decaying initial value. Descent of F_- instead gives
Kdot=η(q−γK), which preserves nonnegative edges by positive forcing.
It supplies a coherent conditional joint propagation/learning model.

Here compatibility means that both channels can use the same instantaneous
Laplacian. It does not equate unitary populations with heat probabilities,
or turn a noncommuting adaptive path into functional calculus of one fixed
matrix.

**Strict stationary obstruction [Proven].** For every real zero-sum v,

\[
 v^Tqv=-v^T\operatorname{Re}\rho\,v\le0.
\]

So q is conditionally negative semidefinite. Stationarity of the active
F_- law would require q=γK. Strict W has a positive Fourier eigenvalue
on the zero-sum subspace (the k=1 value is approximately 0.9061861246,
with a positive exact enclosure). Therefore γW cannot equal such a q.

This no-go holds for arbitrary density rank in this particular linear
Dirichlet-gradient law. It does not exclude other potentials or constraints.
It also shows why changing “correlation” to “distance” in a learning source
is a substantive change, not a harmless interpretation of the same formula.

## 23. ST8643 — a regular positive projected law is viable but still nonselecting

Supply the compact convex set

\[
 \mathcal C_s=\{K=K^T,\ K_{ii}=0,\ K_{ij}\ge0,\ K\mathbf1=s\mathbf1\}.
\]

It contains strict and the uniform positive matrix with edge weight
s/(n−1). On this set A(K)=sI−K, so the old F from round 8 is preserved
by the state step with Hamiltonian A(K). Follow that unitary step by

\[
 K^+=\operatorname{Proj}_{\mathcal C_s}
          [K+h(\Pi\operatorname{Re}\rho-\gamma K)].
\]

The projection is unique. Its optimality inequality and the exact quadratic
Taylor expansion give

\[
 F(K^+,\rho)-F(K,\rho)
 \le-\left(\frac1h-\frac\gamma2\right)\|K^+-K\|_F^2,
 \quad 0<h<2/\gamma.
\]

Thus this fully specified split update preserves density positivity and
unitarity at the state step, kernel nonnegativity and regularity, and a
decreasing learning functional. An interior step was evaluated using the
exact affine row-sum projection; it remains strictly positive and satisfies
the inequality. Boundary projection would require the full convex problem,
not clipping edges after projecting row sums.

The 39-dimensional family from round 16 remains fixed under this law.
Consequently this repair does not select strict. It also supplies s as a
coupling budget: K=0 is not an admissible initial point when s>0. This is
a viable conditional construction, not a derivation of that budget or the
constraint from primordial information.

**Further scope test.** Replacing Π by Π composed with cyclic spatial
averaging gives a different, also explicit escape from the raw pure-state
obstruction. Prepare the pure Fourier witness of round 7 and start K0=0.
Its modal occupations are fixed along the exact trajectory
K(t)=(1−exp(−ηγt))W; cyclic averaging of Reρ(t) is I/n+γW at every
instant. This verifies the changed learning equation directly, without
appealing to a slow-time averaging approximation. The raw unaveraged source
still differs. Thus the pure-state rank/no-go statements must not be applied
to this averaged source map. The projection and the target-encoding modal
occupations are supplied, so this is again reproduction rather than an
intrinsic selection theorem.

## 24. ST8644 — normal ordering must be undone before learning a precision

The older precision proposal uses an inverse of the normal-ordered strict
kernel. W is indefinite, so a positive multiple of W^-1 minus a nonnegative
mass shift cannot be a positive precision matrix: its negative eigenvalues
remain negative. A positive determinant alone is not sufficient either.
The log-determinant loss below is defined on positive-definite matrices,
or on an explicitly fixed positive support; changing support is another
problem.

There is a constructive alternative. Let G=cI+W with c>s. Then G is
positive definite and its inverse has the Neumann expansion

\[
 G^{-1}=c^{-1}I-c^{-2}W+
             \sum_{k\ge2}(-1)^k c^{-k-1}W^k.
\]

Since (W^k)_ij≤w_max s^(k−1), the absolute off-diagonal remainder is
at most w_max s/[c²(c−s)]. Thus every off-diagonal entry of G^-1 is
negative if w_max s/(c−s)<w_min. The strict enclosures permit c=100:
s<5/3, w_max<0.47 and w_min>0.011 give an upper ratio 47/5900<0.011.

Regularity gives G^-1 1=(c+s)^-1 1. Therefore

\[
 m^2=(c+s)^{-1},\qquad L=G^{-1}-m^2I
\]

is a positive graph Laplacian with strictly positive edge conductances,
and (L+m²I)^-1=cI+W exactly. Diagonal deletion recovers W. This is an
exact screened **weighted-graph** Green construction, not an exact claim
about the original nearest-neighbor cycle Yukawa fit. Infinitely many c
are allowed, giving different parents and dimensionless mass parameters.

These changes are not generally one clock rescaling. On a strict A mode
with eigenvalue a>0, the parent L eigenvalue is
a/[(c+s)(c+s−a)], so two distinct modes have c-dependent relative rates.
The factorization does not identify the physical Hamiltonian or clock.

The raw signed legacy kernel cannot be the off-diagonal part of an
ordinary positive screened Markov-graph Green function, which is entrywise
positive on connected support. A signed or lifted operator can have a
different interpretation, as earlier reports establish. No legacy role
transfer or unique physical mass is inferred from this construction.

## 25. ST8645 — the quoted precision flow is ascent, not descent

For a fixed C≻0 on the positive-definite domain,

\[
 \mathcal F(L)=\operatorname{Tr}(LC)-\log\det L,
 \quad\nabla\mathcal F=C-L^{-1}.
\]

The gradient expression in the kernel report is correct, but its displayed
evolution Ldot=C−L^-1 increases this F:
Fdot=||C−L^-1||². If the intended procedure is free-energy minimization,
the dynamical sign is reversed. In the scalar case C=1, L0=1/2, the
displayed ascent reaches the positive-domain boundary in finite time
log2−1/2; this follows by integrating dt=L/(L−1)dL.

The corrected descent is Ldot=L^-1−C. It has the unique equilibrium
L=C^-1. It stays in a compact positive-definite sublevel set: if C≥aI,
then F(L)≥sum_j[aλ_j−logλ_j], which diverges at a zero or infinite
eigenvalue. Strict convexity and the gradient identity yield global
convergence. On a bounded sublevel set with λ_max(L)≤M, the Hessian
satisfies D²F[H,H]≥M^-2||H||², giving a standard exponential contraction
to the optimum in that set.

The noncommuting two-dimensional test converges to C^-1 with final error
about 2.62×10^-12. This corrects a procedure; it does not derive C or a
physical propagation law from the kernel.

## 26. ST8646 — self-consistency requires the full chain rule or a joint functional

If C=C(L), the total derivative of Tr[L C(L)]−logdetL includes
(DC(L))^*[L]. Treating C as frozen is a different algorithm. For
C(L)=L^-1, frozen-covariance descent gives zero at every L, while the
composed functional is n−logdetL and has gradient −L^-1. These are
not the same variational problem.

A valid joint covariance/precision functional is

\[
 \mathcal J(L,C)=\operatorname{Tr}(LC)-\log\det L-\log\det C-n\ge0.
\]

Its nonnegativity follows by diagonalizing the positive matrix
L^(1/2) C L^(1/2) and using x−log x−1≥0. Equality holds exactly when
C=L^-1. It is the dimensionless Gaussian covariance-divergence expression,
not a unit-bearing physical action. It supplies a genuine variational
bootstrap, but every reciprocal positive pair is a minimizer. Minimizing
over C makes the envelope identically zero; the missing chain term is
then handled correctly, without producing a unique L.

This positive construction and its flat family distinguish an inconsistent
gradient claim from the separate, unresolved problem of selecting a source.

## 27. ST8647 — simple Gibbs bootstrap cannot support seven strict sectors

Suppose, as an additional candidate closure, that a stationary learning
state is the Gibbs function ρ=exp(βK)/Z and also satisfies Πρ=γK,
with γ>0 and connected support. Commutation gives uniform diagonal as
in round 12, so each kernel eigenvalue satisfies

\[
 e^{\beta\lambda}/Z=\gamma\lambda+1/n.
\]

For β≠0 the difference of the two sides is strictly convex, so it has
at most two real zeros. For β=0 the nonzero kernel is excluded directly.
Strict has seven distinct eigenvalues, certified in the earlier checkpoint,
and cannot satisfy this closure. The analogous positive-precision equation
L^-1=exp(−βL)/Z has λexp(−βλ)=constant and at most two positive roots
for β>0. It likewise cannot reproduce the seven-sector lifted strict
precision in round 24.

This is a no-go for these specified scalar Gibbs/self-consistency equations,
not for all thermodynamic or nonlinear learning mechanisms. There is no
uniform positive residual gap as γ→0: choosing β=nγ makes the residual
order γ². The numerical near-zero examples are explicitly recorded to
prevent a floating tolerance from being mistaken for an exact fixed point.

For a regular frozen K, exp(−βA)/Tr exp(−βA)=exp(βK)/Z. Therefore
“favor high K” and “favor low A=sI−K” do not by themselves define different
principles. Neither makes these dimensionless spectral quantities physical
energies without an additional interpretation.

## 28. ST8648 — infinitely fast learning is not infinitely fast selection

For the original closed model, set T(ρ)=ΠReρ/γ and E=K−T(ρ).
The exact lag equation is

\[
 \dot E=-\eta\gamma E-\dot T(\rho).
\]

With M=max(||K0||_F,1/γ), the earlier bounds give
||dot T||≤2M/γ. Consequently

\[
 \|E(t)\|\le e^{-\eta\gamma t}\|E(0)\|
                 +\frac{2M}{\eta\gamma^2}.
\]

At fixed γ, the finite-time limit η→∞ has
K=T(ρ) and the nonlinear unitary state equation
ρdot=−i[T(ρ),ρ]. Its functional
E_red=||ΠReρ||²/(2γ) is conserved, since its state gradient is T(ρ)
and the trace pairing with its own commutator vanishes. The reduced vector
field is Lipschitz on the bounded density set, so the lag estimate and
Gronwall's inequality justify the finite-time state limit.

If E(0)=0, the learning dissipation over a fixed interval [0,T] is at most
4M²T/(ηγ²), tending to zero. For off-manifold initialization there may
instead be a finite initial-layer loss; it is not covered by that vanishing
bound. Numerical tests with increasing η confirm the lag, state and
energy estimates. Faster tracking and stronger long-time selection are
not the same property.

## 29. ST8649 — the diagonal algebra is part of the model, not spectral data

Let Δ keep the diagonal. If a unitary U satisfies
Δ(UXU*)=UΔ(X)U* for every Hermitian X, apply it to a diagonal rank-one
projector E_i. Then UE_iU* must itself be diagonal, rank one and trace one,
so it equals some E_(π(i)). Thus U is a permutation times diagonal phases;
the converse is immediate. The same normalizer controls Π=I−Δ.

Requiring every real symmetric kernel to stay real restricts relative
phases to signs. Preserving every connected strictly positive off-diagonal
kernel then permits only a common sign, leaving vertex permutations and a
global phase. An arbitrary unitary change of eigenbasis is not a symmetry
of the unchanged diagonal-learning rule. A Hadamard example explicitly
violates covariance of Π.

This clarifies the 39-dimensional minimum family: it cannot be dismissed
as a gauge of the unchanged vertex algebra merely because its members
are abstractly isospectral. Conversely, the two circulant representatives
can be equivalent under an admitted C12 automorphism if no nearest-neighbor
generator is distinguished. The state space, diagonal algebra and physical
readout interpretation must all be stated; the spectral theorem alone does
not choose them.

## 30. ST8650 — adversarial synthesis and final scope audit

This closes the thirty-round research campaign, not FIN as a physical
theory. Its decisive change from the preceding campaign is that the
investigation began with the actual archived learning update and followed
its consequences, rather than treating arbitrary compatible dynamics as
already sourced by the kernel.

### Strongest positive results

- The projected mixed-state model has exact self-consistent strict states
  and a genuine Lyapunov functional. Strict is a global minimum on its
  encoded state-spectrum orbit.
- There is an exact pure-state time-average witness, and a different
  spatially averaged update can reproduce strict from a suitably encoded
  pure state. The instantaneous raw-rule no-go is not transferred to these
  other source maps.
- In the full real isospectral learning class, the local minimum manifold
  has dimension 50, or 39 with the strict row sum also fixed. Under the
  much stronger real-circulant/positive/given-spectrum specification, an
  exact census gives one C12-automorphism orbit of kernels. This is useful
  conditional identification, not a source for the supplied spectrum.
- Positive regular projected updates can preserve the required kernel
  admissibility and an energy-decrease inequality. They are coherent
  mathematical completions after their constraints are supplied.
- A full positive Green parent of the normal-ordered strict kernel can
  be constructed exactly. A corrected fixed-covariance precision descent
  converges, and a joint Gaussian functional gives a valid bootstrap.

### Refutations and necessary qualifications

- The actual teacher-driven script learns the supplied cosine covariance;
  its verification string is unconditional and its twelve-state entropy
  target is unattainable. Neither is evidence of intrinsic FIN emergence.
- Unprojected PSD-covariance obstructions cannot be applied unchanged after
  diagonal deletion. Conversely, a valid projected mixed or averaged
  witness does not prove a pure instantaneous stationary source.
- The learned strict/legacy off-diagonal data require real covariance rank
  at least six. At strict stationary points the exact complex-rank minimum
  is seven in the interior parameter range and six at its special endpoint.
  These are not counts of particles or physical dimensions.
- Minimum rank, minimum entropy or maximum purity leaves phase choices;
  the unitary state equation conserves its spectrum and cannot dynamically
  impose those criteria on generic initial data.
- The raw learning rule need not preserve positive graph rates, regularity
  or a common instantaneous wave/heat spectral structure. Frozen duality
  remains correct. Stronger projections repair specific failures but change
  the law and do not automatically produce unique selection.
- The written precision evolution is ascent of its stated loss, not descent.
  A source-dependent covariance also requires a chain rule or a valid joint
  variational formulation. Gaussian reciprocal consistency by itself leaves
  an entire family of operators.
- The specified scalar Gibbs bootstraps cannot support seven distinct strict
  spectral sectors. Near-zero numerical residuals as decay tends to zero
  do not refute the exact positive-decay obstruction.
- Faster tracking has a Hamiltonian limiting dynamics and need not mean
  faster selection. The diagonal algebra and source projections are part
  of the mathematical model, not automatically produced by a spectrum.

### Most important conclusion

**Self-consistency, stability and derivation are different claims.**
FIN's strict kernel can be realized by carefully specified states and
learning laws; the statement that it can never self-reproduce is too broad.
But the strongest mechanisms tested here either learn information already
present in the input, retain equally good alternatives, or require new
constraints before preserving the desired dynamics. None establishes an
inevitable physical law from the kernel alone.

The most precise surviving interpretation of the audited learning sector
is a basis-dependent projected Hamiltonian–gradient covariance system,
with conditional Green/precision and propagation constructions. This
is substantive mathematics. It is not a derivation of a neural universe,
physical matter, gauge sectors, spacetime or a unique initial state.

The comparison is not a claim of mathematical priority. Covariance tracking,
spectral pinching, isospectral orbits, projected gradients and Gaussian
precision objectives are established mathematical constructions. The
contribution of this campaign is their source-faithful FIN audit, the exact
finite certificates and the explicit separation of incompatible claims.

### Highest-information next research direction

Seek an independently justified source of the **state spectrum and
coherence law**, together with an admissible spatial/algebraic constraint,
that does not encode the desired W as its input. A useful candidate must:

1. reproduce the seven-sector spectral structure without fitting those
   seven outputs back into its definition;
2. specify whether it is pure, mixed, reduced, time-averaged or spatially
   averaged, and pay for any change of state spectrum;
3. preserve its claimed positivity, geometry and propagation class, or
   explicitly predict their failure;
4. give a distinguishing dynamical or higher-order observable with a
   controlled error bound and stated operational preparation/readout.

A physically meaningful law may have initial conditions and calibrated
parameters; that alone is not a refutation. What is currently missing is
a non-tautological relationship between the proposed primordial dynamics
and its claimed physical observables. A state or field reduction with
derived memory may be an admissible route, but a supplied bath or covariance
must not be relabelled as an intrinsic source. This priority is more
informative than another high correlation against a target-defined kernel.

### Repository and verification scope

The repository-wide source/state map was used to choose this update-law
frontier. Both kernel lineages, the corrected strict gate provenance,
previous covariance, dual-dynamics, information, selector and operational
results are retained in their stated scopes. The actual archived neural
script was read and replayed; the adaptive and precision passages of the
kernel referee report were checked; the later ST2208/ST2209 scope was
preserved; existing exact transcendental enclosures were reused and tested.
The prior FAR and compendium audits remain the source map for action, RG,
gravity, units and physical-role claims, not newly generated closure proofs.

This is not a fresh proof of every historical repository statement. Each
new general claim is supported by its displayed argument; exact determinant
and assignment certificates support the finite claims, and floating
integrations are labelled as checks. Test counts and hashes are evidence
of replay, not proof of physics. No old result is counted as a new round
merely because it was re-executed.

The final scientific suite has 51 tests. All three result sets, including
the original 30,000-step run with a fixed seed, are replayable. The
requirement/evidence audit is in `COMPLETION_AUDIT.md`; execution details
are in `verification.json`. No PDF, external audit, laboratory record,
publication upload, selector closure, legacy physical-role transfer,
Standard Model/GR or ToE closure is claimed.
