# FIN: projected learning, self-consistency and persistent state motion

Krzysztof Żuchowski — research checkpoint, 8 September 2026.  
New campaign after ST8620; rounds **1–10 of 30**, ST8621–ST8630.
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

## Checkpoint and next questions

The actual source learns its teacher covariance; its advertised entropy
target and unconditional pass label do not establish FIN emergence. A
closed projected extension has a real Lyapunov structure, but its
self-consistency depends crucially on pure versus mixed state and on
instantaneous versus averaged updating. The next high-information questions
are minimum covariance-completion rank, stability of these fixed families,
and whether the positive strict cone is preserved by a physically meaningful
learning law. Rounds 11–30 remain open; no physical closure is claimed.
