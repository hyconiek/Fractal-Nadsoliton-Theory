# FIN: transition-resolved composition and the hierarchy problem

Krzysztof Żuchowski — Independent Researcher, Fractal Information Theory Project  
ORCID 0009-0002-0909-3613. Final research report, 7–8 September 2026.

**Thirty-round campaign: ST8591–ST8620.** This is the single research report;
the final synthesis and verification audit are in section 30. No PDF is
generated. Computations and tests are in this directory. “Proven” below
always includes the displayed class and premises; numerical residuals are
checks, not proofs of unqualified physical claims.

The principal result is stronger than the preceding activity audit: even the
**complete two-constituent generator** does not universally determine higher collective
dynamics. This is demonstrated for the actual strict kernel in an explicitly
infinitely extendible, reversible, maximally synchronous class. A second
result shows that a nonnegative Gram transition matrix can fail even the
necessary three-constituent event inequalities.

All theorems below are conditional on their stated mathematical class,
classical or quantum.
The standard finite-difference, moment and contraction mechanisms are not
claimed as globally new mathematics. No physical population, laboratory
record, units, selector or legacy-role transfer is inferred.

## Common scope and input

Use the frozen strict positive matrix W on C12, with cyclic distance,

\[
 W_{ij}=\frac{\cos(0.18575d+0.16250)}{1+d^{9/5}},\quad W_{ii}=0,
 \qquad Q=W-sI=-A,
\]

where s is its common row sum. These weights are positive for distances
1,…,6, and W1>W2 because both the positive cosine numerator and reciprocal
denominator decrease on these distances. The signed legacy profile is not
substituted as positive jump rates. The current AGENTS state map and the
actual replication-consistency package were consulted; this campaign
does not claim to revalidate every historical file.

The state spaces X^N below are supplied classical population spaces.
Replication into X^N is not identified with fractal refinement of X.
The strict dual calculus U_t=exp(-itA), P_t=exp(-tA) is unchanged; this
campaign studies additional collective laws and does not identify their
probability semantics with the unitary channel.

## 1. ST8591 — full transition tensor: necessary constraints and pair reconstruction

Let e=(i,a), a≠i, and f=(j,b), b≠j, label directed transitions. In a strongly
projectively consistent exchangeable Markov population family define

\[
 T_{e f}=G_2((i,j),(a,b)),\qquad q_e=Q_{ia}.
\]

Exchangeability makes T symmetric and entrywise nonnegative. For every
origin j, the block-capacity constraints are

\[
 \sum_{b\ne j}T_{(i,a),(j,b)}\le q_{(i,a)}.
\]

**Pair reconstruction [Proven].** Q and T determine G2 completely:
joint transitions have rates T, and

\[
 G_2((i,j),(a,j))=Q_{ia}-\sum_{b\ne j}T_{(i,a),(j,b)}.
\]

There is the analogous formula for changing only the second coordinate;
the diagonal is minus the total exit rate. Nonnegativity and block capacities
are sufficient for these formulas to produce a valid exchangeable pair
generator with autonomous Q. Reversibility and higher extendibility impose
additional conditions.

**Transition Gram condition [Proven].** Replicate each origin m times and
let Z_(i,a) count coordinates making transition i→a in one event. The
nonnegative Gram matrix m^-2 sum_y g(x,y) ZZ^T equals

\[
 T+\frac1m\operatorname{blockdiag}_i
       \left[\operatorname{diag}(q_{i,\cdot})-T_{ii}\right],
\]

where T_ii is the block whose two initial origins are i. For different
destinations from the same origin, one coordinate cannot make both moves,
so there are m(m−1), not m², contributing coordinate pairs. For an identical
transition the diagonal also has m q_(i,a). For distinct origins there are
m² pairs. These observations prove the formula entry by entry. Unlimited
replication therefore makes T a completely positive **real matrix** by the
closedness argument in ST8582. This is not quantum-channel complete
positivity.

**Next:** test whether these degree-two Gram and capacity conditions suffice.

## 2. ST8592 — completely positive transition data with no triple extension

Take the nonnegative length-12 vector

\[
 v=(0,1,2,3,2,1,0,1,2,3,2,1),\qquad \|v\|^2=38,
\]

let H have the cyclic shifts of v as rows, and set F=HH^T/38. This is an
explicit nonnegative Gram factorization. F is symmetric circulant, has
period six in each index, satisfies 0≤F_ij≤1, and has

\[
 F_{ii}=1,\quad F_{01}=F_{12}=16/19,\quad F_{02}=11/19.
\]

Let f(i)=i+6 mod12 and a=W6>0. At pair state (i,j), add a common antipodal
event (i,j)→(f(i),f(j)) with rate a F_ij. Subtract that rate from each
singleton antipodal transition; leave all other strict rates unchanged.

**Admissibility [Proven].** The resulting pair generator is nonnegative,
has autonomous strict Q, is exchangeable and D12 invariant. Period-six
invariance means F_(f(i),j)=F_ij, so each modified singleton rate equals
its reverse. Common events are also reversible. Non-antipodal single
transitions stay positive and connect the full population space, proving
irreducibility. The equilibrium is the product-uniform law.

Its full T is a F embedded on the twelve antipodal-transition indices and
zero on all other directed-transition rows/columns. Thus it satisfies the
complete-positive transition Gram condition and all block capacities.

**Triple obstruction [Proven].** For three event indicators Y0,Y1,Y2,

\[
 Y_0Y_2-Y_0Y_1-Y_1Y_2+Y_1\ge0
\]

holds at each of the eight Boolean configurations. Prepare origins (0,1,2)
and let each indicator mean the corresponding specified antipodal move.
Summing this inequality over any nonnegative event-rate measure requires

\[
 T_{(0,6),(2,8)}-T_{(0,6),(1,7)}
 -T_{(1,7),(2,8)}+Q_{1,7}\ge0.
\]

The proposed pair data instead give

\[
 a(11/19-16/19-16/19+1)=-2a/19<0.
\]

Hence **no triple generator with these pair projections exists**, regardless
of how its extra events are designed. This is not a failed numerical search.
The exact Boolean facet and rational Gram matrix are checked directly.
Strict marginal/reversal/dihedral replay residuals are below 10^-13.

This is a classical compatibility inequality for simultaneous events, not
evidence for a quantum Bell violation. Its lesson is that Gram positivity
is necessary but weaker than a realizable joint-event law.

**Next:** use genuinely infinitely extendible examples, so that failure of
extension cannot explain any remaining nonuniqueness.

## 3. ST8593 — identical complete pair laws, distinct triple laws

Set R0=W/s. Let C1,C2 be the unweighted cyclic-distance-one and -two
adjacency matrices, and choose

\[
 h=(R_0)_{02}/2>0,\quad D=h(C_1-C_2),\quad R_u=R_0+uD,
 \qquad -1\le u\le1.
\]

Each R_u is symmetric, circulant, stochastic, zero on the diagonal and
strictly positive off the diagonal. Indeed D has zero row sums and only
changes the distance-one/two entries; R0,1>R0,2 and h=R0,2/2 keep both
endpoints of the allowed interval positive. The curve D is a constructed
comparison direction, not a strict-sourced fluctuation law.

For a probability law μ on [-1,1] of mean zero define

\[
 G_N^\mu=s\left(\int R_u^{\otimes N}\,\mu(du)-I\right).
\]

Interpretation: a common Poisson clock of rate s, followed at each tick by
conditionally independent coordinate updates with the sampled R_u. A new
independent u is sampled at every tick; a different temporal sampling law
would be another model.

**Structural properties [Proven].** These generators are Markov,
exchangeable, strongly projectively consistent, reversible with respect
to product-uniform equilibrium, and D12 invariant. Tensor restriction uses
R_u 1=1. Symmetry of each R_u gives reversibility. Positivity of all
off-diagonal entries makes a two-tick transition between any population
configurations possible, proving irreducibility. Mean zero gives G1=Q.
Every coordinate changes at every tick, so the aggregate activity matrix
is the same maximal B=s 11^T for every μ.

Take

\[
 \mu_A=\tfrac12\delta_{-1/2}+\tfrac12\delta_{1/2},\qquad
 \mu_B=\tfrac45\delta_{-1/4}+\tfrac15\delta_1.
\]

Both have m0=1, m1=0, m2=1/4. Their third moments are 0 and 3/16.
The tensor R_u^⊗N is polynomial of degree N in u. Consequently the
**entire matrices G1 and G2 agree**, not just their spectra or aggregate
activities. They therefore give identical one- and two-coordinate path
laws for common initial conditions in this model class.

For the triple transition (0,0,0)→(1,1,1), the rate difference is exactly

\[
 \Delta q=s\,(3/16)h^3>0.
\]

For strict this is about 0.00006021899861. The two rates are approximately
0.03883869026 and 0.03889890926. The rational third-moment difference and
positivity of s,h prove nonidentity; floats only check the application.

**Next:** determine whether some larger finite population cutoff eliminates
the ambiguity.

## 4. ST8594 — no finite hierarchy cutoff determines the law

For any m≥1 let n=m+1, u_k=-1+2k/n, and define probability laws

\[
 \mu_{\rm even}=2^{-m}\sum_{k\text{ even}}\binom{m+1}{k}\delta_{u_k},
 \quad
 \mu_{\rm odd}=2^{-m}\sum_{k\text{ odd}}\binom{m+1}{k}\delta_{u_k}.
\]

**Finite-order no-go [Proven].** These laws both have mean zero and agree
in every moment through degree m. They disagree at degree m+1. Therefore
their full population generators agree for every N≤m and differ for N=m+1,
while retaining every structural property proved in round 3.

**Proof.** The alternating binomial sum annihilates every polynomial of
degree at most m. On the monomial u^(m+1), its value after normalization is

\[
 (-1)^{m+1}\frac{(m+1)!}{2^m}
                  \left(\frac2{m+1}\right)^{m+1}\ne0.
\]

The even and odd binomial masses each total 2^m. Their common first moment
is the centered binomial first moment, zero. Tensor polynomiality gives
equality of all lower generators. The all-zero to all-one transition at
order m+1 has a nonzero difference s h^(m+1) times the displayed moment
gap. This proves the result for arbitrary m; the script checks m=1,…,12
using exact fractions as a replay, not as a substitute for the proof.

This disproves a finite-cutoff identification theorem in the supplied
common-clock class. It does **not** prove that a finitely stated physical
law cannot generate the hierarchy, or that finite-order predictions are
useless. Nor does it assert nonuniqueness for every individual moment
vector: zero variance, for example, fixes a point mass. Those distinctions
must be retained when examining extremal or finitely parametrized sources.
They motivate the next rounds.

## 5. ST8595 — all ideal orders determine the scalar noise law

Let p=(R0)_01 and h=D01≠0. The all-zero to all-one N-coordinate rate q_N
gives z_N=q_N/s=E[(p+hU)^N]. With m0=1,

\[
 m_N=h^{-N}\left[z_N-\sum_{k<N}\binom Nk p^{N-k}h^k m_k\right].
\]

**Identification [Proven, within this supplied curve and clock].** All
ideal population orders determine every moment and hence the probability
law μ on [-1,1]. To prove the last step, equal moments give equal integrals
of polynomials. Polynomial uniform approximation on a compact interval
then gives equal integrals of every continuous function, which determine
the probability measure. The exact triangular inversion is replayed on a
rational example through order eight.

This identifies μ only after the curve R_u and clock have been supplied.
It does not derive them, grant access to every population size, or select
an absolute time unit.

## 6. ST8596 — nonidentification is not identical to operational instability

The even/odd laws in round 4 have disjoint supports, so their latent-law
total variation is exactly one. Nevertheless

\[
 W_1(\mu_{\rm even},\mu_{\rm odd})=\frac2{m+1}.
\]

**Proof.** On the ordered nodes their cumulative signed mass is
(-1)^k binom(m,k)/2^m. The one-dimensional transport cost is the integral
of the absolute cumulative difference. Summing over the intervals of
length 2/(m+1) yields the displayed value. All these identities are checked
with exact fractions.

For arbitrary μ,ν on the supplied noise curve, a single row of R_u and
R_v differs in TV by 2h|u−v|. Product coupling gives an N-coordinate row
bound 2Nh|u−v|. Integrate an optimal coupling of μ,ν and use Duhamel's
identity for the two finite generators to obtain

\[
 \sup_x\operatorname{TV}(e^{tG_N^\mu}(x,\cdot),e^{tG_N^\nu}(x,\cdot))
 \le\min\{1,2sNh\,t\,W_1(\mu,\nu)\}.
\]

Thus exact latent-law identification can fail while selected finite
observations remain close. The bound concerns transition distributions,
not complete continuous-time path records. For N≤m the matched examples
have exact equality, a stronger result than the transport bound.

## 7. ST8597 — inversion conditioning has to be typed

Holding the lower moments exact, an error ε in z_N changes the reconstructed
m_N by exactly ε h^-N. For strict, h≈0.05783373781, so h^-8 is about
7.9901×10^9. More generally,

\[
 m_N=h^{-N}\sum_{k=0}^N\binom Nk(-p)^{N-k}z_k
\]

gives an algebraic error bound ε[(1+p)/h]^N when every raw z_k is known
within ε. This is a conditioning statement about this inversion, not a
minimax impossibility theorem for all estimators. Some arbitrary raw-moment
perturbations leave the positive-measure moment cone and are not physical
alternatives. The feasible even/odd alternatives from round 4 supply the
separate exact nonidentification witness.

**Next:** test whether a finite recursive law can fix the hierarchy without
claiming that large moment tables must themselves be fundamental.

## 8. ST8598 — a finite self-similar law conditionally closes all orders

Supply a contraction 0≤r<1 and independent fair signs ε_k∈{-1,1}. The law
of

\[
 U=(1-r)\sum_{k=0}^\infty r^k\varepsilon_k
\]

is the unique probability law on [-1,1] satisfying

\[
 \mu=\tfrac12(F_+)_*\mu+\tfrac12(F_-)_*\mu,
 \quad F_\pm(u)=ru\pm(1-r).
\]

**Existence and uniqueness [Proven, conditional].** The series converges
absolutely and stays in [-1,1]. Coupling two instances with the same first
sign contracts their Wasserstein distance by at most r. Iteration therefore
has one fixed probability law; equivalently the tail influence of any
initial value is at most 2r^n after n steps.

Its mean is zero and variance is

\[
 v=\frac{1-r}{1+r},\qquad r=\frac{1-v}{1+v}.
\]

Thus a supplied pair variance 0<v≤1 determines r **within this specific
two-map, fair-sign, fixed-amplitude class**. All higher moments follow from
the distributional equation U=rU'+(1-r)ε with independence. For v=1/4,
r=3/5 and m4=35/272. Exact moments through degree ten are computed.

This is genuine finite recursive specification of an infinite collective
hierarchy after the curve and fresh-clock sampling have been supplied.
It is not a strict-derived physical law. At r=3/5 the map images overlap;
no Cantor-dimension formula, absolute continuity or spacetime dimension is
asserted. A stationary self-similar distribution also does not specify the
temporal correlations of successive clock samples.

## 9. ST8599 — the innovation rule is indispensable

Replace the fair binary innovation by ξ taking values -1,0,1 with
probabilities 1/4,1/2,1/4, and use r=1/3. The same contractive construction
has variance 1/4, just like round 8, but fourth moment 11/80 rather than
35/272.

For general centered independent innovation ξ, the exact recurrence is

\[
 m_n=\frac{\sum_{k<n}\binom nk r^k(1-r)^{n-k}
                   m_k\,E[\xi^{n-k}]}{1-r^n}.
\]

It verifies the different fourth moments without numerical fitting.
Both resulting laws are symmetric and self-similar and generate valid
strict population models with identical first and second moments. A claim
that pair information plus unspecified “fractal compression” determines
the full law is therefore **refuted**. The stated binary innovation axiom,
not fractality alone, paid the uniqueness in round 8.

## 10. ST8600 — spatial symmetry does not supply the noise symmetry

Every matrix R_u on the chosen curve is individually D12 invariant.
Reflection preserves the distance-one/two contrast D; it does not send
D to -D. Therefore the action of graph symmetries on the parameter u is
trivial. Every μ, including the asymmetric μ_B in round 3, has the same
spatial graph symmetries. Neither a fair sign nor a zero third moment
follows from those symmetries.

The old mirror-odd carrier concerned a different representation type.
Importing its branch exchange into this even radial noise coordinate would
be a new axiom. The historical radial damping/path-counting formulas also
do not specify binary versus ternary innovations, a contraction r, or
fresh independent sampling at event times. The finite recursive model is
therefore a candidate **source-law form**, not an already sourced result.

## Outcome after round 10

The new strict-compatible no-go is not merely that eigenvalues omit
geometry: even complete finite-order collective Markov laws omit higher
composition. Conversely, all ideal orders identify a supplied scalar law,
and a specified self-similar recursion can compress that entire law into
a finite prescription. The unresolved issue is which prescription, if any,
is actually derived from FIN, including its temporal sampling semantics.

This led to the temporal and clock tests below. The remaining goal is not
redefined as merely completing this composition subproblem.

## 11. ST8601 — the same recursive law supports inequivalent temporal dynamics

Let μ be the binary affine fixed law from round 8, with contraction r=3/5
and variance v=1/4. Both the stationary recursion

\[
 U_{n+1}=rU_n+(1-r)\varepsilon_{n+1}
\]

and independent resampling U_n~μ have exactly this one-time distribution.
Their lag-k covariances are respectively v r^k and zero for k>0.
The covariance formula follows by iterating the affine recursion and using
the independence and zero means of the innovations.

**Not just a clock difference [Proven].** Apply the recursion at rate-one
Poisson ticks. The centered functions u and u²−v have decay rates
1−r and 1−r². This follows from the exact conditional expectations
E[U'|u]=ru and E[U'^2−v|u]=r²(u²−v). The second observable has positive
variance (m4−v²=9/136), so its mode is not vacuous. A rate-ρ full-reset
process has rate ρ for *every* centered observable. Matching the first
mode sets ρ=2/5, but the affine process's second mode has rate 16/25,
not 2/5. No one scalar time calibration matches both.

Thus an invariant self-similar probability law, even if uniquely specified,
does not by itself define its time evolution or the sampling of operator
updates. These are hidden-parameter models, not claims of physical clocks.

## 12. ST8602 — averaging a generator does not average its dynamics

Take Q=-A and C=sD from the supplied radial noise curve. Both Q+C and Q−C
are symmetric Markov generators. Let an independent binary environment
switch between signs at rate κ>0 and use the joint generator

\[
 L=\begin{pmatrix}Q+C-\kappa I&\kappa I\\
                  \kappa I&Q-C-\kappa I\end{pmatrix}.
\]

Prepare the environment with equal probabilities, independent of the
initial vertex. Collapse it after evolving. With J=(I,I) and E=(I,I)^T/2,
the observed transition matrix is S_t=J exp(tL) E.

**Counterexample [Proven].** Direct block multiplication gives

\[
 S'_0=Q,\qquad S''_0=Q^2+C^2.
\]

Since C is a nonzero real symmetric matrix, C²≠0. Therefore S_t is not
exp(tQ), and it is not a time-homogeneous semigroup on the observed vertex
space. If it were a differentiable semigroup, its generator S'_0 would
force S''_0=Q². The joint process and every observed S_t remain positive
and probability preserving. What fails is the coarse homogeneous Markov
description, not probability conservation.

The strict test at κ=0.7 gives ||C²||_F≈0.191638 and a
||S_0.5−S_0.2 S_0.3||_F defect ≈0.00271211. These numbers supplement
the block-identity proof. Matching the instantaneous drift alone is not
matching the complete strict singleton path law used in earlier rounds.

## 13. ST8603 — exact memory from the hidden temporal law

Write p=p_++p_- and m=p_+−p_-. The joint equations give

\[
 \dot p=Qp+Cm,\qquad \dot m=Cp+(Q-2\kappa I)m.
\]

For the specified preparation m(0)=0, elimination is exact:

\[
 \dot p(t)=Qp(t)+\int_0^t C e^{(Q-2\kappa I)(t-u)}C p(u)\,du.
\]

The Laplace response, for real z>0, is

\[
 \widehat S(z)=\left[zI-Q-C(zI-Q+2\kappa I)^{-1}C\right]^{-1}.
\]

The direct 24-state resolvent and this Schur expression agree numerically
to below 10^-13. The analytic derivation does not assume C commutes with
Q. A nonzero hidden initial imbalance would add an inhomogeneous term
C exp[(Q−2κI)t] m(0), another preparation datum.

This constructs a missing temporal-response object from a *supplied*
environment law. Neither the stationary sign distribution nor Q fixes κ.
No physical environment or information substrate underneath the nadsoliton
has been inferred. Positive damping in the legacy distance formula is not
being reinterpreted as this time-memory kernel without a separate map.

## 14. ST8604 — fast switching recovers heat only as a controlled limit

For the chosen circulant Q,C, diagonalize them jointly with modal values
q≤0,c∈R. Put ω=sqrt(κ²+c²). The observed mode multiplier is exactly

\[
 S_{q,c}(t)=\tfrac12(1+\kappa/\omega)e^{(q-\kappa+\omega)t}
          +\tfrac12(1-\kappa/\omega)e^{(q-\kappa-\omega)t}.
\]

This solves the two-dimensional modal system and its initial conditions.
The factor multiplying exp(qt) is at least one: its derivative is
exp(−κt)c² sinh(ωt)/ω≥0 and its initial value is one. It is at most
exp[(ω−κ)t], because the two exponential weights are nonnegative and sum
to one. Since ω−κ=c²/(ω+κ)≤c²/(2κ),

\[
 \|S_t-e^{tQ}\|_2\le e^{\|C\|_2^2t/(2\kappa)}-1.
\]

Thus heat is recovered uniformly on bounded time intervals as κ→∞.
It is not obtained exactly at finite κ for nonzero C. At t=0.7, the
observed norm errors for κ=1,5,20,100 are approximately 0.00465315,
0.00172377, 0.000483691 and 0.0000995145, each below its analytic bound.
The slow modal exponent is q+c²/(2κ)+O(κ^-3); treating that truncation
as a separate Markov generator would require its own positivity check.

## 15. ST8605 — finite data can identify special extremal laws

The finite-cutoff no-go in round 4 is not a theorem of universal
nonidentification at every data point. For U∈[-1,1], E U=0:

- E U²=0 forces U=0 almost surely.
- E U²=1 forces U∈{-1,1} almost surely; the mean fixes equal endpoint
  weights.

These follow from the zero expectation of the nonnegative functions U²
and 1−U². A supplied maximum-variance principle therefore selects a unique
law on this supplied compact one-dimensional curve. It does not select
the curve, its endpoints or the principle itself.

**Finite-support certificate [Proven].** If the known moment matrix has
a nonzero null polynomial p, then E[p(U)²]=0, so any representing measure
is supported on the finitely many zeros of p in [-1,1]. Its weights are
uniquely fixed by the corresponding finite Vandermonde system whenever
those moments are supplied. For the symmetric ±1/2 law,

\[
 E[(U^2-1/4)^2]=m_4-\tfrac12m_2+1/16=0.
\]

The asymmetric law from round 3 gives 9/64 instead; the exact certificate
correctly distinguishes the laws once fourth-order data are included.

For the binary affine law with 0<r<1, the support is infinite: it contains
1 and the distinct points 1−2(1−r)r^k. A nonzero polynomial cannot vanish
on that support. Hence its moment matrices do not acquire such a finite
null-polynomial certificate. Its uniqueness in round 8 came from the
recursive source-class assumption, not flat finite moment data.

## 16. ST8606 — rigidity of an exactly spectrum-preserving positive random clock

Consider the explicitly restricted class of common subordinators with no
killing, drift d≥0 and nonnegative Lévy measure ν on (0,∞), satisfying
the usual finite-exponent condition integral min(1,τ) ν(dτ)<∞. The scalar
decay exponent is

\[
 f(\lambda)=d\lambda+\int_{(0,\infty)}(1-e^{-\lambda\tau})\,\nu(d\tau).
\]

This is the stated class assumption, not a conclusion of the spectral
theorem. Subordinating a supplied heat semigroup gives exp[−t f(A)].

**Two-mode rigidity [Proven].** If f(λ1)=cλ1 and f(λ2)=cλ2 for
0<λ1<λ2 and c≥0, then ν=0 and d=c.

**Proof.** For every τ>0,

\[
 \frac{1-e^{-\lambda\tau}}{\lambda}
       =\int_0^\tau e^{-\lambda u}\,du
\]

is strictly decreasing in λ. A nonzero positive ν therefore makes
f(λ1)/λ1 strictly larger than f(λ2)/λ2. The asserted equality rules
this out, forcing ν=0, after which f(λ)=dλ. The integrability premise
makes both positive integrals finite. No cancellation with negative
clock-jump weights is permitted in this class.

The strict A has at least two distinct positive eigenvalues. Otherwise a
connected symmetric Laplacian with only eigenvalues 0 and λ would equal
λ(I−11^T/12), with all off-diagonal weights equal; strict W1≠W2.
The numerical extremes are about 0.7541211542 and 2.3421820411.
Thus f(A)=A forces the deterministic identity clock in this class.
Allowing an overall calibration gives f(A)=cA and leaves c free.

This does not derive the independent tensor base used in some population
constructions. The same rigidity can be applied to a correlated base with
the same singleton A, and would simply leave that correlated base unchanged.

## 17. ST8607 — robust bounds detect finite clock jumps, not arbitrary tiny ones

For λ1<λ2 define

\[
 g(\tau)=\frac{1-e^{-\lambda_1\tau}}{\lambda_1}
        -\frac{1-e^{-\lambda_2\tau}}{\lambda_2}.
\]

Then g(τ)>0 and g'(τ)=exp(−λ1τ)−exp(−λ2τ)>0. The measured spectral
shape defect is Δ=f(λ1)/λ1−f(λ2)/λ2=integral g(τ) ν(dτ). Therefore

\[
 \nu([\tau_0,\infty))\le\frac{\Delta}{g(\tau_0)},\qquad \tau_0>0.
\]

An atom at τ0 attains this bound, so this statement is sharp for the
declared information. If individual exponent errors are bounded, replace
Δ by its justified upper bound; a sampled near-zero is not an exact
deterministic-clock certificate.

Tiny jumps remain different. The measures ν_h=h^-1 δ_h have diverging
total jump rate but f_h(λ)=(1−exp(−hλ))/h→λ. Their spectral shape
defect tends to zero. Finite precision therefore does not uniformly bound
the total mass of ν near zero. This is the familiar drift limit, not an
experimental detection of a new time scale.

## 18. ST8608 — one nonzero mode is insufficient

For a two-state generator with spectrum {0,1}, choose ν=(1/2)δ_1 and
d=1−(1−exp(−1))/2>0. Then f(1)=1 but

\[
 f(2)=\tfrac32+e^{-1}-\tfrac12e^{-2}\approx1.8002117996<2.
\]

The singleton generator is matched exactly, while its independent
two-coordinate base has a mode at 2 whose subordinated rate changes.
Both generators are genuine positive Markov generators. This is an
explicit premise-removal test for the two-mode requirement; one matched
mode or one fitted spectral gap does not eliminate a random clock.

## 19. ST8609 — commuting stochastic updates need not be heat mixtures

R0=W/s is a stochastic symmetric matrix commuting with A. It has trace
zero and eigenvalue one. Hence at least one other eigenvalue is negative
(numerically the minimum is about −0.4106919069). It cannot be a positive
mixture of exp(−τA), whose eigenvalues are nonnegative. The same trace
argument applies to the zero-diagonal R_u curve in the first checkpoint.

Thus the earlier common-event counterexamples do not contradict
subordination rigidity: their allowed jump maps belong to a larger class.
Entrywise probability positivity, positive-semidefinite operator positivity,
commutation with A, and positive heat-subordination structure are distinct
premises. Requiring the last one is additional mathematical structure.

## 20. ST8610 — source and kernel audit of the new clock theorem

Before this FAR source reconciliation, K1, K2, F2, F3, S2 and
SUMMARY_GROK.md were read completely. Their priorities were compared with
the later current AGENTS state-map rules, not treated as an automatic next
task. Relevant later guards include P298/P300, P432, ST2622–ST2801 and
the current composition frontier.

| Evidence inspected | What remains established | What the new result does not supply |
|---|---|---|
| K1/K2 | Distinct canonical legacy and gate-selected strict provenance | An ontological identification or a source for D, ν or temporal sampling |
| F2 | Explicit FAR input classification and no silent role inheritance | Permission to turn a strict benchmark clock theorem into a legacy physical-role theorem |
| F3 plus later state-map guards | Route-scoped defect/source obligations; old immediate priorities can be stale | A reason to reopen closed m2, selector or generic bridge loops |
| S2 plus current state-map-first rule | Legacy remains an intermediate bridge object, with separate completion and role transfer | A completed bridge or an automatic clock source |
| P298/P300 guards | Specific Bernstein/distance-mixture bridge obstructions | A contradiction with the new self-spectrum clock theorem; its arguments are spectral λ and clock τ, not distance d |
| ST2622–ST2801 | Equilibrium, kinetic activity, clock, circulation and apparatus remain typed separately | A supplied global subordinator or an SI time unit |
| SUMMARY_GROK and unit/selector guards | Dimensionless phase/information identities and open source cuts | Seconds, action units, selected orientation or a unique population law |

The decisive additional premise is now explicit: **all allowed temporal
modifications must be a positive homogeneous common subordinator of a
supplied base semigroup**. None of the inspected source packets derives
that class restriction. If the base is chosen to be sum_k A^(k), its
independent composition has already been supplied. Applying the same
clock argument to the correlated bases constructed earlier leaves their
correlations intact. The theorem removes clock-jump freedom within its
class; it does not pay for the composition class itself.

No old bridge search was rerun, and no absence claim is extrapolated to
every unread historical file. This source audit prevents promoting the
new conditional rigidity result into a false physical closure. The most
useful remaining target is a source-defined event/temporal law that both
passes the hierarchy constraints and fixes its own admissible class.

## Outcome after round 20

Rounds 11–20 add real temporal obstructions, an exact memory model, a
controlled fast-noise limit, finite-identification exceptions and a sharp
positive-clock rigidity theorem with removal and robustness tests. The
suite at that checkpoint had 31 scientific tests. The remaining investigations
below address the full law, quantum alternatives and operational meaning.

## 21. ST8611 — the full common-event law has an invisible identity convention

Let S be the compact space of n×n row-stochastic matrices and let ν be a
finite nonnegative intensity measure on S. Define

\[
 G_N^\nu=\int_S(R^{\otimes N}-I)\,\nu(dR).
\]

These are legitimate finite Markov generators and are strongly consistent
under coordinate deletion. This is a specified model class, not a claim
that every conceivable FIN dynamics must have this form.

**Identity gauge [Proven].** Replacing ν by ν+cδ_I, c≥0, leaves every G_N
unchanged. Such a clock event attempts the identity update and has no
state-path record. Consequently even all population dynamics cannot identify
the *attempted* total event rate without a convention excluding idle events
or an additional record of those attempts. This is an observational gauge
for the declared records, not a physically distinguishable coupling.

For strict, intensity sδ_(W/s) and the same measure plus 5δ_I give exactly
identical generators, despite different total attempted clock rates. The
test checks this at several population sizes.

## 22. ST8612 — all orders identify a finite intensity law modulo that gauge

**Uniqueness [Proven in the finite-intensity class].** If two finite measures
ν,ν' give the same G_N for every N, then

\[
 \nu-\nu'=c\delta_I
\]

for a real constant c. In particular, imposing ν({I})=ν'({I})=0 gives
ν=ν'.

**Proof.** Put σ=ν−ν'. Any monomial M in matrix entries can be written as
M(R)=product_k R_(i_k,j_k) and hence appears in an entry of R^⊗N. Equality
of the corresponding generator entries gives
integral M dσ=σ(S)M(I). This also holds for constants. Polynomials in
matrix coordinates separate points of the compact stochastic-matrix space
and uniformly approximate continuous functions on it. Thus σ and
σ(S)δ_I agree against every continuous function and are equal as finite
signed measures. Excluding the identity atom fixes c=0.

For finite ν the observable event rate away from identity can also be
recovered as a limit of exit rates from configurations with m copies at
each origin:

\[
 \int\left[1-\left(\prod_iR_{ii}\right)^m\right]d\nu
 \longrightarrow\nu(S\setminus\{I\}).
\]

Monotone convergence proves this because a stochastic matrix has every
diagonal entry equal to one exactly when it is I. This is an ideal
all-orders identification result, not a finite experimental protocol or
a derivation of ν from W.

## 23. ST8613 — independent dynamics is a singular limit, not a nontrivial finite common clock

If a common-event measure gives independent nonzero singleton dynamics,
its same-origin pair-departure rates would satisfy

\[
 b_{ii}=\int(1-R_{ii})^2\,d\nu=0.
\]

Positivity forces every row to stay put almost everywhere, hence R=I
almost everywhere and the singleton generator is zero. Therefore a
nonzero independent tensor generator cannot be represented solely by
such an ordinary nonnegative common-event measure with zero pair activity.

It is nevertheless a controlled limit. For a finite Markov generator Q,
choose 0<ε≤1/max_i(-Q_ii) and R_ε=I+εQ. Then

\[
 G_N^{(\varepsilon)}=\varepsilon^{-1}
       [(I+\varepsilon Q)^{\otimes N}-I]
 =\sum_kQ^{(k)}+\sum_{j=2}^N\varepsilon^{j-1}
        \sum_{|J|=j}\prod_{k\in J}Q^{(k)}.
\]

For fixed N the remainder norm is bounded by
sum_(j=2)^N binom(N,j) ε^(j−1)||Q||^j and tends to zero. The common clock
rate diverges, while a two-coordinate specified joint rate is order ε.
The numerical and exact second-order expansion tests agree with this bound.

Thus exact representation type need not be stable under a small observed
generator perturbation. A separate independent-transition component is
needed for a nonsingular general specification that includes the limit.

## 24. ST8614 — a minimal identifiable law in the declared mixed class

Consider the explicitly declared family

\[
 G_N=\sum_{k=1}^NQ_0^{(k)}+
       \int_{S\setminus\{I\}}(R^{\otimes N}-I)\,\nu(dR),
 \qquad \int w(R)\,\nu(dR)<\infty,
\]

where Q0 is a finite singleton Markov generator and
w(R)=sum_(i≠a) R_ia. The positive measure ν may have infinite mass near I,
but has no identity atom. Since the integrand row norm is at most 2Nw(R),
the integral defines a finite generator for each finite N.

**All-orders uniqueness [Proven within this representation class].** The
complete hierarchy determines Q0 and ν uniquely.

**Proof.** Every transition in which at least two specified coordinates
change has rate integral product R_(i_k,a_k) dν, with a_k≠i_k. The
independent term does not contribute. Thus all monomials of degree at
least two in off-diagonal matrix entries are known. Let τ=w²ν, a finite
measure because w≤n and integral w dν is finite. Every polynomial moment
of τ is known: expand w² as a sum of products of two off-diagonal entries.
Off-diagonal entries determine the entire stochastic matrix, and compact
moment determinacy therefore identifies τ. Away from I, w>0, so ν=w^-2 τ
is uniquely determined. Finally

\[
 (Q_0)_{ia}=Q_{ia}-\int R_{ia}\,\nu(dR),\qquad a\ne i,
\]

recovers Q0. The assumed integrability makes these integrals finite.

This is **uniqueness given existence in the declared class**, not an
unproved universal representation theorem for all exchangeable or
non-Markov FIN models. Independent drift and a common-event law are both
needed to include the removal examples; idle-event intensity has already
been quotiented out.

Related exchangeable-process constructions are part of established
probability theory; see the author paper by
[Crane and Lalley](https://arxiv.org/abs/1307.1713). The theorem used here
is the displayed class-scoped uniqueness argument, not an imported claim
that the whole FIN repository belongs to a universal representation class.
Identification uses fixed state labels and full transition data, not only
unlabelled spectra.

An exact finite-support check illustrates identification, not just forward
substitution. On symmetric two-state update matrices with flip probability
x, take common intensity atoms of masses 2 and 1 at x=1/4 and 3/4. Since
w=2x, τ has masses 1/2 and 9/4. Its first five moments are

\[
 (11/4,\ 29/16,\ 83/64,\ 245/256,\ 731/1024).
\]

Exact moment elimination recovers the polynomial
x²−x+3/16, its two roots and the two weighted masses. Its squared integral
is exactly zero, certifying the support rather than merely fitting two
atoms. Dividing by (2x)² recovers ν. A supplied total singleton flip rate
13/8 then recovers the independent rate 3/8. This is a positive finite
identification example of the extremal kind in round 15, not a general
finite-cutoff identification theorem.

## 25. ST8615 — the finite-hierarchy ambiguity also occurs for quantum channels

The classical examples do not alone decide a quantum interpretation. Here
is a separate standard-quantum construction. Let P0 be the spectral
projector onto the strict constant mode, and let

\[
 V_\theta=e^{i\theta P_0},\qquad
 \mu_\pm(d\theta)=\frac{1\pm\epsilon\cos((m+1)\theta)}{2\pi}\,d\theta,
 \quad 0<\epsilon\le1.
\]

These are positive probability densities. Define collective random-unitary
channels by averaging conjugation with V_θ^⊗N against μ±. Such channels
are completely positive, trace preserving, unital, permutation covariant
and consistent under partial traces. P0 commutes with A and is graph
invariant; no directed selector is inserted.

The channel multiplies a coherence between tensor sectors with difference
k in their P0 occupation by the Fourier coefficient of μ± at k, where
|k|≤N. Those coefficients agree for all |k|≤m. The entire channels
therefore agree for N≤m but differ at N=m+1, where their multiplier
difference is ε. A GHZ-like superposition of the all-P0 and all-complement
sectors gives output trace distance ε/2. The two-level subspace test at
m=2, ε=3/4 gives zero channel differences through N=2 and trace distance
3/8 at N=3. No 12^N-dimensional density matrix is needed to witness this.

Adding these channels at Poisson event times to the supplied coherent
Hamiltonian sum_k A^(k) also gives valid open quantum semigroups with the
same lower-order generators. These examples do **not** keep the bare
closed U_t as the complete observed noisy singleton channel; they share
the same *open* singleton completion. Quantum postulates and instruments
are assumed, not derived from FIN.

**Extremal falsification check.** Exact unitary singleton marginals for
every joint input are a much stronger premise and can force tensor-product
evolution. To see this, conjugate away the specified local unitaries.
The normalized Choi state of the joint channel then has a pure maximally
entangled marginal on each local input-output pair. A positive state with
a pure marginal factorizes across that marginal: all orthogonal diagonal
blocks vanish and positivity forces their off-diagonal blocks to vanish.
Iterating gives the product Choi state and therefore the product channel.
An entangling controlled-NOT is not a counterexample: it maps a pure
singleton input into a mixed marginal for a suitable product input.
The test verifies this premise failure explicitly.

Thus the quantum no-go is generic over the displayed open completion class,
not a denial of finite extremal-channel identification. Permutation covariance
here is not a derivation of bosonic or fermionic particle statistics.

## 26. ST8616 — even a fixed complete law needs operational preparation and readout

Let u be uniform, J=11^T/n, and choose nondegenerate detector channels
M_c=cI+(1−c)J with 0<c<1. Prepare p_b=u+b(e_i−u). For strict heat,
M_c commutes with P_t and

\[
 M_cP_t p_b=u+bc\,P_t(e_i-u).
\]

Hence (b,c)=(0.6,0.8) and (0.75,0.64) give identical complete observed
heat trajectories at all times although preparation purity and detector
contrast differ. Both detector matrices are invertible; this is not a
device that reports a constant outcome. An independently established
preparation or detector calibration can break the ambiguity.

The result does not say that physics must derive every initial condition.
It says that a prediction needs specified operational inputs. A finite
matrix law alone does not identify which experimental procedure prepares
e_i or measures its population. This is an implementation/calibration
obligation separate from selecting the dynamics.

## 27. ST8617 — the hierarchy obstruction is not an artifact of choosing strict

Return to the canonical signed legacy V and use the explicit positive
two-sheet cover

\[
 \widetilde W=\begin{pmatrix}V_+&V_-\\V_-&V_+\end{pmatrix},
 \qquad R_L=\widetilde W/s_L.
\]

This is the previously proved signed-cover construction, not a new
legacy-to-strict completion. It retains a signed odd sector and an unsigned
even sector; no legacy physical role is transferred.

Within each sheet the distance-one and antipodal legacy edges are positive.
Choose an equal-sheet perturbation D_L with +h_L on both distance-one
edges and −2h_L on the antipodal edge, where
h_L=min((R_L)_01/2,(R_L)_06/4)>0. Its row sums vanish, it respects the
graph and sheet symmetries, and R_L+uD_L stays stochastic for |u|≤1.

The two laws in round 3 therefore give identical full one- and two-copy
generators on this 24-state cover and different triple generators, with
rate gap s_L(3/16)h_L³>0. Numerically the pair error is below 10^-13 and
the triple gap is about 0.00002369179618. The cover is connected, and the
positive triangle (plus 0, plus 1, minus 3) makes it aperiodic; two-sheet
population realizations are not excluded by a disconnected-support artifact.

Thus simply replacing strict by a positive realization retaining legacy
signs does not eliminate the higher-law freedom. This statement is about
the supplied cover. It is not an identification of either raw kernel with
the other, nor a claim that every possible signed formulation is classical.

## 28. ST8618 — approximate recursive sources give controlled predictions

Let H be the binary affine probability-map operator from round 8, with
W1 contraction r<1 and fixed law μ*. For any proposed probability law μ,

\[
 W_1(\mu,H\mu)\le\varepsilon
 \quad\Longrightarrow\quad
 W_1(\mu,\mu_*)\le\frac{\varepsilon}{1-r}.
\]

**Proof.** The triangle inequality and contraction give
W1(μ,μ*)≤ε+r W1(μ,μ*). Rearrangement proves the bound.
Combined with round 6, this gives the operational transition bound

\[
 \operatorname{TV}\le
 \min\{1,2sNh\,t\,\varepsilon/(1-r)\}.
\]

There is also a finite-information approximation of the exact recursive
law. Truncate its bit series after K digits. Coupling the truncation to
the full series bounds its error by W1(μ_K,μ*)≤r^K. For r=3/5 its
variance deficit is exactly (1/4)r^(2K). The finite rational laws, their
recursion residuals and the bounds were computed for depths 1–6, with
an independent depth-eight comparison.

This is a genuine conditional benefit of a finite self-similar specification:
it controls an infinite hierarchy without listing every high-order law.
The contraction, innovation distribution, embedding R_u and fresh sampling
remain source premises. The theorem does not manufacture them from the
observed kernel or certify a physical fractal dimension.

## 29. ST8619 — relational clocks remove a unit convention, not a dynamical law

An inability to derive an absolute SI second is not, by itself, a proof
that a model cannot yield physical predictions. Consider a supplied
Markov system G and an independent reference clock with Poisson rate ω.
The system transition matrix when the first reference tick occurs is

\[
 K_\omega=\int_0^\infty\omega e^{-\omega t}e^{tG}\,dt
          =\omega(\omega I-G)^{-1}.
\]

This is stochastic, and the record at the nth tick has transition matrix
K_ω^n. Joint rescaling (G,ω)→(aG,aω) leaves every such prediction
unchanged. For strict this is the normalized Green operator
ω(ωI+A)^−1. It gives an exact relational connection between a Green
family and an internal reference-clock protocol.

This is **random observation stopping**, not subordination of the system's
physical evolution. It therefore does not contradict round 16. A periodic
reference clock would give a different sampling law. Which clock exists,
how it couples, and how it is read are operational premises.

Ratios such as G/ω remain consequential: changing ω alone changes the
predictions. Likewise, the two third-order laws in round 3 differ even in
their dimensionless triple/singleton rate ratio, so their ambiguity is not
an absolute-unit gauge. Initial data may be legitimate inputs of a physical
theory, whereas the absence of a specified predictive law is a different
problem. These distinctions prevent overstating the old unit and selector
obstructions as universal impossibility proofs for physics.

## 30. ST8620 — adversarial synthesis and completion audit

The closing investigation cross-checks the principal conclusions against
their strongest counterexamples, scope restrictions and the repository
source map. It is not counted as a new universal theorem.

### What is proved in this campaign

1. Exact positivity and capacity conditions on pair-transition data do not
   suffice for triple compatibility; a strict, reversible, symmetric
   counterexample violates a Boolean event inequality by −2W6/19.
2. Complete lower-order population laws do not universally determine higher
   laws. Positive rational constructions give equality through any finite
   cutoff and a nonzero next-order difference, including on strict and on
   the declared legacy sign cover.
3. Infinite ideal data can identify a supplied event-law representation,
   modulo idle conventions; the independent-plus-common-event class has
   the uniqueness statement in round 24. No universal existence theorem
   for arbitrary FIN dynamics is claimed.
4. Explicit recursive source assumptions can determine and approximate an
   entire hierarchy. Different valid innovation assumptions still give
   different predictions. Stationarity does not supply time evolution.
5. Positive common-clock modifications preserving two distinct strict modes
   are rigid, but that does not select the base composition or physical
   clock units. Finite precision and one-mode fits have explicit limitations.
6. An open quantum version of finite-hierarchy nonidentification exists,
   while exact unitary marginal constraints are an important extremal
   exception. Classicality was not silently assumed to be a proof about
   every quantum completion.
7. Fully specified relational protocols produce dimensionless predictions
   without choosing an absolute unit; dynamics, preparation and readout
   ambiguities must be classified separately.

### What was refuted or narrowed

- Gram positivity as a sufficient hierarchy-extension certificate.
- Pair-law closure as a general implication of strict spectral structure.
- Unspecified fractality, stationary entropy or maximum aggregate activity
  as unique microscopic source principles.
- Equality of averaged instantaneous generators as equality of dynamics.
- The blanket interpretation that *every* finite moment vector or every
  finite quantum marginal constraint is nonidentifying.
- Treating every parameter freedom as physically meaningful: idle attempts
  and common unit changes can be invisible conventions for stated records.
- Treating the absence of an absolute SI scale or a uniquely realized
  initial state as a universal proof against physical predictivity.

### What remains conditional

The state space of constituents, the actual coupling/evolution rule, its
allowed event class, the innovation and time laws, physical encoding of
preparations/measurements, and the law connecting different spatial scales
remain additional inputs. No strict derivation of them was found in the
inspected source chain. This is not a claim that no future theorem can
derive them. None of the new objects completes the canonical legacy-to-
strict map, transfers the historical EW/EM/gravity identities, discharges
QW-2191, creates spacetime, or derives the Standard Model or gravity.

### Most important result and deepest surviving interpretation

**A kernel can fix rich finite spectral geometry while leaving genuinely
different laws of collective change.** This is not merely an eigenvalue-
versus-eigenvector ambiguity: the examples retain the entire kernel and
arbitrarily many complete lower-order dynamics. In the declared classical
class, the missing law is naturally represented by independent transition
rates plus an intensity measure on stochastic transformations, not by
another scalar extracted from A.

That representation is not asserted to be the unique mathematical nature
of the whole FIN repository. The deepest interpretation supported by this
campaign is a finite spectral-information framework with partially
specified composition and operational laws. The exact spectral duality
is real; a unique physical completion does not follow from it. The negative
results are class-scoped, and the positive recursive and relational results
show why a finitely specified predictive extension is still mathematically
possible after its premises are honestly supplied.

### Highest-information next direction

Take one actual, explicitly sourced FIN update/coupling equation—not a
kernel fit—and derive its state space, event or quantum-instrument law,
temporal semantics and observable higher-order prediction. It must decide
between the concrete matched-pair alternatives in this report without
fitting the distinguishing third/fourth-order outcome. If the update has
memory, derive the memory law rather than forcing a Markov representation.
If it is self-similar, prove the innovation rule and a residual-to-prediction
bound of the form in round 28. A theory may allow states and calibrations
as experimental inputs; it must specify how its outputs depend on them.

This is more decisive than generating more kernel fingerprints. Success
would remove an explicit freedom that changes predictions. Failure would
produce a precisely scoped source obstruction. Actual agreement with the
world would still require empirical confrontation; local computations are
not substituted for experiments or external custody.

### Scope of repository use and verification

The campaign used the repository as an indexed research record and read
the current source/guardrail chain, including all six mandated FAR notes,
relevant later source and clock guards, previous strict/legacy proofs,
coupling diagrams, adaptive and neural claims, and the current composition
frontier. It did **not** rerun every historical script or verify every old
referee assertion. That would be a different and much larger audit. Claims
made here are supported by the proofs and checks in this campaign, not by
an invented count of inherited “Proven” labels.

| Repository lane considered | Source evidence used | Consequence for this campaign |
|---|---|---|
| Strict/legacy genealogy | K1/K2, transformation diagram, F2 and S2 | Both branches are retained and tested separately; no silent physical-role transfer. |
| Spectral, Green and variational structure | Current compendium, preceding proofs, strict Laplacian reconstruction | Exact finite calculi remain valid, but they do not source all-order composition. |
| Adaptive and neural dynamics | Actual neural script, adaptive-law passages and corrected mirror report | The activation, update, memory and innovation laws must be specified; static weights are insufficient. |
| Fractal refinement and internal observers | Compendium's refinement/relational-scale sections and current guards | Population replication is not identified with spatial refinement; normalized Green records provide one conditional relational protocol. |
| Selector and dimensional sources | Current AGENTS source cuts and SUMMARY_GROK | No new deterministic selector or unit-bearing source is inferred. Operational gauge freedoms are distinguished from missing predictions. |
| Action, RG and gravity | Actual A1/A4/A8 source files, inspected after the six required notes | Their assumed fields, symbolic shell formulas and open GR obligations do not constitute a sourced event law or a closed physical theory. |
| Quantum and measurement structure | Existing dual-dynamics scope plus the explicit new channel/preparation tests | Open quantum completion is also nonunique; unitary extremality and readout calibration are treated separately. |
| Laboratory/physical claims | Compendium's operational bundle and current custody guards | Local simulations and tests remain mathematical evidence, not an admitted experimental record. |

In particular, A1/A4/A8 serialize declared ansatz and scope summaries.
Their code explicitly retains assumed field content and foundational open
obligations; merely executing those writers would not prove an action,
renormalization law or gravitational source. No interpretation-ranking
percentage from an older summary is used as a calibrated probability of
physical success here.

Detailed requirement/evidence mapping, executed test output and exact replay
comparisons are retained in `COMPLETION_AUDIT.md` and `verification.json`.
The audit distinguishes analytic proofs, exact rational identities and
floating applications. Existing unrelated work was preserved. No PDF,
laboratory record, external audit, DOI or publication upload was generated
for this goal.
