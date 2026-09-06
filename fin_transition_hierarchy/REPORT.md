# FIN: transition-resolved composition and the hierarchy problem

Krzysztof Żuchowski — Independent Researcher, Fractal Information Theory Project  
ORCID 0009-0002-0909-3613. Research checkpoint, 7 September 2026.

**Active goal: rounds 1–10 of 30. Not the final goal report.**
This is a working mathematical report, not a publication-layout exercise.
No PDF is generated. Exact rational calculations are in `research.py` and
the scientific checks in `test_research.py`. Later rounds must continue from
the unresolved source questions, not repeat the counterexample inventory.

The main new result is stronger than the preceding activity audit: even the
**complete two-constituent generator** does not determine higher collective
dynamics. This is demonstrated for the actual strict kernel in an explicitly
infinitely extendible, reversible, maximally synchronous class. A second
result shows that a nonnegative Gram transition matrix can fail even the
necessary three-constituent event inequalities.

All theorems below are conditional on their stated finite classical model.
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
actual replication-consistency package were consulted; this checkpoint
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

## Current outcome and next checkpoint

The new strict-compatible no-go is not merely that eigenvalues omit
geometry: even complete finite-order collective Markov laws omit higher
composition. Conversely, all ideal orders identify a supplied scalar law,
and a specified self-similar recursion can compress that entire law into
a finite prescription. The unresolved issue is which prescription, if any,
is actually derived from FIN, including its temporal sampling semantics.

Rounds 11–30 remain open. Next investigate whether a stationary recursive
law also determines temporal dynamics; then test stronger composition
classes (in particular common random clocks/subordination) and revisit the
actual source obligations across the repository. Do not mark the active
goal complete at this checkpoint.
