# ST8548 — Same FIN heat flow, different fluctuation actions

Krzysztof Żuchowski  
Independent Researcher — Fractal Information Theory Project  
ORCID: 0009-0002-0909-3613  
6 September 2026

## Research outcome

The quadratic logarithmic-mean representation in ST8547 does not equal the
fluctuation action of independent Markov walkers with the same strict rates W.
Both give the same mean heat equation and the same local time-reversal
asymmetry. Their variances disagree away from equilibrium, and their fourth
cumulants disagree even at equilibrium.

The exact jump action contracts under summing child-resolved currents even
for correlated fine states. This is stronger than the corresponding
logarithmic-mean quadratic-form compatibility, whose equality case in ST8547
requires a product state. This distinction sharpens the interpretation of
coarse observational closure in FIN.

The result is conditional on a declared microscopic population model. The
nonquadratic gradient structure is established literature, with relevant work
by Mielke, Peletier and Renger and by Liero and collaborators cited below.
The contribution here is its derivation, numerical application and falsification
consequences for the current FIN programme, not a claim of new universal physics.

## 1. Inputs and microscopic assumptions

Use the positive strict conductances on twelve vertices,

\[
W_{ij}=\frac{\cos(0.18575d+0.16250)}{1+d^{9/5}},\quad
d=\min(|i-j|,12-|i-j|),\quad i\ne j,\qquad W_{ii}=0.
\]

Let Q=W-diag(W1)=-A be a row generator. Take N independent continuous-time
Markov walkers with these rates; p_i is the empirical occupation fraction.
At state p, a jump i to j changes the empirical measure by
(e_j-e_i)/N and has total rate N p_i W_ij.

The exact finite-N density generator acting on a function F is

\[
(\mathcal G_NF)(p)=N\sum_{i\ne j}p_iW_{ij}
[F(p+(e_j-e_i)/N)-F(p)]. \tag{1}
\]

Applying this to coordinate functions gives mean evolution dot p=-Ap.
The additional premises are a population of walkers, their independence and
the interpretation of W as jump rates. The strict spectral theorem alone
does not provide these premises. No physical clock is assigned.

The raw signed legacy kernel does not directly define these rates. No legacy
completion or role transfer is used.

## 2. Exact local current Hamiltonian

Orient an edge i to j and write

\[
a=W_{ij}p_i,\qquad b=W_{ij}p_j,\qquad a,b>0.
\]

The local Hamiltonian for normalized net current is

\[
h(s;a,b)=a(e^s-1)+b(e^{-s}-1). \tag{2}
\]

**Proof.** The two oriented events have generator rates Na and Nb and change
the count by +1 and -1. Applying the nonlinear transform
N^{-1}e^{-Nf}\mathcal G_N e^{Nf} to a linear counting tilt gives (2).
It is also the cumulant generating function per unit exposure of two frozen
Poisson intensities. QED.

For a closed Markov chain (2) is a local/infinite-population action density.
It is not the exact distribution of counts over an arbitrary finite time
window, since occupations change and returns are possible. We independently
check finite-time counts using a tilted full generator in section 8.

Equation (2) gives the instantaneous cumulant rates

\[
\kappa_1=a-b,\quad\kappa_2=a+b,\quad
\kappa_3=a-b,\quad\kappa_4=a+b. \tag{3}
\]

For N independent walkers cumulants of total counts scale by N. Conditional
short-window statements carry the usual leading factor N dt.

## 3. The exact convex fluctuation cost

The Legendre transform of (2) is

\[
\begin{aligned}
\ell(j;a,b)
&=j\,\operatorname{arsinh}\frac{j}{2\sqrt{ab}}
-\frac j2\log(a/b)
-\sqrt{j^2+4ab}+a+b. \tag{4}
\end{aligned}
\]

**Proof.** The maximizing tilt solves
j=a e^s-b e^{-s}; hence

\[
s_*=\operatorname{arsinh}\frac{j}{2\sqrt{ab}}-\frac12\log(a/b).
\]

At that value, a e^s+b e^{-s}=sqrt(j^2+4ab). Substitution gives (4). QED.

An independent equivalent expression is

\[
\ell(j;a,b)=\inf_{u,v\ge0:\,u-v=j}
[\mathfrak h(u|a)+\mathfrak h(v|b)],
\quad
\mathfrak h(u|a)=u\log(u/a)-u+a. \tag{5}
\]

The minimizers are u=(sqrt(j^2+4ab)+j)/2 and
v=(sqrt(j^2+4ab)-j)/2. This proves joint convexity in (j,a,b), by
convexity of relative-entropy perspectives and contraction under a linear
constraint. The cost also obeys
ell(cj;ca,cb)=c ell(j;a,b), c>0.

There is a unique zero at j0=a-b. Strict convexity in j follows from

\[
\partial_j^2\ell(j;a,b)=\frac1{\sqrt{j^2+4ab}}. \tag{6}
\]

## 4. Cosh dissipation and entropy

Put g=2sqrt(ab) and z=(1/2)log(a/b). Define the symmetric potentials

\[
\mathcal R(j)=j\,\operatorname{arsinh}(j/g)-\sqrt{j^2+g^2}+g,\qquad
\mathcal R^*(\zeta)=g(\cosh\zeta-1). \tag{7}
\]

They are Legendre dual and

\[
\ell(j;a,b)=\mathcal R(j)+\mathcal R^*(z)-jz,\qquad
j_0=\partial_\zeta\mathcal R^*(z)=a-b. \tag{8}
\]

The factor one-half in z corresponds to one-half of the negative relative
Shannon entropy in the energy-dissipation identity. It is a normalization
forced by the path-reversal formula, not a new physical constant.

For independent stationary walkers, the empirical occupation probability is
multinomial:

\[
\Pr(Np)=\frac{N!}{\prod_i(Np_i)!}\prod_i\pi_i^{Np_i}.
\]

Stirling's formula yields minus N^{-1}log probability tending to D(p||pi).
Here pi is uniform because the supplied W is symmetric. Thus independence
of microscopic trials supplies a conditional statistical origin for the
relative entropy and the nonquadratic dissipation together.

## 5. Why the quadratic Onsager model has the same mean but different noise

Set m=Lambda(a,b)=(a-b)/log(a/b). The quadratic action with the same mean and
the same antisymmetric cost difference is

\[
\ell_{\rm Q}(j;a,b)=\frac{[j-(a-b)]^2}{4m}. \tag{9}
\]

Indeed, both (4) and (9) satisfy

\[
\ell(j)-\ell(-j)=-j\log(a/b),\qquad \ell(a-b)=0. \tag{10}
\]

The normalization in (9) is explicit. It makes the corresponding symmetric
dissipation j^2/(4m) pair with one-half of negative Shannon entropy; equivalently
it rescales the ST8547 Onsager convention by the compensating entropy factor.

However their local variances are different:

\[
v_{\rm jump}=a+b,\qquad v_{\rm Q}=2\Lambda(a,b). \tag{11}
\]

**Theorem 1 — fluctuation mismatch.** For all a,b>0,

\[
\frac{v_{\rm jump}}{v_{\rm Q}}
=z\coth z\ge1,\qquad z=\frac12\log(a/b), \tag{12}
\]

with equality exactly when a=b (continuous value at z=0).

**Proof.** Equation (6) at j0 gives inverse curvature a+b.
The quadratic inverse curvature is 2m. Substitute a=sqrt(ab)e^z,
b=sqrt(ab)e^-z to obtain z coth z. The inequality follows from
|tanh z|<=|z|, strictly away from zero. QED.

Near equilibrium the ratio is 1+z^2/3+O(z^4). At equilibrium the quadratic
variance is correct, but the exact fourth cumulant rate is 2a>0 whereas a
quadratic current law has zero fourth cumulant rate.

The Onsager identity in ST8547 remains true. Its deterministic metric cannot
be reinterpreted as a complete microscopic noise law without this extra test.
The present result refutes that reinterpretation, not the gradient identity.

**Theorem 1b — an integer-jump consistency bound.**
For any Markov additive count with integer jump sizes k and nonnegative local
rates c_k (finite first and second moments),
mu=sum k c_k and v=sum k^2 c_k obey v>=|mu|.
This follows from k^2>=|k| and the triangle inequality.
But the quadratic prescription gives

\[
\frac{v_{\rm Q}}{|j_0|}=\frac2{|\log(a/b)|}<1
\quad\text{if }|\log(a/b)|>2. \tag{12b}
\]

In that regime it cannot be the exact local law of any such integer-valued
Markov count, even allowing collective integer jumps. Continuous Gaussian
noise and coarse diffusion approximations are outside this impossibility
statement. In the strict example a/b=1/12, the jump ratio is 13/11, whereas
the quadratic ratio is 2/log(12), about 0.805.

### Three different means with distinct roles

- The logarithmic mean gives j0=Lambda(a,b) log(a/b).
- The geometric mean gives the prefactor 2sqrt(ab) in the cosh potential.
- The arithmetic combination a+b gives the instantaneous shot-noise variance.

These coincide in the needed normalizations close to equilibrium; away from
equilibrium they serve different mathematical purposes.

### Whole-graph covariance ordering

The occupation covariance rate of independent walkers is

\[
\Sigma(p)=\sum_{i<j}W_{ij}(p_i+p_j)(e_i-e_j)(e_i-e_j)^\top.
\]

The quadratic Onsager noise convention would use 2K_Lambda(p).
By arithmetic-logarithmic mean ordering,

\[
\boxed{\Sigma(p)-2K_\Lambda(p)\succeq0.} \tag{13}
\]

For connected positive support, equality of the matrices holds precisely at
uniform p. This follows edgewise, since each difference conductance is
nonnegative and vanishes iff p_i=p_j.

## 6. Exact fluctuation contraction through hidden fibers

Suppose one coarse edge represents several fine channels with positive
rates a_b,b_b. Set a=sum a_b, b=sum b_b.

**Theorem 2 — exact current contraction.**

\[
\boxed{
\inf_{\sum_b j_b=j}\sum_b\ell(j_b;a_b,b_b)
=\ell(j;a,b).} \tag{14}
\]

**Proof.** Let s solve a e^s-b e^-s=j and choose
j_b=a_b e^s-b_b e^-s. All fine costs have derivative s at this choice.
Convexity makes it a global constrained minimizer. The Legendre values sum
to sj-sum_b h(s;a_b,b_b)=sj-h(s;a,b). QED.

No condition a_b/a=b_b/b is needed. In the product FIN refinement the coarse
rates are a=W_ij sum_b P_ib and b=W_ij sum_b P_jb. Thus local current
fluctuation contraction holds even for base-child correlations.

The infinitesimal covariance likewise satisfies
C Sigma_fine(P) C^T=Sigma_coarse(CP), because its edge weights are linear
in the fine probabilities and vertical increments vanish after C.

For the full product generator A tensor I+I tensor B, projection of the base
coordinate is an exactly lumpable Markov process for every initial fine
state. This is not the interacting Ising/prism model used in older stationary
RG studies; the generators and lumpability hypotheses are different.

By contrast the quadratic action contracts to

\[
\inf_{\sum j_b=j}\sum_b\ell_{\rm Q}(j_b;a_b,b_b)
=\frac{[j-(a-b)]^2}{4\sum_b\Lambda(a_b,b_b)}. \tag{15}
\]

Since sum Lambda(a_b,b_b)<=Lambda(a,b), it generally differs from the coarse
quadratic action. Away from j=j0, equality occurs precisely when the two
normalized fine-channel profiles coincide. This is the metric defect
classified in ST8547.

Consequently preserving the exact coarse stochastic process does not require
preservation of that particular quadratic Onsager representation.
The entropy-production KL gap in ST8547 remains valid; it cannot be equated
with an unavoidable defect of the exact coarse counting statistics.

## 7. Concrete strict FIN results

Use p_i=(i+1)/78, edge 0 to 11, and
W_0,11=0.4699856726450201. Then

| Quantity | Value |
|---|---:|
| a | 0.006025457341602822 |
| b | 0.07230548809923387 |
| mean current a-b | -0.06628003075763104 |
| exact local variance rate | 0.0783309454408367 |
| quadratic variance rate | 0.05334609311241991 |
| exact/quadratic ratio | 1.4683539294201824 |

For two child channels with forward proportions (0.8,0.2), backward
proportions (0.25,0.75), and j=0.04, the exact contraction error is
2.78e-17. The quadratic fine-contracted cost is 0.1178653956, compared with
the coarse quadratic cost 0.1058694675. These are floating evaluations of
the exact formulas, not interval-certified constants.

The graph covariance projection error for a correlated 24-state model is
1.11e-16; its quadratic covariance projection defect has norm 0.0079424221.

## 8. Finite-time verification and empirical interpretation

Tests use the full 12-state tilted Markov generator, whose off-diagonal entries
on the counted edge receive exp(+s) and exp(-s), while the diagonal remains
unchanged. Derivatives of exp(t Q_s) are computed by a block-triangular
matrix exponential. This is independent of the local Poisson approximation.

For a single walker initially distributed as p, the exact finite-time variance
per time approaches 0.0783309454 as t decreases:

| t | variance/time | fourth cumulant/time |
|---:|---:|---:|
| 0.1 | 0.0745294243 | 0.0717795571 |
| 0.01 | 0.0779222577 | 0.0776108767 |
| 0.0001 | 0.0783268250 | 0.0783236669 |
| 0.00001 | 0.0783305334 | 0.0783302175 |

The local quadratic variance would be 0.0533460931 and its local fourth
cumulant zero. Finite-time cumulants of a state-dependent Gaussian model need
not stay Gaussian; comparisons concern its local generator and the shrinking
window limit.

A possible dimensionless statistic is short-window net-count variance divided
by the absolute mean net count. Away from zero drift it tends to
(a+b)/|a-b| for single jumps, and to 2/|log(a/b)| for the quadratic prescription.
Clock scale and population size cancel in this ratio. At zero drift one can
instead examine the local fourth-to-second cumulant ratio. This requires
repeated prepared states, bidirectional counts and finite-window/detector
controls. It is a proposed conditional experiment specification, not an
existing FIN measurement.

## 9. An explicit correction to an older kinetic counterexample

The three-bit pure-parity model of ST2547–ST2561 has energy
E=-theta x1 x2 x3. Every flip has |Delta E|=2|theta|.
Therefore its three named rate rules are exactly proportional:

\[
Q_{\rm Met}=(1+e^{-2|\theta|})Q_{\rm HB},\qquad
Q_{\rm sym}=2\cosh(\theta)Q_{\rm HB}. \tag{16}
\]

They do not give different normalized spectral shapes in that example.
At theta=0.3 the factors are 1.5488116361 and 2.0906770283, with replay
residuals below 4.45e-16. The earlier statement that this specific example
distinguished more than clock scaling is refuted.

General kinetic nonuniqueness remains true, but needs a graph/state with
multiple energy-difference magnitudes or a nonconstant edge activity
counterexample. Archived files are retained as historical artifacts; this
correction supersedes that interpretation.

## 10. Why the microscopic premise is still necessary

Let walkers instead move in synchronized groups of size s. With N total
particles, each group jump changes empirical mass by s/N; its rate is the
number of groups at the vertex times W. The local Hamiltonian becomes

\[
h_s(\xi)=\frac1s[a(e^{s\xi}-1)+b(e^{-s\xi}-1)].
\]

It has the same drift a-b but variance rate s(a+b). The static empirical
entropy cost per total particle is D/s. Thus the mean FIN heat equation
alone cannot select independence, fluctuation scale or the population model.

This also distinguishes the 12-state walker occupation model from the older
4096-state binary configuration model. Neither model is automatically a
physical realization of FIN, and they must not be identified.

## 11. Source comparison and novelty

Mielke, Peletier and Renger, [On the relation between gradient flows and the
large-deviation principle, with applications to Markov chains and diffusion](https://arxiv.org/abs/1312.7591),
derive generalized gradient structures from empirical Markov processes and
explicitly distinguish the nonquadratic structure from the logarithmic-mean
quadratic structure. See section 4.1 and example 4.2.

Liero, Mielke, Peletier and Renger,
[On microscopic origins of generalized gradient structures](https://doi.org/10.3934/dcdss.2017001),
further develop cosh dissipation and its relation to coarse limits.
The basic microscopic construction is known. This study's role is to repair
and strengthen the FIN interpretation with explicit current and covariance
tests, not to claim priority over those works.

## 12. What has and has not been achieved

Proved, conditional on supplied positive rates and the microscopic model:

- the exact local Poisson current cost and its cosh decomposition;
- same drift and reversal asymmetry but different quadratic fluctuations;
- a matrix ordering of true and quadratic occupation covariances;
- exact contraction of the jump action over hidden child channels;
- the contrast with the ST8547 Onsager correlation defect.

Refuted: the identification of the quadratic Onsager metric with the full
microscopic fluctuation law; and the older pure-parity three-rule
spectral-shape counterexample.

Still conditional: microscopic independence, a physical population, rates,
clock, preparation and the apparatus recording jumps. No Born rule,
unitary measurement dynamics, legacy role transfer, QW-2191, continuum
spacetime, SM/GR, dimensional L_total or ToE closure follows.

The best next research target is a finite-noise validation protocol comparing
second and fourth count cumulants while controlling detector losses and
finite observation windows. More fundamentally, a sourced FIN microscopic
law must specify single events, collective events or another counting law.
That choice is now falsifiable even when the mean operator A is unchanged.
