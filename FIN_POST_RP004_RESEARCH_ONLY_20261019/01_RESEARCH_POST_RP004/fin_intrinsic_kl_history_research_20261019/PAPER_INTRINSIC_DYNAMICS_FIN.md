# Can informational self-consistency select FIN dynamics?
## A Fisher–Rao action obstruction and a spectral incompatibility certificate on the strict twelve-state core

**Research note — 19 October 2026.**  **Scope:** mathematical physics / intrinsic-source question F01. **Status:** original, self-contained conditional proofs and bounded numerical checks; not peer reviewed, not a physical validation, not a FIN law adoption.  The authorial term “new” means new *in this note* relative to the inspected GEN-654 proposal, **not** priority over published information-geometry mathematics.

### Abstract

The FIN programme proposes that observable physics emerges from an intrinsically evolving informational nadsoliton. Historical proposals GEN-532–554 and GEN-650–654 sought a least-change relational update and an action assembled from the Kullback–Leibler (KL) divergence between successive states. We give these suggestions one precise and deliberately narrow mathematical interpretation: a positive twelve-outcome probability simplex, the discrete KL trajectory cost, and its correctly normalized continuum limit. We prove (i) the unnormalized continuum cost vanishes, while a nontrivial scaling produces Fisher–Rao kinetic energy; (ii) with only an initial condition, least Fisher–Rao action uniquely freezes the state; (iii) with fixed endpoints, minimizers are explicitly determined great-circle arcs and are invariant under time reversal; and (iv) every nonstationary strict FIN heat trajectory is **strictly nonminimal** for that kinetic action between its own endpoints. A separate test of the Fisher–Rao natural gradient flow of KL relative to the uniform stationary state fails to recover the strict FIN heat generator even after constant clock rescaling: its nonzero linear relaxation eigenvalues are all equal, unlike the six distinct positive strict eigenvalues. A twelve-state diagnostic computes the corresponding action gap. This is a **class obstruction**, not an impossibility theorem for the full nadsoliton ontology. The missing source is a non-circular kinetic/bracket or history principle plus any relevant state and clock mappings.

### 1. Source discipline

The primordial informational nadsoliton is **not** defined here as an underlying stochastic medium; the simplex is a *diagnostic model-state interpretation*, not a claim that the fundamental nadsoliton consists of twelve classical probabilities. The strict twelve-site radial kernel is distinguished from the independently supplied **A7** interacting mean-field/heat-bath model and from the legacy coupled `Phi–Psi` 1+1 action. No conversion among them is assumed.

The historical `GEN-654` expression `S_FIN=∫D(P_t||P_{t+dt})dt` does not itself define a nontrivial continuum functional: the differential, time partition, and normalization must be declared. The present work tests a **comparison completion** of that proposal, not an accepted source-derived FIN action.

Prior guards in `AGENTS(1).md` already rule out inferring a generator merely from static symmetry and rule out zero-cost identity successors in least-change laws. We *reuse* these facts. The substantive local questions here are exactly how a normalized KL trajectory action behaves, whether its extremals contain the strict heat dynamics, and whether the entropy's Fisher–Rao gradient selects that generator.

### 2. Frozen model and definitions

Let `Δ°_11={p∈R^12: p_i>0, sum_i p_i=1}`. For a smooth strictly positive curve `p:[0,T]→Δ°_11`, let

\[
 D(p\|q)=\sum_{i=0}^{11}p_i\log\frac{p_i}{q_i},\qquad
 S_m[p]=\sum_{j=0}^{m-1}\frac{D(p(t_j)\|p(t_{j+1}))}{\Delta},\quad \Delta=T/m.
\]

Strict kernel weights (off-diagonal and at cyclic distances `d=1,…,6`) are

`w=(0.46998567264502017, 0.19204355169010282, 0.09142861427792497, 0.04702916874565040, 0.02413122336363006, 0.011070817321442113)`.

Define `W_ii=0`, `W_ij=w_min(|i-j|,12-|i-j|)`, `s=sum_j W_ij` and `A=sI-W`. The source documents identify `A` as the **strict** selfadjoint Laplacian. Positive weights ensure a connected complete graph with strictly positive off-diagonal rates under the ***conditional*** Markov interpretation `Q=-A`. `π=(1/12,…,1/12)` is the uniform equilibrium. Model time `t` has no justified physical SI calibration.

### 3. Theorem I — normalized KL histories and the correct continuum action

**Theorem I (finite-positive-path limit).** If `p∈C^3([0,T],Δ°_11)` with `min_{t,i}p_i(t)≥c>0`, then

\[
 S_m[p]=\frac12\int_0^T\sum_i\frac{\dot p_i(t)^2}{p_i(t)}\,dt+O(\Delta).
\]

Without the factor `1/Δ` in each summand, `sum_j D(p(t_j)||p(t_{j+1}))=O(Δ)→0`. The literal expression `∫_0^T D(p_t||p_{t+Δ})dt` is `O(Δ²)→0`.

**Proof.** Taylor-expand `h_i=p_i(t+Δ)-p_i(t)=Δ\dot p_i+O(Δ²)`. Since `sum_i h_i=0`,

\[
D(p\|p+h)=-\sum_i p_i\log(1+h_i/p_i)
 =\frac12\sum_i\frac{h_i^2}{p_i}+O(\|h\|_\infty^3),
\]

uniformly under the positive lower bound. Thus `D/Δ = Δ/2 Σ dot p_i²/p_i + O(Δ²)`. Summation and Riemann convergence establish the statement. **End.**

The right-hand integral is the Fisher–Rao kinetic energy. It depends on the supplied curve parameter `t`: the normalization produces a nontrivial functional but does not derive a clock from the kernel.

### 4. Theorem II — exact minimizers and the frozen-state obstruction

Set `u_i(t)=2 sqrt(p_i(t))`. Then `sum_i u_i²=4`, `u_i>0`, and

\[
 S[p]=\frac12\int_0^T\|\dot u(t)\|_2^2dt.
\]

**Theorem II(a) (free endpoint).** With `p(0)=p_0∈Δ°_11` and `p(T)` unconstrained, the unique global minimizer is `p(t)=p_0` for every `t`. Its action is exactly zero.

**Proof.** The integrand is nonnegative. The constant path attains zero. Any minimizer has `dot u=0` almost everywhere, hence `p` constant by absolute continuity. **End.**

**Theorem II(b) (two fixed endpoints).** For `p(0)=p_0`, `p(T)=p_1` strictly positive, let `B=Σ sqrt(p_{0i}p_{1i})` and `θ=arccos(B)∈[0,π/2)`. The minimum is

\[
 \boxed{S_\min(p_0,p_1;T)=\frac{2\theta^2}{T}},\qquad
 \sqrt{p_i(t)}=\frac{\sin((1-t/T)\theta)}{\sin\theta}\sqrt{p_{0i}}+
 \frac{\sin((t/T)\theta)}{\sin\theta}\sqrt{p_{1i}},
\]

with the constant-path limit when `θ=0`. In the strictly positive orthant the minimizing short great-circle segment is unique. Reversing time and swapping endpoints preserves `S`.

**Proof.** The `u` map embeds the Fisher metric isometrically into the positive orthant of a Euclidean sphere of radius two. The spherical distance between endpoint vectors equals `2θ`. Cauchy–Schwarz gives `S≥distance²/(2T)=2θ²/T`; constant-speed shortest spherical arcs attain equality. The written sine interpolation is the standard spherical geodesic and has strictly positive entries. Since the integrand is quadratic in velocities, reversal leaves it invariant. **End.**

This excludes autonomous time orientation from this action **by itself**. Fixed nonidentical endpoints can force nontrivial motion, but the second endpoint is then an extra preparation/boundary datum, not generated by the starting relational state.

### 5. Theorem III — strict FIN heat is not selected by Fisher–Rao least action

Adopt, for comparison only, the **existing** strict heat semigroup `p(t)=exp(-tA)p_0`, where `p_0∈Δ°_11` and `p_0 ≠ π`.

**Theorem III (strict nongeodesicity for all finite windows).** For each `T>0`,

\[
 S[p_{\rm heat}|_{[0,T]}]
 > S_\min(p_0,e^{-TA}p_0;T).
\]

**Proof.** The heat flow is analytic and strictly positive on `(0,∞)` because `Q=-A` is an irreducible finite-state Markov generator; `p(t)→π`, `dot p(t)→0` exponentially because `A` has a positive gap. Define the nonnegative Fisher speed squared `v(t)²=Σ dot p_i(t)²/p_i(t)`, real analytic on `(0,∞)`, with `lim_{t→∞}v(t)²=0`. Every energy-minimizing geodesic between unequal endpoints must have **constant speed**. If the heat path minimized on any `(0,T)`, `v(t)²` would be constant there, hence (identity theorem for real-analytic functions) constant on `(0,∞)` and therefore identically zero. Then `A p(t)=0`, forcing `p_0=π`, contradiction. **End.**

The statement is exact for the strict positive-weight generator; the numeric check below merely illustrates its strict inequality. It does not prove that the physical nadsoliton follows strict heat, nor does it rule out dissipative laws built using ***different*** state-dependent transport metrics and an explicit thermodynamic potential.

### 6. Theorem IV — a Fisher–Rao KL gradient flow has the wrong FIN spectrum

A different attempt to obtain nontrivial dynamics from information alone is to postulate the Fisher–Rao gradient flow of the relative entropy to uniform equilibrium:

`E(p)=D(p||π)` and `dot p=-grad_FR E(p)`. In full coordinates this gives

\[
 \boxed{\dot p_i=-p_i\left[\log\frac{p_i}{\pi_i}-\sum_j p_j\log\frac{p_j}{\pi_j}\right] }.
\]

**Theorem IV (strict spectral mismatch).** At `π`, the linearization of this normalized Fisher–Rao entropy-gradient flow on the eleven-dimensional tangent space `sum_i δp_i=0` is `dot x=-x`. Under any constant positive rescaling of its comparison clock it remains `dot x=-c x`, with all eleven decay rates equal. In contrast, strict heat has `dot x=-A x` and positive relaxation eigenvalues

`0.754121154207080 (×2), 1.577049514427610 (×2), 1.961406861976446 (×2), 2.199568849333210 (×2), 2.298606272079097 (×2), 2.342182041146301 (×1)`.

Therefore **no single positive clock scale** makes the two dynamics identical in an open neighbourhood of equilibrium.

**Proof.** The Fisher metric is `g_p(v,w)=Σ_i v_iw_i/p_i` on tangent vectors of zero sum. Projecting `dE` to this metric yields the displayed equation. For `p=π+x` with `sum x=0`, expand `log(p_i/π_i)=x_i/π_i+O(|x|²)`; the scalar average is `O(|x|²)`, so `dot x_i=-x_i+O(|x|²)`. Connected strictly nonuniform radial weights generate the listed nonconstant positive eigenvalues of `A`, and a constant rescaling cannot turn this spectrum into one repeated eigenvalue. **End.**

**Important distinction:** Theorem II minimizes *kinetic energy of an entire history*. Theorem IV assumes a *first-order gradient-flow prescription* for entropy. They are **different** attempted completions of GEN-654, and neither is a forced consequence of the original ontological premise.

### 7. Diagnostic verification on the published strict weights

For `p0=(1-0.35)π+0.35 e_0` and `T=1` (arbitrary dimensionless comparison clock), a float64 spectral/matrix-exponential calculation gives:

| Quantity | Value |
|:--|--:|
| strict row sum `s` | 1.660307278766099 |
| strict gap | 0.754121154207079 |
| heat-history Fisher action `S_heat` | 0.243273129725079 |
| absolute geodesic minimum with heat's exact endpoints | 0.186642154969417 |
| ratio `S_heat/S_min` | 1.303420064802297 |

A discrete KL replay with `m=8,16,32,64,128,256` steps approaches the continuous action with decreasing small-step error asymptotically (not monotonically for every step size). The script `verify.py` reproduces the values, checks the kernel, and independently approximates the Fisher-gradient tangent derivative. These floating computations are not a numerical proof of the all-data analytic theorems.

### 8. Control: even equal stationary state and total attempt/escape rate do not select path laws

For a *diagnostic comparison* set of generators only, retain the strict radial graph but change distance-one rates `w_1→(1+α)w_1` and distance-two rates `w_2→w_2−α w_1`, with `α=0.1`, leaving all other rates fixed. This preserves positivity, symmetry, translation/reflection invariance, uniform equilibrium **and the complete row escape rate** because each distance appears twice. Nonetheless the generator and its temporally resolved path law change. The relative entropy rate from the original stationary CTMC to the comparison stationary CTMC is exactly

\[
 \mathcal K(Q\|Q')=\sum_i\pi_i\sum_{j\ne i}
 \left(q_{ij}\log\frac{q_{ij}}{q'_{ij}}-q_{ij}+q'_{ij}\right)>0
\]

and is ≈`0.01821605667` per dimensionless comparison time in this case. **This is not offered as a novel FIN nonuniqueness discovery:** it is a quantitative sanity control agreeing with accepted earlier activity/source obstructions. Neither comparison generator is being adopted as a fundamental law.

### 9. Constructive contrast — a graph-adapted information gradient *can* reproduce strict heat, but does not source it

For all `a,b>0` define the logarithmic mean `L(a,b)=(a-b)/(log a-log b)` if `a!=b` and `L(a,a)=a`. Define a positive state-dependent quadratic mobility operator on tangent vectors:

\[
 (K_p\xi)_i=\sum_{j\ne i}W_{ij}L(p_i,p_j)(\xi_i-\xi_j).
\]

**Theorem V (exact graph transport representation, reused construction).** For `E(p)=D(p||π)` and the strict graph weights,

\[
 \boxed{\dot p=-K_p\nabla E(p)=-A p}.
\]

**Proof.** Since `∂E/∂p_i=log(p_i/π_i)+1`, the gradient difference on edge `(i,j)` is `log p_i-log p_j`. The defining logarithmic-mean identity yields `L(p_i,p_j)(log p_i-log p_j)=p_i-p_j`. Hence `(K_p∇E)_i=Σ_jW_ij(p_i-p_j)=(Ap)_i`. Each quadratic form `ξ^T K_p ξ=Σ_{i<j}W_ij L(p_i,p_j)(ξ_i-ξ_j)^2` is nonnegative. **End.**

This constructive result shows that information geometry need not *prevent* the desired dynamics. But `K_p` is precisely the **additional graph-adapted kinetic geometry / mobility object**: imposing its logarithmic-mean form in order to get `-A` is an exact representation of the *selected heat law*, **not a derivation of its necessity** from KL alone. This distinction is consistent with established discrete transport-gradient research (Maas 2011) and the earlier FIN activity-source obstructions.

**Additional local no-go (static potential):** If one appends any time-independent smooth potential `V_A(p)` to the Fisher kinetic Lagrangian, then with `u_i=2√p_i` and `||u||=2`, the Euler–Lagrange equations have the constrained second-order form `ddot u= -Π_u∇_u V_A - (||dot u||²/4)u` (sign depends on whether `L=kinetic-V_A` is chosen as here). A single initial point `u(0)` still leaves the tangent initial velocity undetermined. The equation and kinetic term are reversible under `t→-t`. Thus merely adding a static kernel-dependent potential does not establish a unique autonomous first-order transition law, a dissipative arrow, or an intrinsic comparison clock.

### 10. Interpretation and shortest legitimate next research question

The current source data admit an intrinsically defined **symmetric information geometry** on positive distributions. This geometry is useful but insufficient to generate a nontrivial dynamical arrow or to reproduce the strict heat operator:

- Least-change action with one initial state selects stasis (Theorem II).
- Endpoint-forced Fisher geodesics do not yield the strict heat semigroup (Theorem III).
- Fisher natural gradient of Shannon KL yields dynamics but has the wrong strictly-resolved relaxation spectrum (Theorem IV).
- The initial kernel and stationary state do not determine path-space activity, even with a fixed total escape clock (comparison in §8).

A successful next candidate would need a **new typed object already supported independently by FIN**: for example a microscopic mobility/current law, a phase/symplectic bracket, or a mathematically defined relational history/composition principle. Merely declaring a new mobilty, adding a forward entropy gradient, or specifying the desired two-time law is not a derivation. A candidate must exclude at least one current same-static/same-clock countermodel *without using its selected generator as hidden input*.

Nothing here derives emergent spacetime, the legacy action, quantum observables, real SI time, photons, matter or gravity. The primordial-nadsoliton ontology remains a guiding **hypothesis**, not an implication of the Fisher calculation.

### 11. References and provenance

**FIN inputs inspected:** supplied `AGENTS(1).md` (especially kernel split, strict state selection, activity-source and FIN intrinsic research guardrails), `GEN001-gen956(1).md` (GEN-521–555, GEN-649–654, GEN-929–944), earlier `FIN_missing_objects_research_20261019/INTRINSIC_HISTORY_PROVIDER_OR_OPEN_LEMMA.md`; FIN Release 10.88 finite spectral core. File hashes in `INPUTS.sha256`.

**Related external scholarship (not evidence for FIN ontology):** Jan Maas, *Gradient flows of the entropy for finite Markov chains*, Journal of Functional Analysis **261** (2011), 2250–2292, DOI 10.1016/j.jfa.2011.06.009; Fischer–Rao simplex square-root geometry and local KL expansion are standard information geometry; see e.g. Byrne (2013), *Geodesic Monte Carlo on Embedded Manifolds*, Scandinavian Journal of Statistics. Maas uses a **different graph-adapted transport metric**, demonstrating that recovery of Markov heat as an entropy gradient flow depends on kinetic geometry, not just the Fisher metric.

**Scope classification:** Theorems I–IV = proved here under stated finite-positive-probability assumptions (`CONDITIONAL_MATHEMATICAL`). Equal-activity path KL diagnostic = bounded control with standard likelihood identity (`REUSED_OBSTRUCTION_QUANTIFIED`). Fundamental source F01, native kinetic source F03, spacetime/composition F06 = `OPEN_SOURCE`. No physical or empirical claim.
