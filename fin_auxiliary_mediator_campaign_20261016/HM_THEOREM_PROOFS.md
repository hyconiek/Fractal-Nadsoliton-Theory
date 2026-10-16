# HM-001–HM-005: new derivations and proof ledger

Date: 2026-10-16. Primary: unchanged FIN q12, LOO heat bath and A7=XX^T. Comparison: AUX-HB-OU-kappa, explicitly **not** the primary process. The following is mathematical model analysis, not a physical mediator source. All limiting claims invoking MT are **conditional** on the admitted exact local-Gibbs/domain estimates. The current checkout/review is not present in this workspace: predecessor material actually examined is bound in INPUT_PROVENANCE.json.

## HM-001: finite-N conditional cavity identity (exact)

Set `x_i ∈ {0,...,11}`, `n=∑_i e_{x_i}`, and `A=XX^T`, with rows `X_j ∈ R^7`, `A_jj=a_0` constant, `g>0`. The genuine primary labelled Gibbs weight is proportional to `exp[g/(2N) n^T A n]`. The MP7 auxiliary joint weight is proportional to `exp[-N||θ||²/(2g)+θ·X^T n]` (for counts multiply by `N!/∏n_j!`). Completing the square gives `θ|n ∼ N((g/N)X^T n,(g/N)I_7)` and the correct Gibbs count marginal. For one chosen old slot, `r=n-e_i`; *not* the full-n conditional, define

`ν_r(dθ) ∝ exp[-N||θ||²/(2g)+θ·X^T r] [∑_j exp(X_j·θ)] dθ`.

Draw fresh θ∼ν_r and then `j` with `p_j(θ)=exp(X_j·θ)/∑_l exp(X_l·θ)`. The joint mass for output j is proportional to

`∫ exp[-N||θ||²/(2g)+θ·X^T(r+e_j)]dθ = C exp[(g/(2N))(r+e_j)^T A(r+e_j)]`.

With constant `A_jj`, this is `C_r exp[(g/N)(Ar)_j]`, hence `Pr(j|r)=softmax_j((g/N)Ar)` EXACTLY. This is the native original-denominator unit-attempt LOO rule. The factor `∑exp(X_j·θ)` within ν_r is essential; using the Gaussian posterior ρ_n or bare Gaussian ρ_r is generally wrong. An old label is discarded *before* the cavity draw. If reflecting at D_N, generate the full proposal and reject if `r+e_j∉D_N`; never renormalize over allowed j. θ is ephemeral: the identity does not imply persistent mediator paths, field clock, source of A7 or an independent heat-bath principle.

## HM-002: genuine persistent diagnostic, detailed balance, and compression rigidity (exact)

Define joint `Π_N(n,dθ)∝Mult(n) exp[-N||θ||²/(2g)+θ·X^T n]dθ`. The AUX-HB-OU-kappa joint generator acts on labelled x,θ:

`L^aux f(x,θ)=∑_{i=1}^N∑_j p_j(θ)[f(x^{i→j},θ)-f(x,θ)] + κ( (g/N)X^T n-θ)·∇_θ f + (κg/N) Δ_θ f`.

For fixed θ, independent slots refresh from the *exact* conditional product law `Π(x|θ)∝∏exp(X_{x_i}·θ)` (symmetric detailed balance). For fixed x, the field OU is reversible against `ρ_n=N((g/N)X^T n,(g/N)I_7)` (symmetric detailed balance). The sum is reversible for Π. Rejecting count jumps exiting D_N preserves the restricted Π_D, and the θ OU remains unchanged: no allowed-label renormalization.

Write `P=E_Π[·|n]`, `Q=I-P` in L²(Π), `m=(1/N)∑ σ_{x_i}`. The actual joint generator satisfies

`L^aux m = φ(θ)-m`, `φ(θ)=∑_j σ_j p_j(θ)`,

on the unrestricted model (away from the reflected boundary in D_N). Since ρ_n has positive Gaussian density on R^7, and along the parity field θ=t e_6, `φ(t e_6)=tanh(b t)` with `b=√(Λ6/12)>0`, φ is not ρ_n-a.s. constant. Therefore

`||QL^aux Pm||²_{L²(Π)} = E_π Var_{ρ_n}[φ(θ)] > 0`.

In a reflected domain, retain any positive-mass *interior* count state; at this state all attempted types remain allowed, so the same strict inequality holds. The reflection transition conventions matter, and a positive interior count type must exist.

For `K(t)=P exp(tL^aux)P` on count functions, Gaussian smoothness and bounded derivatives of softmax make `m∈Dom(L²)` (finite labelled states, OU with Gaussian moments). Self-adjointness gives the exact second-order defect

`<m,[K''(0)-K'(0)^2]m>_π = <m,PL^aux QL^aux Pm>_π = ||QL^aux Pm||²_Π >0`.

Thus K cannot obey the semigroup identity and, in particular, cannot be the native count semigroup. Equality of the *count stationary distribution* is insufficient to equate the *count paths*. This no-go applies to **AUX-HB-OU-kappa** for κ>0, not every conceivable mediator lift. The first instantaneous count derivative usually differs as well, but is not used as the decisive proof. The stationary projection retains memory.

## HM-003: full characteristic calculation, critical soft pair (exact generator, formal scaling coefficients)

Let `p=n/N`, `b=√(Λ6/12)=1/√g_c`, `σ_j=(-1)^j`, `m=σ·p`, `h=b θ_6`; `θ_6` is the X parity column coordinate. The exact mean drift is

`L^aux p_j = softmax_j(Xθ)−p_j`; `L^aux θ=κ(gX^Tp−θ)`;

`L^aux m = φ(θ)−m`; `L^aux h=κ((g/g_c)m−h)`.

For θ-only noise: `d<θ_l,θ_r>/dt=2κg δ_lr/N`. For count jumps: `d<p_j,p_k>/dt=(1/N)∑_i p_i∑_l p_l(θ)(δ_{jl}−δ_{ji})(δ_{kl}−δ_{ki})`; equivalently `(1/N)[diag(p)+diag(s)-p s^T-s p^T]` with `s=softmax(Xθ)`. Count/field cross brackets are **zero** because jumps and Brownian increments have disjoint quadratic variations. In the parity direction,

`d<m>/dt = (2/N)[1−mφ(θ)]`, `d<h>/dt = 2κ(g/g_c)/N`, `d<m,h>/dt=0`.

On the parity axis `Xθ=hσ`, `φ=tanh(h)`; at g=g_c the linearization in (m,h) is `J=[[-1,1],[κ,−κ]]`. Right zero vector (1,1), normalized left zero vector `(κ,1)/(κ+1)`, other eigenvalue `−(1+κ)`. The exact soft chart is `s=(κm+h)/(1+κ)`, with

`L^aux s = rκ [φ(θ)−h+(g/g_c−1)m]`, `rκ=κ/(κ+1)`,

`d<s>/dt=[2κ²(1−mφ)+2κ(g/g_c)]/[(1+κ)²N]`.

Thus projected leading noise is `2rκ/N`, and the candidate gain/cubic drift `rκ[(χλ)z−z³/3]` with `χ=1/g_c`. At κ=1, `r=1/2`; at κ=4, `r=4/5`. Their **static** law can coincide despite these different dynamic coefficients.

The complete 18-dimensional tangent space consists of 11 independent occupation coordinates and seven θ coordinates. Four count modes in ker X^T have rate −1. Each of the six active transverse Fourier coordinates produces a two-by-two stability block with determinant `κ(1−Λ_k/Λ6)>0` and trace `−(1+κ)`, for k=3,4,5 (cos/sin). The parity block has one zero and one negative eigenvalue. Therefore exactly one soft and 17 stable tangent modes, provided the admitted strict spectral inequalities. These are *linearized local coefficients*, not an unconditioned generator convergence theorem. No fitted κ, new noise, or clock is introduced into primary FIN.

## HM-004: conditional stationary critical path theorem with a complete quantitative remainder contract

This is a **stationary reflected AUX-HB-OU** theorem only. Fix κ∈{1,4}, fixed bounded λ, `g_N=g_c+λ N^{-1/2}`, N divisible by 12. Start exactly from `Π_D,N(dn,dθ)=π_D,N(n)ρ_n(dθ)`, and reflect only the count proposals at the admitted local domain D_N. The existing MT-003 local-Gibbs result has the count marginal π_D,N, quartic variable `N^{1/4}m` converging to `ρ_λ(dz)∝exp[χλ z²/2−z⁴/12] dz`, and **polynomial moments of every fixed order**. Its transverse Gaussian moments give `E||w||^{2a}=O(N^{-a})`, with w the transverse occupation vector. The θ|n Gaussian formula adds `E||θ_⊥||²=O(N^{-1})`, `E|h|^q=O(N^{-q/4})`, `E|h−m|^q=O(N^{-q/2})` for fixed q (λ bounded). A sufficient explicit boundary hypothesis is `N^{3/4} Π_D(n in B_D)→0`, where B_D collects types with an allowed one-copy jump leaving D_N; the inherited fixed-barrier local Gibbs bound is exponentially smaller. **This boundary rate and MT uniform moment theorem must be audited against the local review before independent acceptance**. Nothing here assumes a native iid-start theorem for the *new* generator.

The parity softmax expansion is particularly economical. Decompose `Xθ=hσ+ξ` with ξ orthogonal to constants and σ. Because the parity-line derivative of φ in any ξ-direction is identically zero,

`|φ(θ)−tanh(h)| ≤ C||ξ||²` (uniform globally by bounded softmax derivatives),

and `|tanh(h)−h+h³/3|≤C|h|^5` (valid globally by bounded derivatives/Taylor plus polynomial domination). With `s=(κm+h)/(κ+1)` and `z_N(u)=N^{1/4}s(√N u)`, the *integrated* scaled drift discrepancy compared with `rκ[χλ z_N−z_N³/3]` is bounded in expectation on every fixed [0,S] by

`C_S N^{3/4}{ E||ξ||²+E|h|^5+E|h³−s³| + N^{-1/2}E|m−s| + Π_D(B_D) }`.

The terms are respectively `O(N^{-1/4}) + O(N^{-1/2}) + O(N^{-1/4}) + O(N^{-1/4}) + N^{3/4}Π_D(B_D)`. (For the h³−s³ bound use Hölder/Cauchy and stationary `h,s=O_L4(N^{-1/4})`, `h−s=O_L2(N^{-1/2})`, hence O(N^-1).) Boundary drift: at most N rate-one attempts times a 2/N parity increment, so the absolute reflected/full difference in Lm is ≤2 on B_D; its slow scaled contribution is the last term.

The scaled bracket is `N * d<s>/dt`: by stationarity and finite moments its integrated difference from `2rκ` tends to zero in L¹ on fixed horizons, with the reflected-count bracket subtraction ≤C Π_D(B_D). Jump sizes of `z_N` are ≤`2κN^{-3/4}/(κ+1)`; the h component has continuous Gaussian martingale noise. The initial scaled variable converges to ρ_λ since `N^{1/4}(s−m)→0` in L². Uniform stationary moments control drift/BDG increments and Aldous tightness. Any subsequential limit solves the well-posed quartic martingale problem with

`G_κ f(z)=rκ f''(z)+rκ(χλ z−z³/3)f'(z)`, i.e.

`dZ_s=rκ(χλ Z_s−Z_s³/3)ds+√(2rκ)dB_s`, `Z_0∼ρ_λ`.

The quartic confining drift has nonexplosive pathwise unique solutions; hence the entire sequence converges in D([0,S]) (continuous limit). Since κ is frozen, each result is pointwise; it does not cover κ_N, iid initialization for the diagnostic, unreflected model, or any physical time. It gives a *conditional analytical proof*, not an audited theorem independent of MT/OBS premises. If the boundary or polynomial-moment premise is absent from the actual inherited review, this HM-004 result must be downgraded to OPEN_BOUNDARY_OR_MOMENT_PREMISE; it may NOT be rescued with a simulation.

## HM-005: exact hidden-clock nonidentifiability (finite N)

Let `q_D(n,n')` be the ORIGINAL reflected native count rates, `π_D` their reversible count law and `ρ_n(dθ)` the static MP7 Gaussian posterior. For every fixed η≥0 define a **separate exact analysis lift** with transitions

`(n,θ) -> (n', θ')` at rate `q_D(n,n')ρ_{n'}(dθ')` for `n'≠n`, and additionally `θ->θ'∼ρ_n` at hidden rate η with n unchanged.

Joint invariant Π_D=π_D ρ_n. Detailed balance holds edge by edge:

`π_D(n)ρ_n(dθ) q_D(n,n')ρ_n'(dθ') = π_D(n')ρ_n'(dθ') q_D(n',n)ρ_n(dθ)`.

The count generator is exactly `q_D`, completely θ-independent; therefore **all native count path distributions are identical** for every η, given the same initial count distribution. Yet the hidden centred coordinate `v_ℓ(n,θ)=θ_ℓ−(g/N)(X^T n)_ℓ` (conditional variance g/N) has a different clock. On the event of *no count change during [0,t]*, with `d(n)=∑_{n'≠n}q_D(n,n')`,

`E[v_ℓ(t)v_ℓ(0) 1_{no count jump}] = (g/N) e^{-η t} E_{π_D} e^{-d(n)t}`.

This exact bounded-time covariance witness distinguishes hidden rates η while leaving every native count observation unchanged. The field itself is unbounded Gaussian; the product is integrable by its finite second moment. If a bounded witness is preferred, use bounded truncations and dominated convergence. Hence the static auxiliary identity and the full native count path **cannot identify a unique hidden mediator clock or mechanism**. Comparing the three objects gives: (i) ephemeral cavity latent draw EXACT native update; (ii) reset lifts EXACT native count path with many hidden clocks; (iii) persistent feedback AUX-HB-OU-κ SHARED Gibbs statics but NON-Markov projected count path and κ-dependent stationary critical kinetics. None is a discovered physical medium.
