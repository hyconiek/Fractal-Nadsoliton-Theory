# KO-003 — One frozen bounded kernel contrast

**Status:** strict, ideal, asymptotic model contrast under the two KO-001 local-window/KO-002 snapshot hypotheses. Exact monotonicity is analytical. The decimal endpoints in `INDEPENDENT_VERIFICATION.json` use mpmath.iv interval evaluation of the supplied kernel, not independent empirical measurements.

Let B_N have exactly N/2 members, chosen before seeing any labels. The readout is one q12 block snapshot at the common post-MT origin. Freeze the model-independent, no-tuning event

`E_N = { |V_{3,N}| > 1 }`, `V_{3,N}=sqrt(N)(2/N)Σ_{a∈B_N} sqrt(2)cos(2π3 ξ_a/12)`.

By KO-002, under ε∈{0,.01}, `V_{3,N}⇒G_ε` with `G_ε∼N(0,v_ε)`, `v_ε=1+1/γ_3(ε)`, `γ_3(ε)=1−(Λ3+ε)/Λ6`.

The boundary ±1 has Gaussian probability zero, so the portmanteau theorem applies directly to the bounded event (no unbounded moments need TV transfer):

`lim_N Pr_ε(E_N) = p_ε = 2 Φ(−1/sqrt(v_ε)) = erfc(1/sqrt(2v_ε))`.

Differentiating for positive γ, `∂v/∂ε=1/(Λ6 γ^2)>0`. Since `p(v)=2[1−Φ(v^{-1/2})]`, `p'(v)=φ(v^{-1/2}) v^{-3/2}>0`. Hence **p_.01>p_0 strictly**. The distinction is neither a new fitted threshold nor a Monte Carlo fluctuation.

From the supplied exact-decimal kernel with inherited A7 construction, the standalone interval calculation reports nominal centres:

| | ε=0 | ε=1/100 |
|---|---:|---:|
| Λ3 | 1.9614068619764455 | 1.9714068619764455 |
| γ3 | 0.16257283698729818 | 0.15830331402779933 |
| Var(V3) | 7.15108906586978 | 7.31698714674029 |
| p(E) | 0.70844127966426925 | 0.71161592530141653 |

The limiting event-law total-variation gap between Bernoulli variables is

`Δ=|p_.01−p_0|≈0.0031746456371472875`.

Under equal model priors and one ideal two-outcome observation, its optimal binary model-classification error is exactly

`R*=(1−Δ)/2≈0.49841267718142636`.

This is a very small discrimination advantage, 0.0015873 below random guessing. The event distinguishes the two particular *specified mathematical operators*, not an unknown natural A7 source. An experiment with uncertain A, time/readout gain, finite N, sampling bias, instrument noise or unspecified preparation is not covered. No certified finite-N event error or sample complexity is provided. This result does not contradict the earlier RG-F no-go over broad physical uncertainty classes; those allow nuisance variation absent in this ideal comparison.

The observable is kernel-sensitive because λ3 changes while λ6, g_c, local parity process, and Gaussian sampling noise (set by k=N/2) are fixed. Model-specific: γ3,γ4 and variance v3,v4. Symmetry-protected: soft parity drift and bracket. Class-level: quartic universal critical law and fixed finite parity-record predictions.
