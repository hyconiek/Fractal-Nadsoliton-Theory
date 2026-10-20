# ID-001 — quantitative stationary-joining entropy theorem

**Status: analytical conditional theorem; scope limited to the declared source class.** Date: 2026-10-19. All logarithms are natural. This is not a new entropy inequality in general: Fano, data processing, and the entropy power rule are inherited. Their *joint specialization* to the precise N2 reset generator and refinement is the new packet.

## Freeze

Let `(Ω, F, μ, (T_t)_{t∈R})` be a measure-preserving measurable flow on a standard probability space. Define its time-one Kolmogorov–Sinai entropy `H=h_μ(T_1)<∞`. `H` is per UNIT OF THE SUPPLIED MODEL TIME; no SI second or physical clock is derived. At a sampling interval `δ=1/m`, where `m∈N`, choose *any* deterministic measurable map `f_δ: Ω→{0,...,11}` and set `Y_k=f_δ(T_{kδ}ω)` for k∈Z. `f_δ` may vary with δ; the underlying invariant source measure/flow and finite H stay fixed.

Let X be the stationary q12 FIN N2 difference CTMC sampled at intervals δ, defined precisely in NATIVE_RESOLUTION_BOUND.md. Assume a stationary **joining of whole bi-infinite sampled processes** `(X_k,Y_k)_{k∈Z}` with the prescribed marginals and shift-invariant joint law. Let `e_δ=P(X_0≠Y_0)`; permit a new joining for each m. Inversion statements below use `0≤e_δ≤11/12`. No identification of *physical time* or extra random readout is licensed.

Define entropy rate for stationary finite-alphabet sequence Z by `h(Z)=lim_{n→∞}H(Z_0,...,Z_{n-1})/n`, with subadditive existence. Set `h_δ=h(X)`.

## Theorem A: all-resolution source/error budget

Under the above assumptions, for every positive integer m,

`h_δ ≤ δ H + h_2(e_δ) + e_δ log(11)`,   `δ=1/m`, 

where `h_2(e)=-e log e-(1-e) log(1-e)` and `0 log0=0`. It is valid even if `e_δ>11/12`; the restriction is used only to invert the Fano term monotonically.

### Proof by conditional block entropy

For each index k let `E_k=1{X_k≠Y_k}`. Conditional on `Y_k`, identifying `X_k` requires at most the error flag E_k (binary entropy) and, when wrong, one of at most q−1 other symbols. Thus, without any independent-error assumption,

`H(X_k | Y_k) ≤ H(E_k) + e_δ log(q−1) = h_2(e_δ)+e_δ log11`.

The stationary joining gives `P(E_k=1)=e_δ` for every k. For n-blocks,

`H(X_0^{n-1}) ≤ H(Y_0^{n-1}) + H(X_0^{n-1}|Y_0^{n-1})`
`≤ H(Y_0^{n-1}) + Σ_{k=0}^{n-1} H(X_k | Y_k)`
`≤ H(Y_0^{n-1}) + n[h_2(e_δ)+e_δ log11]`.

Divide by n and take n→∞. Since `Y_k` is a finite measurable partition observation of `T_δ`, its entropy rate is at most `h_μ(T_δ)`. The integer power identity `h_μ((T_δ)^m)=m h_μ(T_δ)` and `(T_δ)^m=T_1` imply `h_μ(T_δ)=H/m=δH`. Combine. QED.

### Independent proof check (entropy via flags)

Let `E^n` denote the whole error flag string and note that `E^n` is determined by `(X^n,Y^n)`. Chain rule gives `H(X^n|Y^n)=H(E^n|Y^n)+H(X^n|E^n,Y^n)`. Bound first by `H(E^n)≤Σ H(E_k)`, second by `Σ P(E_k=1)log11=n e_δ log11`; divide by n. This separately checks the handling of dependencies and the 11 rather than 12 factor.

## Theorem B: exact finite-δ inverse budget

Define `F_12(e)=h_2(e)+e log11`. For `0≤e≤11/12`, `F_12'(e)=log[11(1−e)/e]≥0`, and `F_12` maps `[0,11/12]` continuously onto `[0,log12]`. Hence every admitted joining has

`e_δ ≥ F_12^{-1}(max{0,h_δ−δH})`.

This is a **necessary**, not sufficient condition; the inverse applies only when the assumed e lies in the monotonic branch. It gives an exact numerical lower bound once the frozen p and H are specified. An arbitrary e greater than 11/12 already represents poor reconstruction and needs no inverse claim.

## Theorem C: no sublinear-δ mismatch with uniformly finite entropy

Let `p=(p_0,...,p_11)` be any fixed strictly positive N2 distribution and `r̄=2(1−Σp_i²)>0`. From the exact native entropy expansion in NATIVE_RESOLUTION_BOUND.md, `h_δ/δ = r̄ log(1/δ) + c(p) + O_p(δ |logδ|)`. For every fixed `H<∞` and any allowed sequence of stationary joinings,

`liminf_{m→∞} e_{1/m}/(1/m) ≥ r̄`.

Proof. Suppose to the contrary that for some `η>0`, an unbounded subsequence satisfies `e_δ ≤ (r̄−η)δ`. For small δ this upper bound lies below `11/12`; `F_12` is increasing on the interval. Since `F_12(aδ)=aδ log(1/δ)+O_a(δ)` for each fixed a>0 (and equals zero if a=0), Theorem A yields `h_δ ≤ δ H + (r̄−η)δ log(1/δ)+O(δ)`, contradicting `h_δ=r̄ δ log(1/δ)+O_p(δ)`. QED.

This implication is deliberately **not** the claim `liminf e_δ≥r̄`: actual `e_δ` may tend to zero, yet cannot be `o(δ)` under the specified fixed finite-H budget. No lower bound for ordinary one-window weak approximation follows without the joining hypothesis.

## Theorem D: required entropy if H is allowed to grow with resolution

If each δ instead has its own deterministic source with entropy `H_δ` per common time unit, then

`H_δ ≥ [h_δ−F_12(e_δ)]/δ`.

Exact matching `e_δ=0` requires `H_δ≥r̄ log(1/δ)+c(p)+O_p(δ|logδ|)`; if `e_δ=o(δ)`, the necessary leading-order growth remains `H_δ ≥ [r̄−o(1)]log(1/δ)`. The latter follows because `F_12(e_δ)=o(δ log(1/δ))` for e=o(δ). The assertion is a resource bound, not a construction of suitable sources.

## What this proof does not supply

It does not determine the intrinsic FIN generator, explain why heat-bath kinetics should be physical, or derive a relational or calibrated physical time. It does not prohibit stochastic, singular, infinite-entropy, or infinite-dimensional dynamics. It does not imply a laboratory observation or term-by-term entropy production. It does not turn positive KS entropy into thermodynamic irreversibility. It does not rely on the QW-206 histogram.
