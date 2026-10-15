# DKO-003 — Exactly stated clock and amplitude invariance

At fixed α, c=(1−α)/α>0 and h>0, let a common change of model-time unit multiply all relaxation rates by ν>0. With q=exp(−νh), the null Fourier mode γ_1=1 satisfies ρ_1(h)=q. For γ∈(0,1):

ρ_γ(h)=[q^γ+cγq]/[1+cγ]=q+(q^γ−q)/(1+cγ).

The measured ideal sign correlations determine ρ_r=sin(πs_r/2), s_r=E(sign(V_r(0))sign(V_r(h))). Thus obtain q from mode 1, then infer γ from the displayed function. For 0<q<1, 0<γ<1,

∂ρ/∂γ={q^γ log(q)(1+cγ)−c(q^γ−q)}/(1+cγ)^2 <0.

Hence γ is uniquely identified *in the prescribed ideal family* without prior knowledge of the common ν. A fixed nonzero multiplicative gain a_r in each channel (even negative) multiplies both time samples and does not change their sign product. It does **not** eliminate additive offsets, errors in sign, unknown channel-specific clocks, nonlinear distortion, a changed block-sampling convention, unknown label-to-mode mapping or drift/measurement contamination. All formulas concern ensemble probabilities, not a single observed product.

When γ_1=1, null channel ρ_1=e^(-νh) is a legitimate internal clock reference only because the null eigenvalue Λ1=0 in A7 and both operators, under the same microscopic update-clock convention. Rescaling a common clock does not make ν into physical seconds.
