# Conditional moving-cone-to-original-flux gluing lemma (independent derivation)

**Assumptions (not discharged here):** Let the frozen potential, parity, retarded data and coefficient interval `a∈[-3/4,-1/4]` be unchanged. For the zero-data solution `u_tt+L u=-S(1+cos Ωt)` suppose, uniformly in `a` and in the two outgoing spatial bands `T/6≤|x|≤T/3` at observation time `t=T`, that

    sup_(band) (|u-u_out|+|∂_t(u-u_out)|+|∂_x(u-u_out)|) = O(T^(-1/2)),

where `u_out` has outgoing harmonic amplitudes of modulus `A=|M(q)|/(2q)` on both half-lines; any exponentially localized DC/shape part must be INCLUDED in the remainder. Assume the spectral certificate `Gamma=|M(q)|²/(2q)>125/73728` from the frozen parent. No incoming wave or new preparation may be inserted.

**Claim 1 (reference field, conditional on cone estimate):** On the two bands at time T,

    E_band(T) = [Ω² Γ/(12q)] T + O(T^(1/2))+O(1),

uniformly on the fixed compact parameter interval. Derivation: on either band the leading field is `A cos(q|x|-ΩT+θ)`, up to parity sign. Since Ω²=q²+4, the average of `e=(u_t²+u_x²+4u²)/2` is `Ω² A²/2`. The total band length is `T/3`; the oscillatory spatial correction is bounded O(1), and the uniform O(T^-1/2) field/trace remainder integrates to O(T^1/2). The potential `-6sech² x` contributes exponentially small tails. `A²=Gamma/(2q)`.

**Claim 2 (fixed-radius outgoing flux):** For `R=T/6`, `E_ext(T;R)≥E_band(T)` since `4-6sech² x>0` for `|x|≥1`. The exact exterior energy identity for the *forced* reference yields

    F_u(R,T) = E_ext(T;R) + ∫_0^T∫_(|x|>R) S(x)(1+cos Ωt) u_t dx dt.

The second term is bounded below by `-C T² e^{-4R} = -C T² e^{-2T/3}` from the source tail and linear graph energy; its sign must not be discarded. Hence

    liminf_(T→∞) inf_(a∈I) F_u(T/6,T)/T
       ≥ inf_(a∈I) Ω² Gamma/(12q)
       > 5*(125/73728)/(12*sqrt(3))
       =625/(884736 sqrt(3)) ≈ 4.07854906145e-4.

**Claim 3 (nonlinear coupling conditional on bootstrap):** The original pre-repair RP004-B states at `T=12 log(1/eps), R=2 log(1/eps)`,

    F_exact(R,T)=eps^4 F_u(R,T)+o(eps^4 T).

The statement must be **re-established** using the corrected triangular bootstrap, preserving sign and finite family-bound dependence. If that input is fully proved, Claim 2 gives, for all sufficiently small eps,

    F_exact(R,T) ≥ c eps^4 T,  for any fixed c<625/(884736 sqrt(3))

with sufficiently small eps depending on c, M4, M6 and fixed gap constants. In particular, the user-reported `c=625/12386304` is strictly smaller than the conditional limit constant. This does **not** independently verify that the author's claimed `epsilon_*` is legitimate.

**Caution on spatial/time band:** At `t=T`, the whole `[R,T/3]` lies behind the resonant group velocity. For a *fixed* R=T/6 used in `F_u` the late times `t∈[T/2,T]` have `R/t∈[1/6,1/3]`; the early contribution must be covered by the exterior energy identity. `R(t)=t/6` would have a moving-wall term and is NOT interchangeable with fixed `R`.

**Conclusion:** Full finite-window flux is rigorously reduced to TWO as-yet unaudited proof obligations, rather than proved from their declarations.
