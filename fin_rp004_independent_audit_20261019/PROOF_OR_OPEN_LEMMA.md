# First unpaid analytical estimate — `EXACT_UNPAID_ESTIMATE`

**Attempted conclusion:** `F_exact(2 log(1/eps), 12 log(1/eps)) ≥ (625/12386304)eps^4 12 log(1/eps)` for all `0<eps≤eps_*` and every frozen `a∈[-3/4,-1/4]`, with declared finite fixed coefficients `g4,g6`.

**First precise unverified prerequisite:** The actual new `MOVING_CONE_THEOREM.md` is unavailable here. To independently prove it (rather than assume it), fix `I=[-3/4,-1/4]`, `u(0)=u_t(0)=0`, `L=-∂_x²+4-6sech² x`, `S=(9/4)tanh x sech^4 x`, `Omega=2 sqrt(a+2)`, `q=2sqrt(a+1)`. Establish a finite `C,T0`, **independent of `a`**, and the true outgoing pole field `u_out` with sourced Jost residue and exact signs, such that

    sup_(a∈I) sup_(T≥T0) sup_(T/6≤|x|≤T/3)
      sqrt(T)*( |u(T,x)-u_out(T,x)|
               +|∂_t(u-u_out)(T,x)|
               +|∂_x(u-u_out)(T,x)| ) ≤ C.

Require an appropriate formulation of `u_out` also covering the real DC bound and shape mode, either in the reference or controlled remainder. The distributional `i0` resonance must include both ±q poles and a rigorously controlled nonresonant smooth cutoff, high frequency and threshold, no incoming field. A weaker integrated band estimate sufficient for `E_band(T)≥cT` with a uniform c>0 would also pay the logical requirement but not automatically the specific new theorem's stronger claim.

**Next unpaid prerequisite:** The new `BOOTSTRAP_COMPLETION.md`, required to repair the previous wrong constant closure and produce finite constants `C_r,C_w,eps_*` depending only on the fixed family, is also unavailable. Verify precise nonlinear norm and trace estimates through `T=12log(1/eps)` and current error `o(eps^4 T)` uniformly. The new `OUTWARD_FLUX_THEOREM.md` and `COMPLETE_FLUX_BUDGET.json` are required to check the claimed prefactor exactly. Do not accept the pre-repair RP004-B as a substitute.

**Verdict:** `EXACT_UNPAID_ESTIMATE` because the decisive latest proof *materials* are absent in this session. Not a counterexample, not a mathematical claim of impossibility, and not `ADMIT_CONDITIONAL_FINITE_WINDOW_THEOREM`.
