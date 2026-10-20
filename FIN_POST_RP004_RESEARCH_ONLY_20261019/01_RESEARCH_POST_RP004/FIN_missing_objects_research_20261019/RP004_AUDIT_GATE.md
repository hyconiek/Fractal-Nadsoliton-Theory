# V-RP4-001 — limited independent checks; full proof audit not completed

Date: 2026-10-19. Status: **EXACT_MISSING_PROOF_INPUT** (not a counterexample and not ADMIT).

## Available inputs

The supplied current message provides a summary, a stated positive-flux constant, the task's audit checklist, and the earlier frozen shared contract in the user-provided ZIP. The actually claimed new proof files
`fin_rp004_flux_review_20261019/MOVING_CONE_THEOREM.md`, `BOOTSTRAP_COMPLETION.md`, `OUTWARD_FLUX_THEOREM.md`, and `COMPLETE_FLUX_BUDGET.json` are **not present in the locally mounted files or in the uploaded historical strategy archive**. The public GitHub raw path was also not retrievable (404 in this session). Without their proof lines, no honest line-by-line independent validation or exact-hash sign-off is possible.

## Independent controls from the frozen original contract

The original reduced equations and prepared data are

```
phi_tt - phi_xx + 2phi(phi²-1) + 6phi psi² = 0,
psi_tt - psi_xx + (a+6phi²)psi + g4 psi³ + g6 psi^5 = 0,
K=tanh x, psi_initial=eps sqrt(3/4)sech² x,
a in [-3/4,-1/4], phi_initial=K, both velocities zero.
```

From `phi=K+eps²u+...`, `psi=eps f cos(omega t)+...` the O(eps²) scalar forcing is **exactly**

`u_tt + Lphi u= -S(1+cos Omega t),  S=9/4 tanh x sech^4 x`, because `6K f² cos²(omega t)=3K f² (1+cos(2omega t))`.

`omega²=a+2`; `Omega=2omega`; `q=sqrt(Omega²-4)=2sqrt(a+1)`; the resonant group speed is `v_g=q/Omega` and, throughout the frozen a interval,

`1/sqrt(5) <= v_g <= sqrt(3/7)`, i.e. `[0.4472,0.6547]`.

The claimed asymptotic spatial band `x/t in [1/6,1/3]` is separated from `v_g` by at least `1/sqrt(5)-1/3>0`. For a **fixed** observation radius `R=T/6`, the late interval `t in [T/2,T]` lies in exactly this band. One must not substitute a moving boundary `R(t)` in the flux identity.

For all smooth solutions of the original conservative coupled equations, with properly vacuum-subtracted energy density `e` and current `j=-phi_t phi_x-psi_t psi_x`, `partial_t e+partial_x j=0`; therefore

`F_exact(R,T)=integral_0^T [j(t,R)-j(t,-R)]dt=E_{[-R,R]}(0)-E_{[-R,R]}(T)`.

By the frozen parity, phi and phi_t are odd; phi_x even; psi and psi_t even; psi_x odd. Thus `j(t,-R)=-j(t,R)` and `F_exact(R,T)=2 integral_0^T j(t,R)dt`. This is an exact check of sign conventions.

The full result reported by the user, `F_exact >= (625/12386304) eps^4 T`, uses a positive constant `~5.045895854e-5`. The arithmetic sign is consistent with the amplitude and group speed, but **these checks do not certify the hard stationary-phase and nonlinear bootstrap proofs**.

## Exact remaining validation target

For a nonzero outgoing pole amplitude `A_+(a)` derived from the sourced Jost transform, prove uniformly on the entire frozen a interval and on the band `1/6<=x/t<=1/3` that both derivatives of the **zero-data retarded** solution differ from the outgoing pole reference by `O(t^-1/2)` (or another error integrable to o(T)). Prove the necessary constants and separate the DC, shape, threshold, transient and cutoff terms. Independently establish the claimed **corrected triangular bootstrap**, remainder trace control, and final signed flux lower bound with all quantifiers and epsilon_*.

A conditional algebraic gluing argument is straightforward *if* the uniform cone theorem is paid: for `R=T/6`, the late outgoing right pole `u_p=A cos(qR-Omega t+theta)` satisfies `j_p=Omega q A² sin²(...)`; hence its time integral over `[T/2,T]` has positive main term `Omega q A² T/4+O(1)`. A trace error O(t^-1/2) contributes only O(sqrt(T)). The early interval `[0,T/2]` must still be treated with a positive exterior energy plus the exponentially small exterior source-work term, **not discarded**. This is a conditional bridge lemma, not an independent proof of the missing cone bound.

## Honest disposition

Do not mark this independent audit `ADMIT_CONDITIONAL_FINITE_WINDOW_THEOREM` on the supplied excerpt. Preserve the newly claimed theorem as `REPORTED_CONDITIONAL_PROOF_PENDING_INDEPENDENT_FULL_TEXT_AUDIT`. This limits verification, not the historical author's claimed result. No change to F01/F03/F06.
