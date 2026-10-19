# RP-004 — exact open nonlinear outward-energy estimate

**Disposition: OPEN_FLUX_DOMINATION_FOR_NEUTRAL_CAUCHY_TRAJECTORY.** RP-002 and RP-003 are certified conditionally, but no complete controlled finite-interval positive outward-flux theorem has been obtained. Do not relabel Gamma>0 as measured outward energy.

## What is paid

Canonical total local flux J=-phi_t phi_x - Re sum_j (conj(Psi_j,t) Psi_j,x), with the exact energy identity

    integral_0^T [J(t,R)-J(t,-R)] dt=E_[-R,R](0)-E_[-R,R](T),

where E is the complete original canonical energy density (vacuum constants handled consistently). The same full coupled PDE and the NS small-data orbital bound apply to the frozen F_RP family. No separately appended radiation equation or packet is introduced.

For each FIXED finite smooth time interval [0,T], smooth local wellposedness, analytic polynomial dependence, and invariance under Psi -> -Psi yield the controlled formal orders (with constants depending on T):

    Psi_epsilon=epsilon f cos(omega t)+O_T(epsilon^3) in H1 x L2,
    phi_epsilon=K+epsilon^2 u+O_T(epsilon^4) in H1 x L2,

where u is the retarded solution in the phase file. Standard Duhamel subtraction and energy estimates on the linearized self-adjoint operators yield these orders for smooth small data over each fixed T, after accounting for the translation zero mode (here the quadratic source is odd). They do NOT establish a uniform estimate for T or R increasing as epsilon ->0.

## Unpaid dominance when extracting transported radiation

At each boundary x=+/-R the total current includes

    J= -epsilon^2 u_t K' -epsilon^2 (f cos omega t)_t (f cos omega t)_x
       -epsilon^4 u_t u_x + [remaining O_T(epsilon^4) and background-error terms].

The first two are *order epsilon^2 at any fixed R*. The first represents interference with the kink tail; the second is oscillatory near-field flux of the exponentially localized transverse bound mode. Their integrated signs need not equal radiative energy, and a positive J at one early time is not proof of continuum transport. The scalar u includes both discrete and continuum pieces and an unavoidable turn-on/DC response. The radiative quadratic scalar energy flux is order epsilon^4: it cannot dominate those fixed-R terms simply by invoking Gamma>0.

A plausible isolation window must obey at least R_epsilon ~ log(1/epsilon) to make K'(R)~exp(-2R) suppress the O(epsilon^2) cross current beneath O(epsilon^4). Even with that spatial choice, causality requires T_epsilon >= R_epsilon for arrival from the core, while massive continuum group velocity q/Omega <1 makes the on-shell response later. Hence this requires a **long-time** bound (at minimum logarithmic in 1/epsilon), a precise turn-on/wave-front asymptotic, discrete-mode and bound-mode tail payments, and a strictly positive integrated flux margin. None of these is implied by a fixed-T Duhamel bound or NS orbital stability. The present package makes no numerical transient substitution or unproved stationary-flux assumption.

The missing statement is an estimate of the form

    F_exact(R_epsilon,T_epsilon) >= epsilon^4 C_rad(R_epsilon,T_epsilon)
                                    - E_background - E_bound - E_slip
                                    - E_nonlinear >0,

where every term is evaluated for one prescribed neutral Cauchy solution and where C_rad>0 is a **finite-time retarded outgoing contribution**, not the eternal Gamma alone. In particular a nonnegative global energy or a spectral density is not a substitute for the signed canonical flux. One must also avoid a trivially positive tail-induced flux being marketed as radiation.

**Stop:** ON_SHELL_COUPLING_ONLY; CERTIFIED_CONDITIONAL_OUTWARD_TRANSPORT is not established. No asymptotic scattering, second-wall absorption, electromagnetic radiation or photon inference.
