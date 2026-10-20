# Original-current audit and one epsilon threshold

Date: 2026-10-20. Uses the same T,R, action and Cauchy data as the source.

## 1. Reference energy and radius convention

The outgoing harmonic field has average density Omega^2|C|^2/2, not an
assumed energy power at the boundary. On BOTH bands of total length T/3,
the oscillatory spatial term integrates to O(1). The paid pointwise C1
error changes energy by O(sqrt(T)+1). Since |C|^2=Gamma/(2q), the coefficient
is Omega^2 Gamma/(12q). The inherited Gamma floor and q<7/4 give
c0=625/1548288. Outside R>=1 the potential 4-6sech^2 x is positive.

For the forced reference, direct exterior integration gives

    F_u(R,T)=E_ext(T;R)+integral_0^T integral_|x|>R
                           S(1+cos(Omega t))u_t dx dt.

The work term is NOT assumed positive. Its absolute bound is
C T^2 exp(-4R). For each final T the radius R=T/6 is held fixed from
time zero. A truly moving R(t) would add integral R'(t)[e(t,R)+e(t,-R)]dt
to the flux identity. That is a different experiment and is not used.

## 2. Every term of the actual nonlinear current is retained

Expand J=-phi_t phi_x-psi_t psi_x using
phi=K+epsilon^2u+w and psi=epsilon v+r. The source inventory (17) is
complete. At either endpoint use

    |K'|<=4exp(-2R), |v_t|<=6exp(-2R), |v_x|<=8exp(-2R),
    |u_t|,|u_x|<=420tau,
    |r_t|,|r_x|<=Cr epsilon^3tau^2,
    |w_t|,|w_x|<=Cw epsilon^4tau^3.

The trace inequalities in the preceding proof justify these values. After
integrating BOTH endpoint absolute errors, sufficient coefficients for
the five error groups of the submitted budget are respectively

    D1=1680, D2=96, D3=2Cw+10Cr,
    D4=(1680Cw+2Cr^2)/5, D5=2Cw^2/7.

Thus a safe common constant is

    C_J=1+D1+D2+2Cw+10Cr+1680Cw+2Cr^2+2Cw^2.

On the original logarithmic window the complete current error is at most

    C_J[epsilon^6 tau_T^5+epsilon^8 tau_T^7+epsilon^10 tau_T]. (1)

No kink interference, bound-mode tail, r-r or w-w term is post-subtracted.
The exact parent has e_t+J_x=0 for the supplied potential. The forced work
term belongs only to u, not to the parent. Initial excitation is order
epsilon^2, so an order epsilon^4 log(1/epsilon) outward transfer does not
contradict conservation or require a perpetual external drive.

## 3. Uniform quantifiers, not only a sequence of epsilon values

The cone proof supplies a finite UNIFORM C_ref>=1 such that

    F_u(T/6,T)>=c0T-C_ref sqrt(T)-C_ref-C_ref T^2exp(-2T/3).

For T>=6, T<=exp(T/3). A sufficient reference threshold is

    T0=max[6,(4C_ref/c0)^2,4C_ref/c0,3log_+(4C_ref/c0)]. (2)

Then each of the three losses is at most c0T/4, giving F_u>=c0T/4.
No numeric value of C_ref is asserted by finite quadrature.

Let L=log(1/epsilon), T=12L. Choose an integer L_*>=max(3,T0/12)
large enough that

    exp(-2L_*)(1+12L_*)^3<=K0,
    C_J[exp(-2L_*)(1+12L_*)^5
        +exp(-4L_*)(1+12L_*)^7
        +exp(-6L_*)(1+12L_*)]/(12L_*)<=c0/8.         (3)

Such an integer exists. Every left side is decreasing for L>=3: logarithmic
derivatives are bounded by -2+3/L for the bootstrap, and -2+4/L,
-4+6/L, -6 respectively for the error-ratio terms. Consequently
epsilon_*=exp(-L_*) pays all requirements for EVERY 0<epsilon<=epsilon_*,
not merely a selected subsequence. Equations (1)--(3) give

    F_exact(R,T)>=(625/12386304)epsilon^4T>0.

Constants depend only on the declared finite coefficient bounds and fixed
spectral/Sobolev/source data, uniformly in a. This is a sufficient threshold
recipe with an uncalibrated finite C_ref, not useful physical timing or
an all-time scattering theorem.
