# A controlled Volterra equation for the actual bound-mode coordinate

Date: 2026-10-21. Conditional on the actual source model, invariant class,
Darboux resolution and accepted strong remainder proof. Same-author work.

## 1. Source equations and the inherited strong window

Use the contract's affine eta and real psi. Exact expansion gives

    eta_tt+H eta=-6K psi^2-6eta psi^2-6K eta^2-2eta^3,
    psi_tt+L_a psi=-12K eta psi-6eta^2 psi-g4 psi^3-g6 psi^5. (1)

The accepted RP bootstrap supplies constants C,K0, uniform in a0 and
the declared finite M4,M6, such that on epsilon^2(1+T)^3<=K0

    ||psi-epsilon f cos(omega t)||_(H2 x H1)<=C epsilon^3(1+t)^2,
    ||eta||_(H2 x H1)<=C epsilon^2(1+t).                    (2)

These are reused actual proof bounds, not inferred from a numerical flag.
On odd eta, H has gap3; on the prepared even psi sector L_a has gap>=5/4.
Their spectral strong Duhamel estimate maps L1_t H1 forcing to H2 x H1
with a uniform constant. The background K itself is not called an H2 state.

Let z=<f,psi> and w=psi-zf. Projection of (2), including velocities, gives

    |z|+|z_t|<=C epsilon, ||w||_(H2 x H1)<=C epsilon^3(1+t)^2. (3)

Shrink an ANALYSIS constant kappa<=K0 if needed; no parameter, preparation
or generator is changed. All product estimates below use the1D Sobolev
algebra and the smooth fixed profiles K,f.

## 2. Eliminate eta with its ACTUAL causal data

Put F=K f^2 and

    eta_z(t)=-6 integral_0^t sin((t-s)sqrt(H))/sqrt(H) F z(s)^2 ds.

This is defined on the odd sector; no inverse of the translation zero mode
is inserted. Both eta_z initial data are zero. The difference xi=eta-eta_z
is driven by

    -6K(psi^2-z^2 f^2)-6eta psi^2-6K eta^2-2eta^3.

By (2)–(3), this H1 forcing is O(epsilon^4(1+t)^2). For example
psi^2-z^2 f^2=2zf w+w^2, and the higher term is absorbed using
epsilon^2(1+t)^2<=kappa. Hence

    ||xi||_(H2 x H1)<=C epsilon^4(1+t)^3.                 (4)

Define the scalar causal kernel

    k(t)=<F,sin(t sqrt(H))/sqrt(H)F>, t>=0.                (5)

Projection of (1), using L_a f=omega^2 f and (4), now gives

    z_tt+omega^2 z+(18/35)g4 z^3
       -72 z(t) integral_0^t k(t-s)z(s)^2 ds=R_epsilon(t),
    |R_epsilon(t)|<=C epsilon^5(1+t)^3.                   (6)

The cubic self coefficient is integral f^4=18/35. The retarded coefficient
is (-12)*(-6)=+72 on the RHS, hence -72 on the displayed LHS. Its sign
and source S=3F are not adjustable. Terms eta*w and z*xi have the maximal
epsilon^5(1+t)^3 size; eta^2 psi, psi^3-z^3 f^3 and psi^5 are also paid.
No signed g6 is set to zero. In fact at t=0 the residual is exactly
R_epsilon(0)=-(24/77)g6 epsilon^5, since integral f^6=24/77.

## 3. A genuine scalar prediction, not only a projected identity

Let zbar solve (6) with R=0 and zbar(0)=epsilon, zbar_t(0)=0. The bounded
continuous kernel and polynomial history map give local existence and
uniqueness by a Volterra/Picard argument. Since H>=3 on the actual odd
profile, |k(t)|<=||F||^2/sqrt(3)<=2/35 for all t.

Fix the scalar amplitude bootstrap constant before choosing kappa: take it
larger than twice the uniform free-oscillator bound. On a bound
|zbar|+|zbar_t|<=L epsilon, its nonlinear force is bounded by
C_L epsilon^3(1+t). Oscillator Duhamel changes the free amplitude by at
most C_L epsilon^3(1+t)^2. Now choose kappa small enough to improve that
bootstrap using epsilon^2(1+T)^2<=kappa. Continuation gives zbar on the
whole stated window. This constant ordering is not a circular bootstrap.

For D(t)=sup_(s<=t)(|z-zbar|+|z_t-zbar_t|), cubic/history Lipschitz bounds
give

    D(t)<=C epsilon^2 integral_0^t (1+s)D(s)ds
                          +C epsilon^5(1+t)^4.

Gronwall and epsilon^2(1+T)^2<=kappa therefore prove

    |z-zbar|+|z_t-zbar_t|<=C epsilon^5(1+t)^4.             (7)

Constants depend on the fixed source bounds and chosen analysis kappa,
not epsilon,t or a0 in its frozen interval. This supplies a nonlinear
causal reduced predictor derived from the same full trajectory.

## 4. Scope, not an asymptotic damping promotion

The window includes the original logarithmic RP flux interval and the
stated source smallness regime. It does NOT reach the natural epsilon^-2
amplitude-damping scale. No cubic averaged envelope, all-time relaxation,
white noise, Born statistics, particle law or physical time follows from
(6)–(7). The derived memory has the actual shape and continuum components
in [SPECTRAL_MEMORY_KERNEL.md](SPECTRAL_MEMORY_KERNEL.md); neither can be
dropped because a convenient instantaneous friction law is desired.
