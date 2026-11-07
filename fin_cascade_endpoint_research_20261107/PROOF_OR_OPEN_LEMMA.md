# CE - full cubic endpoint control and the combined threshold functional

Date 2026-11-07. Parent c6111d1638e5020ffee9f7421703b952879beb2c.
Self-audited conditional analysis using the supplied RP model and the
recorded QP, CN, CS, FH and TB estimates.

## 1. Endpoint estimate for the quadratic multiplier

Use CN's cross-channel symbol b(p,q)=1/[Omega_psi(p+q) Phi(p,q)],
with all eight choices of output/input signs. The predecessor quadratic
phase theorem supplies uniform bounded derivatives of b on the mass band.
Only derivatives through order two are needed in this section.

Let a(p,q)=b(p,q)/[<p><q>]. Every derivative of a through order two is
in L2(R^2), uniformly in the signs and mass parameter. This follows from
the bounded derivatives of b and the integrability of <p>^-2<q>^-2.
Let K_a(y,z) be its inverse Fourier kernel. Plancherel and Cauchy-Schwarz
with the integrable weight (1+y^2+z^2)^-2 imply

    ||K_a||L1 <=C sum_(|beta|<=2) ||partial^beta a||L2.

For f,g in H1, put F=<D>f and G=<D>g. The bilinear operator has the
kernel representation

    B_b(f,g)(x)=integral K_a(y,z)F(x-y)G(x-z) dy dz,

with the Fourier normalization constant included in K_a. Applying
Cauchy-Schwarz in x and then the kernel L1 bound proves

    ||B_b(f,g)||L1<=C ||f||H1 ||g||H1.                (CE1)

The representation follows first for Schwartz data and then by density.
This estimate pays the L1 output without demanding a spatial moment of
either input or a derivative beyond the already available H1 norms.

## 2. The actual complete cascade in physical phase

Retain CN6 and CS1. For each output sign sigma set

    S_sigma=sum_(tau,rho) [
        B_(sigma;tau,rho)(F_G,phi^tau,h_psi^rho)
       +B_(sigma;tau,rho)(h_phi^tau,F_G,psi^rho)],
    R_C^sigma=-sigma a K S_sigma, a=-3.

The paid norms are ||h_i^sigma||H1<=C epsilon and
||F_G,i^sigma||H1<=C epsilon^2. Thus CE1 gives

    ||S_+||L1+||S_-||L1<=C epsilon^3.                 (CE2)

The real physical phase before projection is exactly

    A_C=-(a/2) K(S_+-S_-),
    B_C=-(a/(2i)) Omega_psi[K(S_++S_-)].              (CE3)

The momentum component involves a free multiplier of order one. We
estimate its low-frequency Darboux transform by duality. The H1/L2 phase
regularity of CE3 is already supplied by CN and CS.

## 3. Bounded distorted test functions pay the momentum conversion

For |k|<=kappa<=1 let e_k(x)=Wcal exp(ikx). The exact formula is

    e_k(x)=[3K^2-1-k^2-3ikK]exp(ikx)
                          /sqrt((1+k^2)(4+k^2)).

The functions e_k and their first two spatial derivatives are uniformly
bounded on this compact k interval. The heat representation

    Omega_psi^-1=(1/sqrt(pi)) integral_0^infinity
                         t^-1/2 exp(-m^2 t) exp(t partial_x^2) dt

has an L1 convolution kernel of norm 1/m. Consequently

    Omega_psi e_k=Omega_psi^-1(-e_k''+m^2 e_k)

is bounded uniformly in x and |k|<=kappa. This identity holds for the
bounded smooth functions here by the distributional Fourier calculus,
and the last expression is a convergent convolution with an L1 kernel.

Since Wcal-adjoint P_e=Wcal-adjoint, the correctly projected and converted
cascade source has free Fourier coordinates

    Fourier(Wcal-adjoint P_e A_C)(k)
                             =(2pi)^-1/2 integral conjugate(e_k) A_C,
    Fourier(Wcal-adjoint P_e B_C)(k)
       =-(a/(2i sqrt(2pi))) integral conjugate(Omega_psi e_k)
                                                    K(S_++S_-).

The second identity transfers Omega onto the test function. The L1
bound CE2 justifies both pairings. Density yields the identities for the
actual H1 input class. Therefore the free phase energy density satisfies

    e_C(s,k)<=C epsilon^3, |k|<=kappa.                (CE4)

All physical phase projections are retained through the exact identity
Wcal-adjoint P_e=Wcal-adjoint. The localized potential and projection
sources assigned in CN9 to R_2loc and R_3loc retain their existing bounds.
The cubic cascade here is precisely CS1's F_G contribution.

## 4. Full actual cubic trace and local-continuum improvement

TB8 bounds the direct component's low-output density by C epsilon^3.
Adding CE4 proves for the complete converted CS cubic forcing

    e_3(s,k)<=C epsilon^3, |k|<=kappa.                (CE5)

The source also has the previously paid global phase norm C epsilon^3.
Use TB1--TB5 in the low output and FH10 in the fixed high output to get

    integral_0^(S epsilon^-2) |q_3(t)|^2dt
                         <=C_S epsilon^2[1+log(1/epsilon)]. (CE6)

Here and below the cutoff kappa is fixed, and its dependence is absorbed
in the constant. Insert CE6 into CS4. This gives the improved actual
whole local-continuum estimate and its compact trace estimate

    integral_0^(S epsilon^-2) E_local(t)^2dt
                         <=C_S epsilon^2[1+log(1/epsilon)]. (CE7)

The source-history bounds remain valid with all interactions and signed
fixed coefficients. The logarithmic loss remains in this upper bound.
No temporal cancellation between different source times was used.

## 5. The complete threshold coefficient

At k=0 the distorted test function is r0. Define the two threshold
coordinates of the complete free cubic phase forcing by

    gamma_1=-(a/(2 sqrt(2pi))) integral r0 K(S_+-S_-),
    gamma_2=(1/sqrt(2pi)) [integral r0 G_D
           -(a/(2i)) integral (Omega_psi r0)K(S_++S_-)],
    G_D=-6v^2u-g4u^3.                                (CE8)

The threshold energy density is sqrt(m^2|gamma_1|^2+|gamma_2|^2).
The term Omega_psi r0 is bounded by the preceding argument, making both
formulas meaningful under the paid assumptions. These are instantaneous
functionals of the actual continuous phase, including its velocities.
The output half-wave branches are conjugate for real physical input,
so the two coordinates in CE8 are real.

The cascade in CE8 uses only G_phi=-6K(v^2+u^2), G_psi=-12Kvu and the
instantaneous phase. It has no explicit g4 dependence. Hence at a fixed
phase snapshot

    partial gamma_1/partial g4=0,
    partial gamma_2/partial g4
                         =-(2pi)^-1/2 integral r0 u^3. (CE9)

This derivative fixes the input snapshot. An actual trajectory also
changes with g4, so CE9 supplies no parameter derivative of that trajectory.
TB's projected snapshot already has a nonzero direct threshold pairing.
Thus a cancellation identity valid for every phase snapshot and every
g4 would require an additional restriction. The numerical comparison
below computes the full CE8 functional at fixed admissible snapshots.

## 6. Next research obligation

CE7 improves the actual whole-continuum bound from O(epsilon) to
O(epsilon^2[1+log(1/epsilon)]). Removing the logarithm requires further
information about the complete source. The exact threshold functional
CE8 specifies what a cancellation argument must control.

Next study the temporal modulation of gamma_1,gamma_2 along the actual
prepared evolution and the low-frequency remainder around k=0. Retain
velocities, source-time phases and cascade terms. A zero instantaneous
coefficient and sufficient spectral regularity are one possible route.
A time-integrated cancellation estimate is another. Critical O(epsilon^2)
closure, three-channel modal-action errors and independent proof review
remain open.
