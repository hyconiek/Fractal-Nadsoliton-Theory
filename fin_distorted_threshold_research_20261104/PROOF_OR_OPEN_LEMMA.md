# DT: threshold-reference cascades with the actual potential and projections

Date: 2026-11-04. Base: 85fb974fdb4e7dc7962addea13973e85ed98e247.
Exact operator/reference calculations, not an actual nonlinear trajectory
or a completed cubic trace estimate. Self-audited; no formal verification.

## 1. The step beyond CS's homogeneous far field

CS isolated the actual cubic trace gate, but its coefficient A_self was
computed in a homogeneous far-field reference. Here retain the exact
operators H=-partial_x^2+4-6sech^2(x), L=H+w, K=tanh(x), and the original
continuous projections. We test the exact even threshold profile

    r0=(3K^2-1)/2=1-(3/2)sech^2(x),
    L r0=m^2 r0, m^2=4+w in [21/4,23/4].              (DT1)

r0 is bounded, not L2, and <f,r0>=0. It is an analysis profile, NOT a
new initial field. Substituting it into a stationary response cannot
replace the prepared finite-energy solution. The purpose is to expose
potential/projection corrections before any actual cubic-history claim.

Let s=sech^2(x), z=m^2, and S=K r0^2. The odd quadratic source generated
by the threshold trial u_1=r0(b exp(-imt)+conj(b)exp(imt))/2 is -6S times
its scalar squared amplitude. It contains a nondecaying part S~K.
Resolvents applied below to S are defined by an explicit subtraction,
not by pretending S belongs to L2.

## 2. Exact static response in the full odd sector

A bounded odd solution of H U_0=S is

    U_0=(5/8)K-(3/8)K^3-(1/2)x sech^2(x).            (DT2)

Verification uses

    H K^n=-n(n-1)K^(n-2)+(2n^2-2)K^n
                                      +(6-n(n+1))K^(n+2),
    H(x sech^2 x)=4K sech^2 x.

They give H U_0=(1/4)K-(3/2)K^3+(9/4)K^5=S exactly.
The asymptotic particular part is K/4; U_0-K/4 is exponentially decaying
up to the harmless x factor. Uniqueness among bounded odd solutions
follows from the zero-energy homogeneous equation: its decaying solution
sech^2(x) is even, while the independent odd solution grows at infinity.
This is an operator identity, not smallness of the actual odd field.

## 3. The second harmonic needs a localized outgoing correction

For the selected second harmonic set E=4z and

    a_z=1/(4-4z)=-1/[4(z-1)],
    J_z=[-3-1/(z-1)]K s+(9/4)K s^2.                  (DT3)

Since (H-4z)K=(4-4z)K-4Ks, direct substitution gives

    S-(H-4z)(a_z K)=J_z.

J_z is exponentially localized and odd. Define the generalized response
unambiguously by

    U_2=a_z K+(H-4z-i0)^-1 J_z.                       (DT4)

The inverse on J_z is the full odd outgoing boundary value, including
its bound component. It is a weighted local response, never a bounded
inverse on global L2 at this energy. Equation (H-4z)U_2=S follows from
the definition. Its outgoing correction is not assumed to decay in space.

Use ET/CD's exact transforms with D(k)=(1+ik)(2+ik) and
S2(k)=pi k/sinh(pi k/2):

    M_(K s)=ik(k^2+4)S2/(8D),
    M_(K s^2)=[(k^2+4)/6]M_(K s).

At the open odd momentum q_z^2=4z-4 in [17,19],

    M_Jz(q_z)=C(z) M_(K s)(q_z),
    C(z)=-3-1/(z-1)+(3/2)z.                          (DT5)

C'(z)=3/2+1/(z-1)^2>0, and C(21/4)=631/136>0. All other factors in
M_(K s)(q_z) are nonzero. Thus the exact supplied potential produces a
strictly nonzero on-shell localized residual throughout the parameter band:

    d_Jz=Im<J_z,(H-4z-i0)^-1J_z>
                                  =|M_Jz(q_z)|^2/(2q_z)>0. (DT6)

This uses CD's both-root dk/(2pi) normalization. It is a spectral density
for a reference source, not emitted power of the actual nonlinear state.
The one-dimensional outgoing Green representation gives wave amplitudes
proportional to the corresponding Jost transforms at +/-q_z. Since these
are nonzero, U_2-a_z K has oscillatory outgoing tails and is not L2.
The source J_z is localized; its on-shell response is not a localized
correction to the homogeneous reference. This distinction is essential.

## 4. Keep the shape projection explicitly

The normalized odd bound profile is g=sqrt(3/2)sech(x)K. Its overlap with
the generalized source is an ordinary convergent integral:

    d_S=<g,S>=17pi sqrt(3/2)/64.                      (DT7)

With y=tanh(x), divide by pi sqrt(3/2) and use the arcsine moments to get
(1/4)[9(5/16)-6(3/8)+1/2]=17/64. Since H g=3g,

    P_o U_0=U_0-(d_S/3)g,
    P_o U_2=U_2-[d_S/(3-4z)]g.                        (DT8)

Integration by parts is legitimate here: g and its derivatives decay
exponentially, while the reference responses and derivatives are bounded.
The subtraction affects only a localized bound profile; it cannot remove
the outgoing continuum tails or change the on-shell density. The projected
equations are H(P_o U_0)=P_o S and (H-4z)(P_o U_2)=P_o S.
Shape terms remain separate rather than being silently omitted.

## 5. Direct cubic profile and finite localized threshold defects

The even continuous projection is P_e h=h-f<f,h>. Direct beta integrals
with f=(sqrt(3)/2)sech^2 give

    <f,r0^3>=2sqrt(3)/35,
    P_e r0^3=r0^3-(2sqrt(3)/35)f.                    (DT9)

The direct cubic profile is not just r0 plus the bound mode. Its remaining
localized defect has a nonzero threshold overlap:

    integral r0(P_e r0^3-r0)dx=-27/35.                (DT10)

For verification, integral(r0^4-1)dx=-132/35 and
integral(r0^2-1)dx=-3; subtract these and use <r0,f>=0.
These are convergent integrals of DECAYING DIFFERENCES. In contrast,
integral r0^4 dx diverges. No finite total threshold coefficient is
inferred by integrating an unsubtracted non-L2 profile.

## 6. A full-potential stationary jet, with its exact limitation

For the selected threshold trial, the quadratic static and second-harmonic
responses have coefficients -3|b|^2 U_0 and -(3/2)b^2 U_2. The coefficient
of b|b|^2 exp(-imt) in the physical stationary cubic source is therefore

    A_full(x)=18K r0 U_0+9K r0 U_2-(3g4/8)r0^3.       (DT11)

This keeps the full odd response. Replacing U_0,U_2 by their projected
versions gives the continuum-intermediate reference; their difference is
an explicit localized multiple of K r0 g. Applying the outer P_e changes
only a localized f component. These operations preserve the outgoing
oscillatory tail in DT11. All such subtractions are defined by convergent
pairings with the exponentially localized bound profiles.

Inserting DT2/DT4 into DT11 separates

    A_full=A_self(w,g4)+9K r0 (H-4z-i0)^-1 J_z+D_loc,  (DT12)

where D_loc is exponentially localized up to a polynomial x factor and
A_self=9/2-9/[4(w+3)]-3g4/8 is exactly CS's homogeneous coefficient.
For example the constant parts use K^2->1, r0->1, U_0->K/4 and
U_2's chosen particular part a_z K. The differences from these constants
are localized. DT12 is a decomposition, NOT a claim that A_full has a
pointwise constant limit: the middle term has nonzero oscillatory tails.

Thus the full-potential threshold reference is not the homogeneous
phase coefficient plus only localized errors. A resonance analysis must
retain the radiating intermediate harmonic as well. No choice of g4 can
cancel that oscillatory reference tail by the spatially constant or
localized direct cubic terms alone. This statement concerns DT11, not
a theorem about cancellation in the actual finite-energy trajectory.

DT11 is NOT asserted to be CS1's complete actual distorted cubic symbol.
It selects a generalized threshold input and stationary harmonics, while
the actual problem contains finite packets, continuum momenta, time
histories and CN's normal-form/projection corrections. Shape-mediated
localized terms also belong to the established source ledger. Matching
this reference calculation to that full symbol remains a required step.

## 7. Next actual-history obligation

Use the explicit subtraction DT3--DT4 to organize a finite-energy packet
calculation near threshold, with outgoing harmonic history retained.
Keep the packet width, frequency support and time window explicit and
bound errors before identifying a stationary coefficient with the actual
compact observation. A formal replacement u(t,x)=b(t)r0(x) is not allowed:
it changes the finite-energy preparation and has no paid global remainder.

The actual O(epsilon) whole local-continuum estimate from ML is unchanged;
the O(epsilon^2) cubic trace and all three modal action errors remain open.
No new initial radiation, drive, g4 retuning, physical source or empirical
FIN claim is introduced. Historical proofs/hashes are preserved.
