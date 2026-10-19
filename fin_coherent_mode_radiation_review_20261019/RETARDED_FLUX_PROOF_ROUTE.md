# A concrete continuation of the still-open RP-004 gate

Date: 2026-10-19. The reduction below is exact for the frozen submitted
family and preparation. The logarithmic-window and flux estimates remain
proof targets, not results or new microscopic laws.

## 1. Pay the prepared invariant reduction

Use the displayed bound radial action, not an arbitrary even R4 extension.
All K_1j vanish, Y is diagonal, self terms depend on each |Psi_i|, and other
species and velocities start at zero. Uniqueness keeps them zero. Complex
conjugation symmetry keeps the active neutral Psi real for its real initial
data. This is an exact invariant Cauchy class, not deletion of arbitrary
imaginary modes in the stability theorem.

Write the active field as psi and its original signed self coefficients as
g4,g6. The full equations on this particular class are exactly

    phi_tt-phi_xx+2phi(phi²-1)+6phi psi²=0,
    psi_tt-psi_xx+(a+6phi²)psi+g4 psi³+g6 psi⁵=0.

No h,d_j or other spectator coefficient enters. The existing datum is
phi=K=tanh x, phi_t=0, psi=epsilon sqrt(3/4)sech²x, psi_t=0. Reflection
equivariance and uniqueness preserve phi odd and psi even. Hence the exact
prepared scalar perturbation is odd and has no even translation-zero-mode
component. This is initial symmetry, not a pinning controller or a theorem
about general perturbations. It simplifies the remaining time estimate.

## 2. Exact remainder equations to use, not an appended wave law

Set eta=phi-K, L_phi=-partial_x²+4-6sech²x and
L_a=L_phi+a+2. Then

    eta_tt+L_phi eta=-6K psi²-6eta psi²-6K eta²-2eta³,
    psi_tt+L_a psi=-12K eta psi-6eta²psi-g4 psi³-g6 psi⁵.

Let v=f cos(omega t), omega²=a+2, and u solve the already frozen retarded
second-order equation u_tt+L_phi u=-S(1+cos(2omega t)), zero Cauchy data.
For r=psi-epsilon v and w=eta-epsilon²u the exact sources are

    w_tt+L_phi w=-6K(psi²-epsilon²v²)
                  -6eta psi²-6K eta²-2eta³,
    r_tt+L_a r=-12K eta psi-6eta²psi-g4 psi³-g6 psi⁵.

There is no independent incoming scalar packet. All free initial r,w data
are zero. Positive operator energy norms on the prepared parity class and
the NS small-data bounds are available inputs. A suitable future bootstrap
should bound r,w in a norm strong enough for the selected boundary flux,
uniformly up to a prescribed logarithmic time. The source's O_T statements
do not prove such a bound.

## 3. Freeze a target window only after validating its estimates

Here q in[1,sqrt3], Omega=sqrt(q²+4), so the harmonic group velocity
q/Omega is at least1/sqrt5>2/5. A **candidate**, not a certified window,
is R_epsilon=2log(1/epsilon), T_epsilon=12log(1/epsilon). It leaves a
strict harmonic travel-time margin, but massive retarded tails do not have
one sharp group-delay front. Full light-cone propagation, turn-on, dispersive
error and bound-state contributions must still be paid. The pair must not
be called licensed merely because its geometric inequality is positive.

Prove one finite-time retarded outward-flux lower bound and one uniform
nonlinear error estimate that jointly dominate kink/bound-mode tails,
initial slip and any discrete shape contribution. If the proposed pair
fails, report the failed estimate; do not optimize a sweep of windows.
No actual outward-flux result is added by this route document.

## 4. Exact success and failure semantics

The actual current is J=-phi_t phi_x-psi_t psi_x in this invariant class.
Positive integrated J at the two endpoints must be proved for the actual
nonlinear solution, not only the forced comparison or a post-subtracted
bound-mode signal. A new logarithmic-time remainder theorem by itself is
useful but leaves the flux sign open. Gamma>0 supplies a coefficient, not
this missing finite-time inequality. Preserve NS and the spectral RP result
even if the present route fails.
