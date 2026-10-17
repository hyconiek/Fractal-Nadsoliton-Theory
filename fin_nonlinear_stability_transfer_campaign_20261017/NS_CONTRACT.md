# NS-001 — conditional radial-action contract, not a FIN-source promotion

Date 2026-10-17. All mathematical statements below are about the *already-posited* canonical radial legacy PDE, conditional on verification of its exact displayed coefficient conventions. No new fundamental action, field, fit, dimensional constant, bath, or physical interpretation is adopted.

## Degrees of freedom and signature

Space-time: supplied 1+1 Minkowski R_t x R_x, wave speed normalized to 1, with no boundary flux at spatial infinity. Scalar phi is real; Psi=(psi_1,...,psi_12) is complex. Write z=(Re Psi, Im Psi) in R^24; z and its 24 velocity components must never be discarded. Canonical real kinetic terms and spatial gradients are assumed, not derived from the strict nadsoliton. The affine finite-energy phase space for a kink K is

X_K = ((K,0)+(H^1(R;R) x H^1(R;R^24))) x (L^2(R;R) x L^2(R;R^24)).

Here the first pair are positions, the second pair independent velocities. Translation a acts on every position/velocity field; K_a(x)=K(x-a). No physical units are inferred.

The conditional action density is L = (|phi_t|^2+|z_t|^2-|phi_x|^2-|z_x|^2)/2 - V(phi,z), with V the *supplied* radial legacy polynomial, not a polynomial reconstructed here. We retain its sign conventions symbolically. Write exactly its Taylor decomposition at z=0:

V(phi,z) = U(phi) + z^T A z/2 + phi^2 z^T Y_24 z + R_4(phi,z),
A=M_2+K_s real symmetric in the 12-complex basis, Y real symmetric in that basis, A_24=A direct-sum A, Y_24=Y direct-sum Y; R_4=O(|z|^4) at fixed bounded phi. If the supplied sextic includes terms in phi, they remain in U. If R_4 includes mixed higher terms, keep them unmodified. If actual action contains additional *quadratic* phi-dependence, mixed phi-z terms at z=0, derivative couplings or non-canonical kinetic metrics, this decomposition is inapplicable: an explicit corrected Hessian is required. All signed self-masses and sextic coefficients are retained as given; no positivity of a sextic sign is imposed to secure a result. Concretely, when U(phi)= -m_phi^2 phi^2/2 + lambda_phi phi^4/4 + g_phi phi^6/6 + C, the signs in this formula are definitions and m_phi^2, lambda_phi, g_phi are NOT assigned physical numbers. Likewise any signed Psi quartic and sextic remain in R_4.

The full PDE is phi_tt-phi_xx+V_phi(phi,z)=0 and z_tt-z_xx+grad_z V(phi,z)=0. Differentiating V twice at (K,0) gives L_phi=-d_x^2+U''(K) and L_perp^24=[-d_x^2 I_12+A+2K^2Y] direct-sum the same real 12-block. Cross Hessian vanishes only because the fixed radial action is even in z.

## Kink stationarity and exact wall relation

K is an actual monotone classical heteroclinic solving -K''+U'(K)=0, K(-infty)=-v, K(+infty)=v, K'>0, U(+-v)=0 and U(eta)>0 for -v<eta<v. Thus K' = sqrt(2U(K)) and K''=U'(K) after the first integral. Assume K,K',K'',U''(K) smooth/bounded, K' in H^1 and both limiting vacuum curvatures positive. The kink is found by this quadrature, not postulated as tanh if a phi^6 term is present.

Important exact obstruction: if U(phi)=-m^2 phi^2/2+lambda phi^4/4+g phi^6/6+C and K=v tanh(kx), coefficient comparison in K''=U'(K) forces m^2=2k^2, lambda=2k^2/v^2, g=0. Hence a genuine nonzero phi sextic cannot share the exact phi4 tanh wall without additional terms cancelling its quintic force. The inherited tanh/Poeschel--Teller calculation is a bound/control, not automatically the actual sextic-wall spectrum.

## Premise/source dispositions

Premises: canonical 1+1 space-time and signature; 12 complex-plus-1 real fields and momenta; exact radial polynomial and its signed parameters; K_s real symmetry and kernel/legacy interpretation; positivity of kink potential along the heteroclinic; constant coefficient/energy convention; Cauchy class and uniqueness; no-flux at infinity. None is a strict-source physical derivation. The previous UE source register explicitly said full QW action arrays were not verified in that run. The newly asserted source-access binding is not independently materialized in the present container: no claim to having checked actual QW bytes or computed FIN-specific gap is made, and no repetitive source inventory is issued. A typed nonlinear coefficient/action map excluding fixed-K/same-Hessian alternatives is still required for a strict claim.
