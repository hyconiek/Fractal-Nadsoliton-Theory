# MT-003 — Joint local Gibbs critical/Gaussian limit

**Status:** asymptotic argument for conditional finite-N Gibbs types in the fixed MT-001 connected sublevel; no replacement by the full Gibbs measure. Model q=12, N∈12ℕ, bounded fixed λ, g_N=g_c+λ/√N. Assumptions and local geometry are inherited from MT-001, not new physical preparations.

## The actual count law

For count vector n, p=n/N,

 π_{D,N}(n)=[1/W_N(D_N)] [N!/(12^N∏ n_i!)] exp[(g_N N/2)p^T A p] 1{n∈D_N}.

The multinomial factor is indispensable. With p_m=u+mσ/12, w=p−p_m, fix an orthonormal q×10 matrix E spanning {1,σ}⊥ and define

 Z=N^(1/4)m,   Y=√N E^T w.

The **joint** candidate density relative to dz dy is

 ρ_λ(z,y) = Z_*^{-1} exp[χλ z²/2−z⁴/12−(1/2)y^T H_fast y],
 H_fast=E^T(12 I−g_c A)E.

H_fast≻0 from the inherited gap; the ten-dimensional Gaussian has covariance H_fast^{-1}, including the four A-null Fourier modes, which have Hessian eigenvalue 12. The marginal of Z is the **quartic** invariant law of the UF critical diffusion. The product here is asymptotic only; finite-N p/detector label laws need not factorize.

## Taylor calculation without invented independence

Because ∂_w F_g(p_m)=0 exactly for all m and g, the free-energy expansion has no term linear in w or mixed m²w, including at criticality. On bounded z,y,

 N F_{g_N}(p_m+w)
 =z⁴/12−χλ z²/2 +(1/2)y^T H_fast y+O_{K,Λ}(N^{-1/4})

(the displayed O can be improved; it is deliberately conservative). Explicit inputs: I(m)=m²/2+m⁴/12+O(m⁶), g_Nχ=1+χλ N^-1/2, H_W(m,w)=H_fast+O(|m|+||w||+N^-1/2), and ||w||=O(N^-1/2). Thus m⁶N=O(N^-1/2) and N||w||²|m|=O(N^-1/4), while all possible linear-fast terms vanish by the parity fixed-subgroup identity.

## Lattice, normalization and uniform integrability

Interior Stirling, uniformly on shrinking central domains with every p_i≥1/24, gives

 μ0_N(n)=(2πN)^(-11/2) (∏p_i)^(-1/2) exp[−N D(p||u)] (1+O(1/N)).

The affine lattice {n:Σn=N} has fixed covolume in the 11-dimensional tangent. Rescaling along the parity mode by N^(1/4) and along ten fast coordinates by N^(1/2) creates mesh tending to zero in all 11 dimensions. The Stirling prefactor tends uniformly to its constant at u on bounded scaled sets. Consequently Riemann sums there approach the displayed joint continuous integral, after the common mesh/covolume normalization cancels.

For tails, MT-001 gives *globally throughout D_N*

 N F_{g_N}(p) ≥ (1/24)Z⁴ +.04||Y||²−C_Λ.

There is also a matching local Taylor upper bound around u, hence a lower bound on W_N(D_N) from any fixed nonempty compact box in (Z,Y). For the remaining tail, bound each rescaled lattice cell by a constant times the integral of exp[−cZ⁴−c||Y||²] over its neighbouring continuum cell: the exponent is uniformly locally Lipschitz at scaled mesh size and polynomial bounds away from the central box are controlled by its tails. All p in D_N stay in the fixed positive interior cylinder, so Stirling prefactors are bounded above/below relative to their central value. The resulting dominating function is integrable, with arbitrarily many polynomial moments. Standard dominated Riemann-sum comparison therefore gives both weak convergence and uniform integrability of any fixed polynomial moment in (Z,Y).

Explicitly, for every continuous φ of polynomial growth, under π_D,N,

 E φ(N^(1/4)m,√N E^T w) →∫φ(z,y)ρ_λ(z,y) dzdy.

This includes Var_π_D(m)=N^-1/2 (Var_ρ(Z)+o(1)), with Var_ρ(Z)>0, which MT-004 uses. No finite-N exact independence, tagged-copy reconstruction, unconditional full-Gibbs quartic marginal or global mixing assertion is made.

## Component and conditioning payment

The canonical-path connectivity result MT-001 covers the full local sublevel for N large, in particular all count types in any fixed compact (Z,Y) box. The fixed outer boundary carries weight at least exp[−Nb] whereas the central normalization is only polynomial in N; it has vanishing conditioned mass. This pays the component/restricted-domain issue. The Gibbs normalization can differ exponentially from the FULL Gibbs normalization because an independent global competitor has lower F. This proof never substitutes one for the other.
