# LAW-03 — geometry coordinate versus internal excitation

**Disposition:** `PASS_WITH_ONE_MISSING_INTERNAL_COORDINATE`.

For phases `(phi3,phi4,phi5)`, define

- `alpha = phi4-phi3`,
- `beta = phi3-2 phi4+phi5`,
- `gamma = 4 phi3-3 phi4`.

The integer transform has determinant 1 and sends a carrier translation `(phi3,phi4,phi5)->(phi3+3d,phi4+4d,phi5+5d)` to `(alpha,beta,gamma)->(alpha+d,beta,gamma)`. Thus alpha is exactly translation-covariant while beta and gamma are internal invariants. The inverse is

`phi3=3 alpha+gamma`, `phi4=4 alpha+gamma`, `phi5=5 alpha+beta+gamma`.

The all-orders FIN locked manifold `phi_k=k a` gives **beta=gamma=0**. Consequently beta is a legitimate transverse/hard-mode deformation, but it is not a second Goldstone mode. Also, alpha+beta alone are not a complete coordinate chart for these three phase variables: an independent invariant gamma remains unless another theorem freezes or slaves it.

There is a separate exact warning supporting the geometry/excitation split. If alpha itself defines positive edge lengths `ell=Delta alpha` and one reuses alpha as the field in the already-derived `c(ell)=kappa/ell` energy, then each edge contributes `(kappa/2) ell`; at fixed circumference the total is constant. The identity map carries geometric length/tension, not an ordinary propagating field stiffness.

The refined architecture should therefore be `geometry alpha + transverse deformation vector r`, not yet `geometry alpha + one scalar beta`, unless FIN supplies a mode-selection/slaving law.
