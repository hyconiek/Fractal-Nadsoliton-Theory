# Synthesis — tested missing-law architecture

The four proposals are **not equally strong**, but none was killed outright.

1. **Mediator neutrality** is a mathematically clean conditional rank-source: it reproduces A7 exactly on the current circulant input and makes distinct predictions away from it. Its remaining arbitrary datum is now sharply isolated: why neutral through multipole order 2?
2. **Geometric measure + face transmission** is the strongest concrete correction. It exactly repairs the prior rectangular GEO-002 anisotropy counterexample without requiring ergodicity of the global aspect variable. General FIN-derived irregular geometry remains open.
3. **Geometry vs excitation** is confirmed, with a correction: `beta` is one transverse invariant, not the whole excitation sector. Already `(phi3,phi4,phi5)` needs a second invariant `gamma`; both vanish on the locked manifold.
4. **Stored flow / reversible full state** is structurally compatible with the forced memory results. Coarse memory appears automatically by exact elimination of hidden conservative modes. This removes the specific dissipative no-drift theorem from the full-state category, but does not supply localized matter.

## Revised architecture

A better working object is

`intracell mediator + cell geometry/measure + translation coordinate alpha + transverse deformation vector r + conjugate flow pi`.

The word **vector** matters: current evidence does not justify collapsing the internal deformation sector to one scalar beta.

A natural conditional graph Hamiltonian is

`M=diag(mu Vol(V_i))`,
`K=B diag(kappa Area(F_e)/ell_e) B^T`,
`H=1/2 pi^T M^-1 pi + E_FIN(r,alpha) + intercell coupling`.

Only the geometry/transmission form and the algebraic Hamiltonian consequences have been tested here. `E_FIN`, the kinetic bracket, the multipole cutoff, and physical interpretation remain source obligations.
