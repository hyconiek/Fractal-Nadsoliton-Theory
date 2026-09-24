# LAW-04 — stored relational flow and reversible full dynamics

**Disposition:** `PASS_AS_STRUCTURALLY_COMPATIBLE_KINETIC_CANDIDATE__SOURCE_OPEN`.

Use an extended state `(beta,pi)` with

`H = 1/2 pi^T M^-1 pi + 1/2 beta^T K beta + U_FIN(beta,...)`,
`dot beta=M^-1 pi`, `dot pi=-K beta-grad U_FIN`.

For fixed `M,K` this conserves `H` exactly. It therefore lies outside the dissipative gradient class killed by DYN-003/004. The old front identity remains important: with positive damping and equal wells it still forces zero speed. Setting the full law to reversible (`gamma=0`) only removes that obstruction; it does **not** prove a FIN soliton.

More importantly, this candidate naturally explains the memory found under dynamic refinement. Partition a conservative second-order system into visible `u` and hidden `v`. In Laplace frequency `s`, exact elimination gives

`D_eff(s)=s^2 M_u+K_uu-K_uv(s^2 M_v+K_vv)^-1 K_vu`.

This is precisely a frequency-dependent response object. For a two-oscillator control (`omega=1.0`, `Omega=3.0`, `g=0.8`) the reduced term is `g^2/(s^2+Omega^2)`, whose time-domain kernel is `(g^2/Omega) sin(Omega t)`. The normal-mode frequencies are 0.959575 and 3.013174. Thus a reversible full state can reduce to non-Markovian memory without inserting viscous damping.

This unifies the earlier refinement-memory obstruction with the proposed kinetic ontology: hidden stored flow/modes are one mathematically coherent *explanation* of the memory object. The remaining hard question is sourcehood of the reversible/symplectic law and of `M`, not algebraic consistency.
