# UE-003: one exact conditional local nonlinear transfer theorem

## Action class, conditional and source-bound

Let `Phi` be real scalar and `Psi` an n-component real field (or a rigorously licensed invariant real slice of a complex action) in 1+1 dimensions. Assume **the same pre-existing smooth action** has canonical positive kinetic terms and the decomposition

`V(Phi,Psi) = U(Phi) + W(Psi) + V_int(Phi,Psi)`.

Do **not** infer this precise separability or any special `V_int` from the supplied transverse Hessian alone. The theorem applies only after it is verified for the *actual* coupled action. Let fields solve

`Phi_tt - Phi_xx + partial_Phi V = 0`,
`Psi_tt - Psi_xx + grad_Psi V = 0`.

Assume sufficient smoothness, local Cauchy wellposedness, decaying perturbations/finite relative energy, and no incoming/outgoing boundary flux over the interval considered.

Define the diagnostic energies

`E_Phi = ∫[Phi_t²/2 + Phi_x²/2 + U(Phi)] dx`,
`E_Psi = ∫[|Psi_t|²/2 + |Psi_x|²/2 + W(Psi)] dx`,
`E_int = ∫V_int(Phi,Psi) dx`.

The full conserved quantity is `E_tot=E_Phi+E_Psi+E_int`; splits are diagnostic, not independent charges. With `F=partial_Phi V_int` and `G=grad_Psi V_int` the exact identities are

`dE_Phi/dt = -∫F Phi_t dx`,
`dE_Psi/dt = -∫G·Psi_t dx`,
`dE_int/dt = ∫(F Phi_t+G·Psi_t) dx`.

## Exact conditional theorem (nonzero Psi, no added interaction)

Fix **any** compatible smooth finite-energy Cauchy profiles `Phi_0,Psi_0`, where `Psi_0` is not identically zero and `V_int` is the *existing* action's interaction. Suppose

`F_0=partial_Phi V_int(Phi_0,Psi_0)` and `G_0=grad_Psi V_int(Phi_0,Psi_0)`

are smooth compactly supported, with `A=∫F_0² dx>0` and `B=∫|G_0|² dx>0`. Choose compatible initial velocities

`Phi_t(0)=delta F_0`, `Psi_t(0)=-delta(A/B)G_0`, `delta>0`.

These are Cauchy data of **the same theory**, not an external forcing, controller, bath or modified law. Direct substitution proves

`E_Phi'(0)=-delta A<0`,
`E_Psi'(0)=+delta A>0`,
`E_int'(0)=0`, `E_tot'(0)=0`.

Continuity of the derivatives for the locally wellposed solution implies a **strict increase of E_Psi and decrease of E_Phi over a sufficiently short positive time**. This is an existence theorem for a carefully specified initial state, not a statement about all data or persistent emission/absorption. Crucially, the interaction energy is stationary at the initial instant, so instantaneous transfer from the Phi diagnostic sector to the Psi diagnostic sector is genuine under the chosen split.

### Exact synthetic control (not FIN coefficients)

For `U(Phi)=(1-Phi²)²/2`, `W(Psi)=m²Psi²/2`, `V_int=Phi²Psi²`, let `Phi_0=tanh x`, `Psi_0=sech x` and define `F_0=2 tanh x sech² x`, `G_0=2 tanh² x sech x`. These decay rapidly (not compactly supported, but are in the required Sobolev spaces); the theorem also applies with rapid decay, by density and vanishing boundary terms. The exact integrals are

`A = ∫F_0² dx = 16/15`, `B = ∫G_0² dx = 8/5`, `A/B=2/3`.

With `Phi_t=delta F_0`, `Psi_t=-(2 delta/3)G_0`, the three exact derivatives are `-16delta/15`, `+16delta/15`, `0`. The local semilinear system has standard smooth local Cauchy evolution. None of `m`, this illustrative potential, or these profiles is advertised as sourced by FIN; they serve only to **replay the structural theorem**.

## What is NOT concluded

No spontaneous radiation under generic Cauchy data, particle formation, asymptotic scattering operator, stable two-wall molecule, absorption by a second isolated wall, probability of interaction, photon mass, or cross section follows. Nor is the diagnostic decomposition unique. No full-complex-action parity or source bridge has been checked.
