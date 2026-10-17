# UE-003: energy, currents, two-configuration composability

For a **verified** canonical local field action with `q=(Phi,Psi)`, energy density and flux are

`e= (|q_t|²+|q_x|²)/2 + V(q)`, `j= -q_t·q_x`.

The field equation `q_tt-q_xx+grad V(q)=0` gives **pointwise** `partial_t e+partial_x j=0`, hence energy preservation if the spatial boundary flux vanishes. For a kink, use finite *relative* vacuum energy; the absolute spatial integral of a nonzero constant vacuum potential is divergent. A constant internal kernel matrix is a local species coupling, not an assumed physical spatial geometry or source of propagation speed.

If the same coupled action and its topology admit two compatible localized configurations, one may prescribe a **single pair of total field Cauchy profiles** (with the necessary global kink topological charge/boundary behavior) and evolve the same nonlinear PDE. Do not superimpose two solutions and call their sum a solution: nonlinearity generally prevents this. In 1+1, incompatible endstates or kink charges obstruct a naïve two-wall configuration. Smooth patches joined by cutoff and satisfying the same end conditions are admissible initial profiles, not exact bound states.

Under a local strictly hyperbolic canonical equation and standard smooth wellposedness, the domain of dependence gives finite propagation relative to the declared coordinates, with speed 1 in these *dimensionless PDE conventions*. It does **not** derive physical spatial locality or light speed. Boundary influx, ghost kinetic signs or nonlocal *spatial* kernels would invalidate the stated energy/causality application.

For the action split in `COUPLED_PROCESS_THEOREM_OR_OPEN_LEMMA.md`, `E_Phi`, `E_Psi` and `E_int` separately vary by work terms; only their sum is fixed. Projecting onto wall shape, discrete spectrum and continuum can be a useful *coordinate choice*, but those projections depend on background and linearization. Their share of the energy is not an additional canonical conserved charge. A local Taylor coupling/transfer witness is NOT general emission -> propagation -> absorption, even if an outgoing continuum projection is nonzero.

## Exactly paid same-field two-configuration cross control

In the **illustrative**, not source-bound scalar action with `V_int=Phi² Psi²`, set initially `Phi_0=1` and `Psi_0(x)=f_+(x)+f_-(x)`, where `f_±(x)=exp(-(x±a)²)`, `a>0`. Take any compatible smooth finite-energy velocities. Both excitations are profiles of the **same single Psi field**, not independent particle species and not an added two-body potential. Because the unchanged action is quadratic in this field at fixed `Phi_0`, its interaction energy includes the exact cross term

`E_int[f_++f_-]-E_int[f_+]-E_int[f_-] = 2∫f_+ f_- dx = sqrt(2*pi) exp(-2a²)>0`.

The identity is checked symbolically. It does **not** imply a force law between persistent solitons, a stable composite, or scattering: the full action's remaining kinetic and potential cross terms matter, and the profiles need not persist under evolution. It is only a transparent same-field composition compatibility test.
