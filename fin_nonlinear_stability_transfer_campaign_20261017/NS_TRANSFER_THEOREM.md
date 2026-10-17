# NS-004 — a single stable trajectory with strictly nonzero local sector transfer

## Sector diagnostics taken from the same fixed action

In the already supplied radial polynomial, split the potential exactly as V(phi,z)=U(phi)+T(z)+W(phi,z), where T(z)=V(0,z)-U(0) (constant convention adjusted once), and W=V-U-T. The quadratic Phi-Psi mixing must obey W(phi,z)=phi^2 z^T Y_24 z + O(|z|^4) modulo any phi-independent quadratic moved to T. All parts are the unchanged terms of V, not new charges. Define:

E_phi = integral [phi_t^2/2 + phi_x^2/2 + U(phi)]dx,
E_Psi = integral [|z_t|^2/2 + |z_x|^2/2 + T(z)]dx,
E_int = integral W(phi,z)dx,
E_total = E_phi+E_Psi+E_int.

Subtract the vacuum constants consistently when defining finite energies (E_phi is finite in the kink class, and T(0)=W(phi,0)=0). Integration by parts using the full PDE, vanishing boundary flux and finite-energy regularity gives, *exactly*,

E_phi'=-int W_phi phi_t dx, E_Psi'=-int <grad_z W,z_t> dx,
E_int'=int (W_phi phi_t+<grad_z W,z_t>) dx,
E_total'=0.

The split is analytic bookkeeping, not a unique conserved particle energy, physical detector output or a new invariant.

## Nonzero transfer within a stable orbit: constructive theorem

Take any frozen coefficient tuple in P_NS with Y != 0 (as a real symmetric matrix), and choose real unit u in R^24 with c=u^T Y_24 u != 0. Choose smooth real compactly supported h,p with h not zero and support of h inside a region where K != 0, and set

a = int K(x) p(x) h(x)^2 c dx > 0,
b = int K(x)^2 h(x)^2 c dx != 0.

This is achievable without changing the potential: fix h near a point K != 0, then choose the *initial velocity profile* p with the sign of c K. Velocity is permissible free Cauchy data, not a source of coefficients or an external force.

At t=0 choose the full compatible Cauchy data

phi=K, phi_t=epsilon p, z=epsilon h u, z_t=epsilon q_epsilon h u,

where q_epsilon= - [int W_phi(K,epsilon h u) p dx] / [int <grad_z W(K,epsilon h u),h u> dx]. The denominator is 2 epsilon b+O(epsilon^3) and nonzero for all sufficiently small positive epsilon; q_epsilon=-epsilon a/b+O(epsilon^3). Hence every position and velocity field, including both velocities, is retained and small; z_t is NONZERO, of size O(epsilon^2), for sufficiently small epsilon. This ratio is selected solely to make the interaction-energy derivative zero initially, not to fit a field-law parameter.

The three work identities now give exactly

E_int'(0)=0,
E_phi'(0)=-epsilon int W_phi(K,epsilon h u)p dx = -2 a epsilon^3+O(epsilon^5)<0,
E_Psi'(0)=+2 a epsilon^3+O(epsilon^5)>0.

The common nonlinear orbital-stability theorem applies to the **same solution for all t**, with sup_t dist(X_epsilon(t),O_K)<=C epsilon. For sufficiently regular compact data, the energy-work derivatives are continuous, so there exists T_epsilon>0 such that both strict opposite inequalities persist for all 0<=t<=T_epsilon, and E_total(t)=E_total(0). For a uniformly bounded H2 x H1 local solution of norm O(epsilon) on [0,T0], the work-derivative Lipschitz estimate |d/dt(E_sector')|<=C_work epsilon^2 permits the explicit lower time choice

0<T_epsilon<=min(T0, a epsilon/(2 C_work))

when the leading coefficient is bounded below by a epsilon^3. Therefore |E_Psi(T_epsilon)-E_Psi(0)| >= a epsilon^3 T_epsilon/2 and the Phi diagnostic sector falls by at least the same lower bound. This is a positive *finite* interval, typically order epsilon, not an asymptotic emission/absorption statement.

## Strict limitation and adverse outcomes

The theorem needs nonzero actual Y in the same action, the radial-even leading interaction, actual admissible K, the P_NS spectral inequalities and the canonical energy-class PDE. If actual Y=0 or the signed sextic forces no kink/negative mode, this particular proof cannot conclude. An arbitrary earlier UE transfer witness does not prove stable transfer without those paid hypotheses. No particle, photon, quantum measurement, pair scattering, bound state, numerical energy scale, ToE or empirical validation follows.
