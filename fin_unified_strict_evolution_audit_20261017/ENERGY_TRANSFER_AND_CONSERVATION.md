# One action, explicit energy exchange, and the closed-world boundary

Date: 2026-10-17. Conditional mathematical balances, not SI energy calibration.

## 1. The actual13-coordinate action has an exact model Hamiltonian

For fixed coefficients in L=||q_dot||²/2-V, the source momenta are p=q_dot
and H=||p||²/2+V. Along the source equations,

    H_dot=grad V.p+p.(-grad V)=0.

Let V_a be the self potentials and U_ab the pair potentials:
`U_ij=(K_s)_ij psi_i psi_j`, `U_i,Phi=y_i Phi²psi_i²`.
Assign e_a=p_a²/2+V_a+1/2 sum_b U_ab. Direct differentiation gives

    e_a_dot=sum_b J_ab,
    J_ab=(p_b partial_b U_ab-p_a partial_a U_ab)/2=-J_ba.

Thus internal coordinate energy exchange is derived from one supplied
action, with no separate force fitted to each phenomenon. Equivalently,
keep a separate E_interaction: self-energy losses plus its derivative
cancel exactly. That allocation is not unique; moving portions of U among
subsystems changes the allocated currents without changing total H.
No physical object identity or locality is selected by this bookkeeping.

## 2. The already written spacetime action supplies a conditional local flux

For a smooth solution of the legacy real scalar sector on its **supplied**
flat Lorentz domain,

    e=1/2 sum_a(q_a,t²+|grad q_a|²)+V(q),
    S=-sum_a q_a,t grad q_a,
    partial_t e+div S=0.

These follow by multiplying the Euler--Lagrange equations by q_a,t, not
by introducing Maxwell equations. Total energy is conserved when decay,
periodicity or a paid zero-flux boundary removes the surface term. Open
boundaries require the actual flux integral. The same component splitting
as above adds internal J_ab on the right of each component balance.
Geometry, boundary conditions and the action-to-physical-unit map remain
inputs. Covariant FAR scaffolding is not a source theorem for them.

## 3. Exact continuum-excitation certificate in the existing Phi4 sector

Use the legacy double-well sector and its kink from the structures document.
For odd eta write eta=a f+r, with normalized shape f, r perpendicular to f.
Translation is even and is excluded by this declared odd data class. With
omega1²=3mu²/2 and N(eta)=3lambda Phi_K eta²+lambda eta³,

    a_ddot+omega1²a=-<f,N(eta)>,
    r_tt+H r=-P_c N(eta).

The exact relative energy decomposes into

    E_shape=(a_dot²+omega1²a²)/2,
    E_cont=(||r_t||²+<r,Hr>)/2,
    E_int=lambda integral Phi_K eta³+lambda/4 integral eta⁴.

Its derivatives are respectively `-a_dot<f,N>`, `-<r_t,N>` and their
opposite sum. This is a diagnostic decomposition with a declared spectral
projection, not a unique measured split between matter and photons.

The continuum source is demonstrably nonzero. Put z=mu x/sqrt2,
f=c sech(z)tanh(z), c²=3mu/(2sqrt2), and S2=Phi_K f². Direct integration
after y=tanh z yields

    ||P_c S2||²=(9sqrt2 mu³/(8lambda))
                  [4/63-3pi²/512] >0.                (1)

The elementary integrals are
`integral_[-1,1] y^6(1-y²)dy=4/63` and
`integral_[-1,1] y^4 sqrt(1-y²)dy=pi/16`; f normalization supplies3/2
in the subtracted projection. Positivity needs only pi²<10, giving the
rational lower bracket79/16128 before the positive prefactor. One proof
of the pi bound is the positive integral
`integral_0^1 x^4(1-x)^4/(1+x²)dx=22/7-pi`.
No floating quadrature or assumed Fermi-golden-rule constant is used.

For initial eta=A f, eta_t=0, the continuous acceleration is
`-3lambda A² P_c S2+O(A³)`, nonzero for sufficiently small nonzero A.
Smooth local evolution then creates positive E_cont at order A^4 t².
For suitably small odd data the published local asymptotic-stability
theorem, with its hypotheses paid, further implies local relaxation to
the kink. Conservation forces positive excess energy outside fixed compact
neighborhoods. No complete scattering profile or decay rate is asserted.

Time reversal gives corresponding reverse-transfer histories in the same
Hamiltonian model. This is not a claim of generic perfect absorption,
irreversibility, a photon or a physical detector. The sector is legacy/
scaffold and kernel-blind; it does not meet the strict-source success gate.

## 4. Why the q12 parent alone is not a conserved-energy closed universe

For the actual finite native count chain, every legal single-copy move has
positive rate and the state graph is connected. If a state function E is
pathwise conserved for every start, it is constant across every positive
edge, hence constant. Even the weaker condition L E=0 forces a constant
on a finite irreducible chain, by the maximum principle. Reversibility also
gives

    <E,-L E>_pi=1/2 sum_(x,y)pi_x q_xy(E_y-E_x)².

It is strictly positive for every nonconstant E. Thus Gibbs stationarity
or E_pi[L E]=0 does not furnish a nontrivial conserved energy on the bare
native count state. This is not a no-go for a larger underlying Hamiltonian
carrier; it identifies exactly the missing state-map/environment/momentum
bridge. Reducible diagnostic classes are not silently covered by this lemma.

The reaction-trace construction in the ancillary package illustrates a
conditional open comparison, not a derivation of a closed-world total energy.
No statement here turns model-time conservation into joules/seconds or
exports physical gravity, gauge charges or L_total closure.
