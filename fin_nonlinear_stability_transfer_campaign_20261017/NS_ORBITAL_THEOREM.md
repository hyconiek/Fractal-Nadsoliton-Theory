# NS-003 — quantitative nonlinear orbital stability, all fields and both Cauchy coordinates

**Conditional theorem; not verified for the actual FIN coefficient tuple.** Adopt exactly the frozen action/phase space and *entire* P_NS region of NS_CONTRACT and NS_GAP_THEOREM. Suppose V is a real polynomial with its inherited coefficient signs, U has the stated genuine kink, and the nonlinearity has local H1 -> L2 regularity (automatic for a polynomial potential in 1+1). Assume the unique energy-class solution enjoys the standard continuation alternative, i.e. it can be continued while its H1 perturbation norm and L2 velocity norm remain bounded. For smooth data, persistence of higher regularity is used in NS_TRANSFER_THEOREM.

## Modulation instead of freezing translation

For positions (phi,z) near the translated kink manifold, the implicit function theorem applied to F(a;phi)=<phi-K_a,K_a'> gives a unique nearby continuously varying a=a(phi) with

<eta_phi,K_a'>=0, eta_phi=phi-K_a, eta_z=z.

At eta=0, the derivative of F with respect to a is +||K'||_2^2 (up to a sign depending on F convention), so it is nonzero. The modulation map is locally Lipschitz H1 -> R. The norm of the modulated full phase displacement and the infimum orbital distance are equivalent for a sufficiently narrow tube. All 24 real z components and both velocities remain in the norm. Nothing constrains the kink's centre or requires its velocity to vanish.

The orthogonality condition also controls its *time derivative* without pinning the centre: differentiating it gives

  a_dot [||K'||_2^2 - <eta_phi,K_a''>] = -<phi_t,K_a'>.

In the fixed tube |<eta_phi,K_a''>|<=||K'||_2^2/2, so |a_dot|<=2||phi_t||_2/||K'||_2. The centre may move, but its speed is controlled by the full retained Cauchy velocity. No artificial constraint a_dot=0 is applied.

## Expansion and a paid local nonlinear remainder

Translation invariance and the static equations make the first variation of E vanish. For eta=(eta_phi,z) and pi=(phi_t,z_t), with <eta_phi,K_a'>=0,

Delta E = E(phi,z,pi)-E(K,0,0) = (||pi||_2^2+Q[eta])/2 + R_3[eta],
Q[eta]=<eta_phi,L_phi eta_phi>+<z,L_perp^24 z>.

For ||eta||_{H1}<=r, the one-dimensional H1 embedding and the polynomial Taylor formula yield

|R_3[eta]| <= C_3(r) ||eta||_{H1}^3,

where C_3(r) can be computed as (C_S/6) times the supremum of the third potential derivative on the bounded segment between (K_a,0) and (K_a+eta_phi,z), with the finite integral of |eta|^3 controlled via L-infinity x L2^2. This is valid despite negative quartic/sextic coefficients: it is LOCAL, not a claim of globally bounded potential. The coefficient dependence is NOT replaced by an attractive potential fit.

Let c_H=delta/(delta+1+C_V) from NS_GAP_THEOREM. Choose r>0 small enough that C_3(r) r <= c_H/4. In the modulation tube, the full relative energy is controlled from both sides:

c_E (||eta||_{H1}^2+||pi||_2^2) <= Delta E <= C_E (||eta||_{H1}^2+||pi||_2^2),

where c_E=min(1/2,c_H/4)>0 and C_E is an explicit finite upper-form/Taylor bound from the frozen coefficients. This includes all complex phases, 24 real Psi velocities, the phi velocity and no artificially pinned zero mode.

## Global orbital estimate conditional on the supplied Cauchy theory

For sufficiently small initial full orbital distance d_0, specifically any d_0 whose modulated representative satisfies C_E d_0^2 < c_E r^2/16 (adjusted by the fixed modulation norm-equivalence constants), the energy identity and continuity prevent exit from the r-tube: an exit would require Delta E>=c_E r^2 whereas initial Delta E is strictly smaller. Local wellposedness, the continuation alternative and translation-invariant coefficients extend the solution for all t in R, staying in a finite-energy small tube. Consequently

sup_{t in R} d(X(t),O_K) <= C_orb d(X(0),O_K),

with C_orb <= C_mod sqrt(C_E/c_E), with C_mod the measured modulation/orbital norm equivalence. The bounds are conditional on the *actual* action satisfying P_NS; no claim of universal orbital stability for all legacy coefficient choices.

The kink position may drift. Translation is controlled by the orthogonality condition, not an added controller or conserved charge. This gives stability, not decay, radiation asymptotics, relaxation, scattering or convergence of the centre.

## Precise blocker if region is not certified

The above is a full conditional PDE proof. For the particular FIN legacy action, the still-open, typed input is a certificate of the exact heteroclinic wall K, the two wall-error suprema and the eigenvalue bounds of the actually supplied M2,K_s,Y. If delta_Psi<=0, e_scalar>=3k² or K fails to exist, this proof does NOT classify the model as unstable. One must either supply a true negative direction (including all real/imaginary blocks) or retain OPEN_FIN_SPECIFIC_ORBITAL_GAP. Exact kink substitution and full coefficients are more important than additional spectral digits.
