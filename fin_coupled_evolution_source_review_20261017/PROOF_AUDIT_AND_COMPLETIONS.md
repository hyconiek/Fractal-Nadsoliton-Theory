# UE proof audit: conditional admission and actual local source binding

Date: 2026-10-17. Analytical model mathematics, not mechanized proof,
independent peer review or physical validation. Source bytes are unchanged.

## 1. The operator theorems survive the audit

For real symmetric constant A,B and k>0, use the closed form on H1(R;R^n)
of L=-partial_x² I+V_inf-B sech²(kx), V_inf=A+B. Bounded symmetric
multiplication gives the H2 self-adjoint realization. Interpolating between
A and V_inf proves nonnegativity when both are positive semidefinite.
No simultaneous diagonalization is needed.

Let mu=lambda_min(V_inf), b=max(0,lambda_max B). Pointwise Loewner order
gives L >= [-partial_x²+mu-b sech²(kx)]I. Put b=k²s(s+1),s>=0. Factoring
with D_s=partial_x+sk tanh(kx) gives

    -partial_x²-b sech²(kx)=D_s^* D_s-k²s².

For b>0, sech(kx)^s attains the scalar lower energy. For b=0 the lower
edge is mu but the constant formal mode is not an L2 eigenfunction. The
source bound remains correct in this endpoint; do not claim an attained
ground state there. Failure of a lower sufficient bound is inconclusive.

The trial f=sech(kx) gives exactly

    ||f||²=2/k, ||f'||²=2k/3, int sech²(kx)f²=4/(3k),
    R(u)=k²/3+u^T(V_inf-2B/3)u.

A strictly negative value is a true negative form direction. When V_inf>0,
the decaying well is relatively compact with respect to the constant
operator: truncate to a compact interval, use compact H2-to-L2 restriction,
and make the exterior multiplication norm small. The essential lower edge
is lambda_min(V_inf)>0, so the negative spectrum contains an isolated
eigenvalue. This is linear instability of the wave linearization, not a
standalone nonlinear-instability theorem for every full action.

The submitted2x2 example has positive far masses and Rayleigh -4/3 exactly.
The intake recomputes the integrals and uses a separate3x3 noncommuting
positive-far control with a negative trial. Neither matrix is FIN data.

## 2. What local access closes, and what it does not

Actual QW2163/2165/2166, their archived reports, the seven predecessor audit
documents, model contract and the legacy TeX are available now. Their
independent previous44-check record and exact original hashes remain valid.
This pays today's structural dependency, not the producer's historical
access or an independently sourced value for any free coefficient.

For the actual QW polynomial, at Psi=0 the12-real transverse block is

    M2+K_s+2Phi²Y,

with M2=diag(m2_i), Y=diag(gY_i), K_s=(K+K^T)/2 and zero self diagonal.
The archived QW code does not supply a numeric strict K_s or select M2,Y.
Nor does its one-coordinate derivative supply a physical spatial metric.
The conditional field interpretation still comes from the separately
written legacy action. Its explicit signed/sextic convention maps remain
required: m2_i=-m_i² and g6_i=delta or3delta/4 in the two legacy sections.

For a nonconstant Phi=v tanh(kx) to solve the displayed Phi equation at
Psi=0, its parameters are not independent: m2_Phi=-2k² and
lambda_Phi=2k²/v². These are necessary coefficient identities obtained by
equating the independent tanh and tanh-cubed terms. They are conditions on
the supplied family, not new natural values. The generic matrix-operator
bounds are valid without a source wall; a full-wall claim needs this check.

## 3. A bounded completion: the full complex transverse blocks

For the **displayed radial legacy action**, with real symmetric K_s and
real coefficients, write Psi=u+i v and keep both velocities. At Psi=0 all
quartic/sextic Hessian terms vanish; the mixed Phi-Psi block also vanishes.
The24-real transverse Hessian is exactly diag(L_perp,L_perp). The Phi block
is its existing scalar Hessian. Thus the operator bound applies to both
copies under this explicit radial/real-kernel contract, without deleting
imaginary modes. This is a structural theorem, not a numerical FIN-wall
verdict or a theorem about unspecified phase-sensitive completions.

For a Hermitian kernel C+iD with D real antisymmetric, realification instead
gives [[C,-D],[D,C]]. The real slice need not be invariant. New phase terms
would likewise require their actual Hessian. The source's broad warning
remains correct; today's completion is scoped to the actual displayed
radial polynomial with real mixing, not every possible full theory.

## 4. Local transfer: correct initial-data existence statement

With V=U(Phi)+W(Psi)+V_int, the work identities follow by integration by
parts. Vacuum constants must be subtracted before defining finite energies.
For F=partial_Phi V_int, G=grad_Psi V_int and A=||F||²,B=||G||²>0,
the compatible initial velocities (delta F,-delta A G/B) give

    E_Phi'(0)=-delta A, E_Psi'(0)=delta A, E_int'(0)=0.

Continuity preserves the strict signs for some positive time. This proves
existence for chosen Cauchy data; it is not generic emission or a uniquely
defined physical energy partition. Rapidly decaying F,G work directly by
vanishing boundary terms; no unproved density upgrade is needed.
Without positivity conditions on the self potentials, these diagnostic
sector energies need not themselves be nonnegative. The derivative signs
do not supply a physical particle/radiation decomposition or an entropy law.

The actual displayed coupling V_int=Phi² sum_i gY_i|Psi_i|² supplies these
forces. For any nonzero y_j, take Phi=v tanh(kx), Psi_j=epsilon sech(kx),
other Psi initially zero. Then

    A=16 y_j² v² epsilon⁴/(15k),
    B=8 y_j² v⁴ epsilon²/(5k), A/B=2epsilon²/(3v²).

This binds a **symbolic nonzero-coupling family**, not a natural value of
y_j. Off-diagonal K can excite the other Psi coordinates after time zero;
the one-active-component initial datum is not claimed to define an invariant
subsystem. The full equations and total Psi-sector energy are retained.

For the canonical1+1 polynomial field family, the local regularity premise
can be paid on an affine background space: write q=q_ref+h, with a smooth
bounded reference having decaying derivatives and an H^{s-1} equation
residual, and h in H^s, h_t in H^{s-1}, integer s>=3. Polynomial nonlinearities
are locally Lipschitz there by the1D Sobolev algebra estimate. The inhomogeneous
linear-wave energy estimate and a short-time contraction give local unique
solutions. Higher smoothness pays the time derivatives used below. This is
conditional PDE wellposedness on supplied geometry, not a derivation of it.

## 5. Required correction: affine coarse state and the Taylor payment

The source writes C(Phi,Phi_t,...) in H^s x H^{s-1}, but its tanh example
has nonzero limits and is not in L2, hence not in that H^s. Use the fixed
affine state Phi-Phi_ref in H^s and retain Phi_t. This is an invertible
relabeling by a known background, not a new observable or law. Differences
between the two Phi paths are in L2, as required by the proof.

For two sufficiently smooth solutions, Taylor with the integral remainder
and a common third-derivative bound gives ||R_i(t)||<=C t³. Their initial
Phi and velocity agree and their acceleration difference has norm a>0.
Then separation is at least a t²/2-2C t³, hence at least a t²/4 for
0<t<=min(T,a/(8C)). The triangle inequality gives worst-case predictor
error at least a t²/8. The conclusion even holds for a deterministic
time-dependent predictor of those same initial coarse data; it need not
be an autonomous semigroup. Memory or a restricted preparation class is
not refuted. For kink data, smallness/regularity and this common T,C must
be paid, not inferred from a trajectory plot.

The Gaussian overlap control is a correct composition-energy identity,
not a soliton force, bound state or scattering result. Exact symbolic
checks and the source's11 flag labels do not replace any of these proofs.
