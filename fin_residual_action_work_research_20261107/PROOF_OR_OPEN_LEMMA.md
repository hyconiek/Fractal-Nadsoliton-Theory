# AW - residual work closes the conditional three-channel action law

Date 2026-11-07. Parent 400c26f5d7cfe4703797d900d636ba35b5e9ff14.
Self-audited conditional analysis using the supplied RP model and the
recorded CM, QP, CN, ML, CE, TP2, MA and RW estimates.

## 1. Target and result

Keep RW9's exact residual works D_I,D_J and its zero-data sources RW8.
For each fixed S, on T=S epsilon^-2, we prove

    sup_(0<=t<=T)(|D_I(t)|+|D_J(t)|)
                          <=C_S epsilon^3[1+log(1/epsilon)]. (AW1)

Thus the unpaid work is o(epsilon^2). With RW10 the actual corrected
canonical actions satisfy

    I(t)=I(0)-integral_0^t[k_a I^2+k_x I J]ds+r_I(t),
    J(t)=J(0)-integral_0^t[k_b J^2+k_x I J]ds+r_J(t),
    sup_(t<=T)(|r_I|+|r_J|)
                          <=C_S epsilon^3[1+log(1/epsilon)]. (AW2)

The coefficients and physical reconstruction retain CD and MA's exact
normalizations. The proof controls the work pairing through adjoint
propagation. Every constant is uniform on the original compact mass band,
with dependence on fixed S and finite signed model coefficients allowed.

## 2. Uniform spatial bound for the modulated adjoint primitive

Let A be H or L on its odd/even continuous subspace, with far mass m_A.
Let J_B be the continuous projection of one of Fa,Fb,G used in RW9.
Its Darboux transform is Schwartz. This follows from the explicit adjoint
formula and smooth massive Fourier multipliers, applied to the exponentially
decaying profiles and their finite-rank projections.

For a carrier nu define the finite primitives

    C_nu(t)=integral_0^t exp(i nu s) cos(s sqrt(A))J_B ds,
    S_nu(t)=integral_0^t exp(i nu s) sin(s sqrt(A))/sqrt(A) J_B ds.

The work weights have carriers +/-2omega for PZ, +/-2lambda for RC,
and +/-lambda+/-omega for PC and RZ. Each carrier is separated uniformly
from both threshold frequencies +/-m_A. The self-channel gap satisfies
2omega-2>=sqrt(5)-2>1/5, and QH4 pays the mixed gaps. The continuous
projection removes all discrete eigenmodes before forming these primitives.

For spatial derivatives j=0,1,2,

    ||partial_x^j C_nu(t)||Linfinity
       +||partial_x^j S_nu(t)||Linfinity<=C log(2+t). (AW3)

To prove AW3, use the exact Darboux eigenfunction e_k=Wcal exp(ikx).
Its j-th spatial derivative is bounded by C_j<k>^j uniformly in x.
The source transform is Schwartz, so the spectral amplitude, including
these polynomial factors and optional (m_A^2+k^2)^-1/2, is integrable.
Each time integration supplies a factor bounded by

    min(t,2/|nu +/- sqrt(m_A^2+k^2)|).

If a denominator has a real zero, its roots are simple and uniformly
separated from zero momentum. In a fixed small neighborhood of either
root the change to that denominator has a bounded Jacobian. Integrating
min(t,2/|y|) over that neighborhood costs at most C log(2+t). Away from
the roots the bound is uniform by the gap and Schwartz decay. Closed
carriers have no roots and satisfy the same bound. This proves AW3 for
all x without assigning a global energy norm to a stationary wave.

The free multiplier Omega_A=sqrt(-partial_x^2+m_A^2) also satisfies

    ||Omega_A S_nu(t)||Linfinity<=C log(2+t).          (AW4)

Indeed Omega_A S=Omega_A^-1(-S''+m_A^2 S), and the inverse multiplier
has an L1 kernel of norm 1/m_A by CE's heat representation. Its order
relative to the full-potential propagator is preserved in this argument.

## 3. The backward work tests have size epsilon^2 log(2+T)

Fix any terminal time t<=T. Each modal work weight in RW9 is a finite
sum of A(s)exp(-i nu s). MA14 and the algebraic amplitude reconstruction give

    ||A||Linfinity<=C epsilon^2,
    ||A'||L1(0,T)<=sqrt(T)||A'||L2<=C_S epsilon^2.     (AW5)

Fubini transforms a zero-data response work integral into source pairings
with the backward tests

    T_C(s)=integral_s^t weight(r)cos((r-s)sqrt(A))J_B dr,
    T_S(s)=integral_s^t weight(r)sin((r-s)sqrt(A))/sqrt(A)J_B dr.

Integration by parts using the finite primitives AW3, at the actual
endpoints s and t, gives

    ||T_C(s)||Linfinity+||T_S(s)||Linfinity
            +||Omega_A T_S(s)||Linfinity
                             <=C_S epsilon^2 log(2+T). (AW6)

For example integrate A(s+tau) against the primitive of
exp(-i nu tau)cos(tau sqrt(A))J_B. The boundary costs ||A||infinity
and the integral costs ||A'||L1. The estimate includes every endpoint
term and uses neither an infinite-time steady-state substitution nor
an independent slow-envelope hypothesis. It is uniform in s and t.

For a physical phase source (F_1,F_2) whose components are integrable,
AW6 bounds its work by epsilon^2 log(2+T) times their time-integrated
L1 norms. More generally for a structured source (A_0,Omega_A B_0+C_0)
it bounds work by the same factor times

    integral_0^T (||A_0||L1+||B_0||L1+||C_0||L1) ds. (AW7)

The momentum pairing transfers Omega_A onto T_S, as in AW4. Sources are
projected in physical phase. Since the adjoint tests are already continuous,
that projection drops out of the pairing before transferring Omega_A.
No commutation of Omega_A with the physical projection is assumed.

## 4. Pay every term of the exact residual source

Use physical eta=cg+v and psi=zf+u reconstructed from the canonical chart.
Split RW8 into the following exact source classes in each output.

The bound reconstruction differences are

    -6[(c^2-C^2)Fb+(z^2-Z^2)Fa] in the odd field,
    -12(cz-CZ)G in the even field.

CM gives |c-C|+|z-Z|<=C epsilon^2. These sources have every fixed positive
weighted L2 norm O(epsilon^3). Their time-integrated L1 norms are O_S(epsilon).

The mixed quadratic mode/continuum sources are

    -12K(cg v+zf u) in the odd field,
    -12K(cg u+zf v) in the even field.

Their L1 norms are at most C epsilon E(s), because f,g localize and E is
MA2's physical local continuum norm. Hence their time-integrated L1 norms
are O_S(epsilon). The full physical cubic and quintic forces have global
L1 size O(epsilon^3), by two L2 factors and the remaining global Linfinity
factors. Their time-integrated L1 norms are also O_S(epsilon).

The remaining quadratic continuum-only forces are

    G_phi=-6K(v^2+u^2), G_psi=-12Kvu.                  (AW8)

Their unweighted size alone is insufficient. Apply the exact CN
continuum-input normal form, extended to the odd output with QP's two
channels (phi;phi,phi) and (phi;psi,psi), each with coefficient -3/2.
The even channel is (psi;phi,psi) with coefficient -3. For each channel
retain CN6's actual input forcing F_V+F_G+F_L+F_C, and the conjugated
physical phase projection. The algebra of CN8--CN9 applies term by term
to these three channels, since their symbols obey the same QP bounds.
Here the auxiliary response contains G_i alone: set N_X=0 in CN9,
while retaining the full actual F_C in the input equations. The direct
physical cubic force has already been counted separately above. For an
output i and channel (i;j,k) with coefficient a, the explicit ledger is

    Q_i^sigma = sum_channels sigma a K sum_(tau,rho) B_b(h_j^tau,h_k^rho),
    R_2^sigma = sigma i Omega_i^-1 V (Q_i)_pos
      +sum_channels,branches [sigma i a [Omega_i^-1,K](h_j h_k)
        -i a [Omega_i,K]B_b(h_j,h_k)
        -sigma a K(B_b(F_V,j,h_k)+B_b(h_j,F_V,k))],
    R_3loc^sigma = -sum_channels,branches sigma a K
                       (B_b(F_L,j,h_k)+B_b(h_j,F_L,k)),
    R_cascade^sigma = -sum_channels,branches sigma a K
         (B_b(F_G,j+F_C,j,h_k)+B_b(h_j,F_G,k+F_C,k)).

Input branch indices in the last three lines are those in the first line.
The potential term is counted once after summing Q_i. Apply the conjugated
physical projection only after this identity. This displays both odd
channels explicitly, without identifying a projected free half-wave with
a componentwise projected position or velocity.
The exact projected response to AW8 equals its quadratic correction Q_c
plus a zero-data history of the following residuals.

(a) Localized quadratic residuals, including potential and commutator
terms, have positive-weight phase norm at most C epsilon E(s). This is
ML1's kernel and mixed-norm proof for each of the finite QP channels.
All mass multipliers and symbol derivatives used there are uniform on
these same triples. Its time-integrated component L1 norm is O_S(epsilon).

(b) Terms with F_L have positive-weight phase norm O(epsilon^3), by
QP's one-localized-input estimate. The L_i terms include the retained
bound-profile forcing and finite-rank quadratic subtraction. Their
integrated component L1 norms are O_S(epsilon).

(c) Terms with F_G or F_C are finite sums of K B_b(F_G+F_C,h).
CE1 gives their half-wave L1 norms at most C(epsilon^3+epsilon^4).
After physical phase conversion they have the form (A_0,Omega_A B_0)
with ||A_0||L1+||B_0||L1<=C epsilon^3. These bounds use the global H1
sizes of h,F_G,F_C from CN6. Their time-integrated preimage L1 norms
are O_S(epsilon), and AW7 pays their work.

The quadratic correction has local phase size at most C epsilon E(s)
by the mixed estimate used in ML. Its initial value is exactly zero,
since the physical continuous phase starts at zero. This is CN's
continuum-input correction. Its work is bounded directly by

    C epsilon^2 integral_0^T epsilon E(s)ds<=C_S epsilon^3.

All positive-weight phase bounds in (a),(b) imply component L1 bounds
by weighted Cauchy-Schwarz. Physical phase conversion and finite-rank
projection preserve the recorded weighted norms. The two odd channels
have exactly the same potential, commutator and input-forcing ledger
as CN's even channel, with the corresponding mass multipliers.

Thus every residual history in RW8 has time-integrated L1 source or
preimage norm O_S(epsilon), after the exact normal form of AW8. Applying
AW6--AW7 and adding the correction work proves AW1. The Fubini identities
hold first for smooth finite-time data. Approximation in the paid energy,
weighted and preimage norms extends the estimates to the actual solution
class. All sources have zero-data retarded responses as in RW8.

## 5. Compare the actual actions with the dissipative ODE

Let (i,j) solve

    i'=-k_a i^2-k_x i j, j'=-k_b j^2-k_x i j,
    (i(0),j(0))=(I(0),J(0)).                          (AW9)

The reference solution stays nonnegative and O(epsilon^2). The actual
actions are nonnegative and O(epsilon^2) by their definition and global
energy. On that box the polynomial vector field has Lipschitz constant
C epsilon^2. Subtract the integral equations AW2 and AW9 and apply
Gronwall on T=S epsilon^-2 to obtain

    sup_(t<=T)(|I(t)-i(t)|+|J(t)-j(t)|)
                           <=C_S epsilon^3[1+log(1/epsilon)]. (AW10)

This is a conditional critical-window averaged action law for the supplied
prepared PDE. It is an approximation with a quantified error, which does
not assert pointwise monotonicity of the exact actions.

For the preparation retained in MA, I(0)=omega epsilon^2/2+O(epsilon^4)
and J(0)=O(epsilon^4). Since j(t)<=j(0), the mixed term changes i by at
most O_S(epsilon^4) on the critical window. Therefore

    I(t)=(omega epsilon^2/2)
                    /[1+(k_a omega/2)epsilon^2 t]
               +O_S(epsilon^3[1+log(1/epsilon)]),
    J(t)<=C_S epsilon^3[1+log(1/epsilon)].             (AW11)

MA's O(epsilon^3) physical action reconstruction error is absorbed by
this bound. AW11 describes actions. A pointwise field waveform and
all-time scattering would require separate estimates.

## 6. Review status and next research step

AW1 completes RW's residual work estimate and AW10 gives the conditional
three-channel action approximation. The proof depends on the predecessor
global existence/smallness, QP symbol bounds, exact canonical and continuum
maps, ML mixed kernels and TP2 critical radiation estimate.

The next step is a focused proof audit of this dependency chain, especially
the odd-channel extension of CN/ML in section 4, physical phase projections,
and the uniform-in-space logarithmic primitive in section 2. Independently
reproduce the residual ledger and action reconstruction before treating this
self-audited theorem as reviewed. Natural source selection, units, empirical
validation and physical interpretation remain distinct research obligations.
