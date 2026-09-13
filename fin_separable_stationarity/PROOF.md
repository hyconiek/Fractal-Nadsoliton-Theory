# Stationary separability does not imply classicality

This is a new research goal after the completed Hartree-equivalence report.
The canonical interaction V=(|Omega><Omega|+S-2D)/2 is supplied. No
independently sourced FIN preparation, bath, clock or particle interpretation
is assumed to have been derived by the construction below.

## 1. Exact separable stationary completion of the declared strict marginal

Let C=I/12+W/20 and p=5/11. Set

    C'=I/12+(C-I/12)/p=I/12+11W/100.

The exact strict eigenvalue enclosures prove C'>I/125. Its diagonal is
uniform. A Paley Hadamard matrix H of order 12 is constructed explicitly
over the quadratic residues modulo 11; HH^T=12I. Its eleven nonconstant
columns z_r have six +1 and six -1 entries. Include both orientations.
In the resulting 22 cuts every label occurs in 11 positive halves and
every distinct pair in 5 positive halves.

For A_r={i:z_ri=1}, B_r its complement, define density matrices

    sigma_r=2P_A C' P_A,  tau_r=2P_B C' P_B,
    R_sep=(1/22) sum_r sigma_r tensor tau_r.

Each term is a product of positive trace-one states, proving separability.
The marginal diagonal remains 1/12, and each off-diagonal is multiplied
by 2*(5/22)=5/11. Both marginals are therefore exactly C.

The coordinate supports of sigma_r and tau_r are disjoint. Hence
DR_sep=R_sep D=0 and also |Omega><Omega|R_sep=0. Complementary cuts make
S R_sep S=R_sep. Therefore [V,R_sep]=[S/2,R_sep]=0 exactly.
For this real program, R_sep^(T_B)=R_sep as well, but this PPT fact is
not used as the proof of separability: the 22-product decomposition is.

Since every ordered unequal pair lies in 6 opposite oriented cuts,

    R_sep >= (12/11) lambda_min(C')² (I-D).

Thus its rank is exactly 132. Its coincidence probability Tr(D R_sep)
is zero, whereas Tr(D(C tensor C))=1/12. In particular its trace distance
from the product is at least 1/12; the numerically computed value is about
0.12141. This is a property of the constructed state, not an improved
universal correlation floor for all stationary completions.

For any uniform-diagonal program the same construction works. For a real
circulant kernel its guaranteed loading range is
gamma<=(5/11)/(12*(-lambda_min(W))). At strict this extends above 0.05.
This is a sufficient construction range, not a separability boundary.

## 2. Explicit preparation and its supplied resources

Choose one of the eleven unoriented cuts uniformly. On two independent
copies of C', perform the local projective measurements {P_A,P_B}.
Accept opposite outcomes and forget their orientation. The two success
Kraus operators are P_A tensor P_B and P_B tensor P_A. Both equal-half
probabilities are 1/2, so total success probability is exactly 1/2.
The normalized accepted state for that cut is the symmetric mixture of
sigma_r tensor tau_r and tau_r tensor sigma_r. Averaging cuts gives R_sep.

The two failure Kraus operators complete a trace-preserving instrument.
The construction uses local quantum operations and classical communication,
not entanglement. Conditioning on the declared success flag is explicit;
it is not a deterministic nonlinear channel on one unknown state.

The amplified C', shared cut choice, comparison/retention of outcomes and
fresh input copies are resources. C' already contains W. Forgetting the
preparation label gives a stationary two-body state; keeping that label
can reveal nonstationary conditional branches. No claim is made that the
whole system including every preparation record is an autonomous equilibrium.

The cut count 11 is minimal only for the declared isotropic balanced-cut
second moment: its matrix has rank 11, while each cut contributes rank one.
This does not prove minimal shared randomness or minimal physical preparation
cost over all protocols. Hadamard existence does not select the number 12
as a physical constant.

## 3. A general obstruction to zero one-sided discord at equilibrium

Let a bipartite Hamiltonian H have an invertible operator block
A=<u|H|v> on the second factor. Suppose its stationary state's first
marginal C has u,v as simple eigenvectors with unequal eigenvalues c_u,c_v.
Then

    [R,C tensor I] != 0.

Proof: if that commutator vanished, the simple u and v sectors of R would
be decoupled from every other first-factor sector. Write R_u=<u|R|u> and
R_v=<v|R|v>. The u,v block of [H,R]=0 then gives

    A R_v=R_u A.

Invertibility implies Tr R_v=Tr R_u, contradicting c_v!=c_u.

A classical-quantum state R=sum_i p_i |a_i><a_i| tensor rho_i necessarily
commutes with its own first marginal tensor I. Hence the preceding theorem
excludes zero discord on that side. This is a statement about the precise
classical-quantum form, not a claim that separability is equivalent to it.

For canonical V, the strict uniform u and alternating v give

    A=(|u><v|+|v><u|)/2-J/12,
    J=diag((-1)^i), J²=I,
    A²=(1/144)(I-P_u-P_v)+(25/144)(P_u+P_v).

Thus A is invertible with ||A^-1||=12. The real canonical C=I/12+gamma W
has simple u,v eigenvalues that differ for every gamma>0 in its PSD range.
EVERY stationary state with this marginal has nonzero one-sided discord.
If both marginals are C, the conclusion holds on both sides. In particular
the explicitly separable R_sep is not a classical-quantum state on either side.

This conclusion is canonical-interaction scoped. The exchange-shifted
family V+aS changes the block; at a=-(n-2)/(2n)=-5/12 it is singular.
The invertibility proof alone makes no verdict at that exceptional value,
and arbitrary pure-flow-invisible blocks are not covered automatically.
For any Hermitian perturbation ||E||<1/12, however, V+E retains the block's
invertibility, so the qualitative no-classical-quantum statement survives
provided the same marginal and stationarity premises hold.

## 4. Quantitative and definition-sensitive checks

For the explicit R_sep, coincidence exclusion gives a simple fixed-marginal
test. Every classical-quantum state with the SAME first marginal C must
contain its unique uniform eigenvector with probability c0. Its coincidence
probability is at least c0/12. Thus its trace distance from R_sep is at least
c0/12. This is NOT a bound on ordinary discord distance when the comparison
state's marginal is allowed to vary: at gamma=0, R_sep becomes a diagonal
classically correlated state, while C becomes degenerate. The apparent
fixed-marginal discontinuity must not be called a finite quantum resource
at zero loading.

A continuous, unrestricted trace-distance bound can instead be obtained
from the nonzero marginal commutator. If Q is classical-quantum and
delta=||R-Q||_1, then

    ||[R,C tensor I]||_1 <= 2(1+||C||) delta.

This follows by adding and subtracting Q and its marginal, using partial-
trace contraction and the commutator norm inequality. Therefore

    inf_(Q classical-quantum) D_tr(R,Q)
       >= ||[R,C tensor I]||_1/[4(1+||C||)].

For R_sep all entries are nonnegative. Its unequal-label diagonal entries
are 1/132. The matrix element of [R_sep,C tensor I] from |00> to |10>
has magnitude at least C_01/132, with no cancellation. Strict outward
bounds give C_01>1/50 and c0<1/6. Consequently

    inf_(Q classical-quantum) D_tr(R_sep,Q) > 1/15400.

This is a trace-distance-to-classical-quantum bound, not the numerical
value of entropic quantum discord or an entanglement measure.

For a universal but weaker bound on every canonical stationary R with
marginal C, pinch R in the distinct eigenspaces of C and call the result
R_tilde. The invertible-block equation and ||A^-1||=12 imply

    ||R-R_tilde||_1 >= (c0-c6)/(2*||V||*12).

If d_min is the smallest nonzero eigenvalue gap of C, the Frobenius norm
of [R,C tensor I] is at least d_min ||R-R_tilde||_F. Since the full matrix
dimension is 144, ||R-R_tilde||_F>=||R-R_tilde||_1/12. These bounds yield
a positive quantitative distance to the classical-quantum set using the
inequality above. Any numerical version must pay the strict spectral gap
with outward bounds rather than rounded eigenvalues.

## Outstanding source obligation

The canonical microscopic equilibrium can be separable but cannot have
zero discord for the stated strict marginal. Its conditional LOCC
preparation is explicit. An internal FIN source of the amplified coherent
program and the preparation mechanism remains absent; no entanglement,
clock, selector, legacy bridge or ToE closure is inferred from the construction.

## 5. A minimal cut frame is still not a unique correlation law

Permuting rows 1 and 2 of the Hadamard matrix preserves all exact cut counts,
positivity, separability, stationarity and the target marginal. Keep the
same labelled W and program C'. The two cut frames have different fourth
moments: the oriented count for labels (0,1) in A and (3,4) in B changes
from 1 to 2. Accordingly,

    Delta R_(03,14) = (2/11) C'_(01) C'_(34) != 0.

Both factors are positive for the declared strict program. An exhaustive
integer comparison of fourth-moment tensors over all 24 dihedral relabellings
shows that the frames are not equivalent under those kernel symmetries.
Their one- and two-point cut data and their minimal unoriented cut count
are nevertheless identical. Thus even this preparation architecture requires
a higher-order source/frame choice to predict the full joint state.

The separate canonical legacy cycle also admits the separable construction
at loading 1/1000, using amplified loading 11/5000. The bound |W_ij|<3
gives a positive input floor 1/12-33*(11/5000)>0. Its uniform and alternating
Fourier modes are simple and unequal, verified exactly over
Q(sqrt(2),sqrt(3)) after removing the common factor 4 ln(2). Hence the
qualitative invertible-block nonclassicality theorem applies there too.
The strict numerical discord-distance constants are NOT transferred to
legacy, whose entries are signed. No legacy physical role is inferred.
