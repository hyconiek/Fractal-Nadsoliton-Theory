# Conditional quantum lift: proof checkpoint

This is the new goal after completed ST8651/ST8652. Tensor replication,
pair interactions, and their normalization are explicit premises. Nothing
here derives an apparatus, dimensional clock, particle statistics or the
finite-speed dissipative controller in the original learning law.

## 1. What is actually lifted

The earlier fast-learning limit is rho_dot=(i/gamma)[T(rho),rho], with
T(rho)=Pi Re(rho). In the real symmetric off-diagonal orthonormal basis
X_ij=(|i><j|+|j><i|)/sqrt(2), T(rho)=sum X_ij Tr(X_ij rho).
Thus the supplied full source map has the two-body Hermitian lift

    V=sum X_ij tensor X_ij=(|Omega><Omega|+S-2D)/2,
    Tr_2[V(I tensor rho)]=T(rho),

where Omega=sum_i |ii>, S swaps the factors and D=sum_i |ii><ii|.
Its three orthogonal spectral projectors are

    P_Omega=|Omega><Omega|/n, eigenvalue (n-1)/2;
    P_+=(I+S)/2-D,         eigenvalue 1/2;
    P_-=I-P_Omega-P_+,     eigenvalue -1/2.

For n>=3 the first eigenspace is isolated and one-dimensional. Also
Tr_1 V=Tr_2 V=0. The antisymmetric space is contained in P_-, but is
not all of P_-: traceless diagonal vectors belong to it too.

The declared finite-copy Hamiltonian is

    H_N=-g/(N-1) sum_(a<b) V_ab,  g=1/gamma, N>=2.

Its ordinary unitary dynamics is linear and quantum-compatible on the
FULL N-copy input state. Nonlinear dependence on a common one-copy
preparation is not the earlier single-copy affine-map claim.

## 2. A controlled short-time mean-field limit

For initially factorized R_N(0)=C^tensorN let rho_N^(k) be its k-body
marginal and C_t solve C_dot=i g[T(C),C]. Write
f_k=||rho_N^(k)-C_t^tensor k||_1 and v=||V||=(n-1)/2 for n>=3.
The finite BBGKY hierarchy and unitary Duhamel formula give

    f_k(t) <= 3gv k(k-1)t/(N-1) + 2gv k integral_0^t f_(k+1)(s) ds.

The second term is absent at k=N. The three inhomogeneous units come
from one internal-pair commutator bound and two units from the finite
(N-k)/(N-1) versus infinite hierarchy coefficient. This uses only
||[V,R]||_1<=2v||R||_1 and contraction of partial trace for Hermitian
differences (including i times the anti-Hermitian commutators).

Iterating all the way to k=N, with f_k(0)=0, proves

    f_1(t) <= [3/(2(N-1))] sum_(r=1)^(N-1) r q^(r+1)
            <= 3 q²/[2(N-1)(1-q)²],  q=2gv t<1.

This is an explicit finite-N estimate on the stated short-time interval.
It is not promoted to a uniform-in-time or long-time stability bound.
No claim is made that the independent-copy initial condition is sourced
by FIN rather than supplied.

## 3. Full-rank product stationarity is much more restrictive

**Collective-commutant theorem.** For n>=3 and Hermitian h,

    [V,h tensor I+I tensor h]=0  iff  h is scalar.

Proof: the isolated Omega eigenspace forces
(h tensor I+I tensor h)Omega=vec(h+h^T)=2vec(Re h) to be proportional
to Omega. Thus h=cI+iA with A real antisymmetric. For every traceless
diagonal M, vec(M) is a -1/2 eigenvector of V. Its image under the
collective h has an off-diagonal symmetric component i vec(AM-MA),
which belongs to the +1/2 eigenspace. It must vanish. Choosing M with
distinct diagonal entries forces every off-diagonal A_ij to vanish.
The converse is immediate. The n=2 exception is real: h proportional
to X, besides I, commutes collectively with V=X tensor X/2.

**Product no-go.** A full-rank density C satisfies
[H_N,C^tensorN]=0 for some finite N>=2 iff C=I/n.

For N=2, functional calculus turns this commutation into
[V,log C tensor I+I tensor log C]=0. The preceding theorem finishes
the proof. For general N take log(C^tensorN)=sum_a (log C)_a.
Every pair commutator has zero one-body partial traces, since both
partial traces of V vanish. Taking the partial trace of the sum of
pair commutators onto any chosen pair therefore retains only that
pair, with a nonzero scalar factor. The same theorem applies.

This is not a blanket no-go for correlated stationary states. It also
does not classify rank-deficient product inputs. The actual strict
C=I/12+gamma W is full rank at gamma=1/20, so its exclusion is paid.

If C is a Hartree equilibrium, one-body stationarity is initially true,
but the exact first nonzero reduced acceleration is

    rho_N^(1)''(0) = -g²/(N-1) Tr_2[V,[V,C tensor C]].

The external-three-body terms vanish because [T(C),C]=0. In particular
the microscopic correction is order 1/(N-1) at this order; extrapolating
the Taylor expression to times growing with N would require a new proof.

## 4. Exact two-copy strict marginal trajectory

For real symmetric C with diag C=I/n and diag C²=(Tr C²)I/n,
which includes every real symmetric circulant strict reference,
the two-copy evolution at phase tau=g t has the exact marginal

    C(tau)=C+A(tau)(C-I/n)+B(tau)(C²-(Tr C²)I/n),
    A(tau)=2(cos tau-1)/n,
    B(tau)=2[cos((n/2-1)tau)-cos tau]/n.

Proof: S commutes with C tensor C. The remaining unitary factor is
I+(alpha-1)D+(beta-alpha)P_Omega, alpha=e^(-i tau),
beta=e^(i(n/2-1)tau). Expand it on both sides of C tensor C and take
the partial trace. The constant diagonal and constant square diagonal
give the displayed expression; no numerical spectral fit is used.

For n=12 the frequencies are 1 and 5. At tau=pi,
C(pi)=(2/3)C+(1/3)I/12, while at 2pi it returns to C exactly.
Uniform vertex populations do not reveal this change; the local spectrum
and entropy do. The product Hartree fixed point is not an exact finite
microscopic fixed point, and the trajectory is not dissipative relaxation.

Global entropy is conserved. Exchange symmetry and subadditivity imply
S(C(t))>=S(C) for a product input at any time. More precisely the total
correlation D(R_N(t)||rho_N^(1)(t)^tensorN) equals
N[S(rho_N^(1)(t))-S(C)]. This does not imply monotonic entropy in time;
the explicit revival refutes monotonicity.

## 5. Same marginal, stationary correlated completion

Let A_-=(I-S)/2, n>=3, and assume lambda_min(C)>=1/[2(n-1)]. Then

    R_A(C)=2/(n-2) A_- [C tensor I+I tensor C-I/(n-1)] A_-

is PSD, trace one, has BOTH marginals C, and is exactly stationary for V.
In an eigenbasis of C it is the mixture of normalized Slater projectors
with weights 2/(n-2)[c_i+c_j-1/(n-1)], i<j. These weights are nonnegative,
sum to one, and their incident sums are 2c_i. This proves all marginal
and positivity claims. V is -I/2 on the entire antisymmetric subspace.
The strict gamma=1/20 density satisfies the sufficient floor by the
inherited exact spectral enclosure, not just floating eigenvalues.

The state is necessarily entangled in this construction. Indeed
<Omega/sqrt(n)|R_A(C)^(T_2)|Omega/sqrt(n)>=Tr(S R_A)/n=-1/n.
Thus its partial transpose has a negative eigenvalue and negativity
at least 1/n. For the numerical strict instance it is about 0.09553.
This is a mathematical entanglement witness, not evidence of fermions
or a laboratory preparation. The antisymmetric code can be used for
distinguishable tensor factors; no particle-statistics identification follows.

This construction is target dependent. It is not a universal broadcasting
channel taking one unknown C to a two-copy state with marginals C. Such
an inference would silently erase the source/preparation problem.

## 6. Falsifying the apparent completion

R_A(C) is a stationary MARGINAL completion, not a dynamical FIN realization.
Every state in its support lies in the same flat energy band. Any common
unitary pulse U tensor U keeps it in that band, so R_A(U C U*) also stays
stationary even when the corresponding Hartree equation would move C.
The construction stores the prescribed kernel but does not generate its
nontrivial propagator. Inserting a separate local W Hamiltonian would be
additional source data, not an internal consequence of this flat band.

For two copies V has only three energies, so fixed linear preparation and
readout can carry at most the positive Bohr frequencies 1, n/2-1, n/2,
times the global rate. Full strict propagation has more frequencies.
P512/O216 already proved a strict Laplacian frequency ratio transcendental
for the declared decimal benchmark. This is not rediscovered here.
Since every finite N Hamiltonian above has rational entries up to one
overall scale, every nonzero Bohr-frequency ratio is algebraic. Applying
the OLD P512 result therefore excludes exact reproduction of the complete
frozen strict unitary time law by this NEW finite-N class with fixed linear
encoding/readout. Approximation, an infinite limit, extra independently
sourced couplings, time-dependent controls, or changed exact parameters
are not excluded. Arithmetic nonidentity alone is not a robust experimental
error bound, and it is not a universal no-physics theorem.

## Outstanding work

Exact structural tests, independent finite-N hierarchy checks, careful
comparison with primary mean-field/no-broadcasting literature, and the
source significance audit remain to be completed. No final breakthrough
or physical completion is claimed merely from this checkpoint.
