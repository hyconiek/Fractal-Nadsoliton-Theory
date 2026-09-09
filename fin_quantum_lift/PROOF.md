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

There is a stronger version that does not assume identical copies or
all-to-all geometry. For n>=3 and Hermitian h,k,

    [V,h tensor I+I tensor k]=0  iff  h and k are both scalar.

The isolated Omega first forces h+k^T=cI. On traceless diagonal vectors,
the imaginary antisymmetric part of h would move the -1/2 band into the
+1/2 band, so that part vanishes. Write h real symmetric and k=cI-h.
On every off-diagonal symmetric matrix M in the +1/2 band, the remaining
commutator [h,M] is real antisymmetric in the -1/2 band, hence must vanish.
Commutation with all E_ij+E_ji forces h scalar for n>=3: use a third index
to kill each off-diagonal entry, then compare diagonal entries.

Consequently, on ANY finite graph with nonzero canonical V couplings on
its edges, and with arbitrary additional one-site Hamiltonians, a stationary
FULL-RANK product state must be maximally mixed at every non-isolated site.
Proof: take the logarithm of the product state. Each edge commutator has
zero one-body partial traces. Projecting the total commutator onto its
two-site, zero-partial-trace component isolates each nonzero edge; the
one-site Hamiltonians cannot cancel it. Apply the two-field theorem above.
For a connected graph the whole product is maximally mixed.

Full rank is essential. The product |+_ij> tensor |-_ij>, with balanced
real superpositions on two labels, lies in the -1/2 energy band and is a
rank-one product stationary state. With zero or compatible local fields,
bipartite graphs admit corresponding alternating product eigenstates.
This does not realize a full-rank strict
marginal. For the canonical strict C_gamma even the boundary gamma_star
has nonzero reduced acceleration: its seven distinct eigenvalues cannot
all solve the nonconstant quadratic acceleration polynomial in section 4.

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

The arithmetic obstruction is kernel-sensitive. In the separately supplied
canonical legacy cycle, all weights equal the common factor 4 ln(2)
times algebraic numbers: its phases pi*d/4+pi/6 are rational multiples
of pi and its attenuation denominators are rational. Its nonzero Fourier
gap ratios are therefore algebraic. P512's strict obstruction cannot be
transferred to that reference. This does not establish actual realization
of legacy by V_N, bridge completion, or any physical-role transfer.
The product-stationarity and correlated-completion arguments DO apply
separately to both references. For legacy loading 1/1000, the coarse bound
|W_ij|<3 gives C>= (1/12-33/1000)I > I/22. Signed legacy weights are not
interpreted as positive classical rates in this calculation.

## 7. A different, operationally valid sample-programmed construction

Supply a probe state sigma and N fresh independent program states C,
independent of the unknown probe. At each collision apply
U_delta=exp(i g delta V), trace out that program, and use delta=t/N.
The resulting channel Phi_delta^N is exactly CPTP. Its limiting Hamiltonian
is g h, h=T(C). For C=I/n+gamma W and g=1/gamma the target density-channel
is conjugation by exp(iWt), equal to the frozen Laplacian channel up to
a scalar phase. This is not the original state-dependent single-copy law.

This architecture is a variant of known sample-based Hamiltonian simulation
and density-matrix exponentiation, not a new general quantum algorithm.
In particular a SWAP interaction is another option for canonical C, since
C and T(C)=C-I/n differ only by a scalar. The canonical lift studied here
is chosen because it reproduces the FULL supplied map T, not because all
other quantum implementations have been excluded.

Here is a finite, reference-assisted error bound for this specific lift.
Let V0=V-h tensor I and

    B_C=Tr_2[V0²(I tensor C)]
       =[I+(n-2)C^T]/4-h² >=0.

The equality follows from V²=I/4+(n-2)|Omega><Omega|/4. Let
E0(sigma)=Tr_2[V0(sigma tensor C)V0], a CP map with E0*(I)=B_C, and
D_C(sigma)=E0(sigma)-{B_C,sigma}/2. The exact second-order comparison is

    Phi_delta-U_h(delta)=g² delta² D_C + O(delta³),
    ||D_C||_diamond <= 2||B_C||.

Unitary-channel Taylor remainders are bounded without an exponential:
their third derivatives have diamond norm at most (2||Hamiltonian||)^3.
Thus the ONE-step difference is bounded by

    2g² delta² ||B_C|| + (4/3)g³|delta|³(v³+||h||³).

Telescoping the N channels, whose diamond norms are one, gives the
rigorous bound

    ||Phi_(t/N)^N-U_h(t)||_diamond
      <= 2g² t² ||B_C||/N
         +(4/3)g³|t|³(v³+||h||³)/N²,

capped by two. This norm covers inputs entangled with a reference.
For strict n=12, gamma=1/20, ||W||<3 gives ||h||<3/20,
||C||<7/30 and ||B_C||<=5/6. The bound becomes

    (2000/3)t²/N + (5324108/3)|t|³/N².

At t=1, N=1,000,000 it is exactly bounded by
501331027/750000000000 < 0.0007. This is a mathematical resource/error
certificate, not an implemented laboratory system or an SI time estimate.
Fresh samples and their encoded state, interaction timing and the overall
rate are real supplied resources.

## 8. An actual finite-copy leakage obstruction, not just an upper bound

Let u be the uniform Perron vector and P=|u><u|. For ANY program sigma
with T(sigma)=gamma W, c0=<u|sigma|u>=1/n+gamma s is fixed, where s is
the strict row sum. A direct calculation of <u|U_tau|u> on the program
space yields A_tau I+B_tau P. Consequently the exact survival loss of
the probe P in ONE collision is

    L(tau,c0)=(1-c0)(1-4/n²) sin²(tau/2)
              +4(n-1)c0/n² sin²((n-2)tau/4).

It is independent of diagonal or imaginary program data hidden from T.
The coefficient of tau² is

    f0=(n²-4)/(4n²)+(n-2)(n-4)c0/(4n).

For n=12 and the strict feasible programs c0 lies strictly between zero
and one, so this coefficient is positive. The n=2 exception has zero
leakage and is retained in the tests.

For each FIXED program and gamma, analyticity near delta=0 gives
log Phi_delta=delta L_h+g²delta²D_C+O(delta³). Expanding its Nth power
at delta=t/N and using [P,h]=0 gives

    1-Tr[P Phi_(t/N)^N(P)] = g² t² f0/N + O(N^-2).

The coefficient is the same for all programs with the declared T(sigma),
although their other finite-copy noise can differ. Hence the trace-distance
error, and half the diamond-norm error, are at least this positive survival
loss. This is a 1/N lower witness for this processor class, not a lower
bound for every quantum simulation algorithm or correlated program supply.

For real C the leading D_C is a unital, self-adjoint Lindblad dissipator
with covariance K_ab=Tr(C X_a X_b)-Tr(C X_a)Tr(C X_b). If C>=c_min I,
then K>=c_min I on the real off-diagonal observables. The reference sum
sum D[X_a] has eigenvalues -n/2 on traceless diagonals and imaginary
antisymmetric matrices, and -(n-2)/2 on real symmetric off-diagonal
matrices. This supplies a positive leading dissipative gap for n>=3.
It does NOT assert that the exact finite-step channel is unital or has
the same stationary states. Complex programs with the same T can already
change the leading noise: D_sigma(I)=(n-2)(sigma-sigma^T)/4.

## 9. Conditional program-loading optimization and its limitation

For any density sigma satisfying T(sigma)=gamma W, dihedral twirling
gives C_gamma=I/n+gamma W. This is an actual mixture of real permutation
unitaries, not an antiunitary transpose operation. Therefore feasibility
of ANY such program implies

    0<gamma<=gamma_star=1/(n ell),  ell=-lambda_min(W),

and C_gamma witnesses sufficiency. Convexity of the largest eigenvalue
and covariance of B_sigma under the same group show that C_gamma
minimizes ||B_sigma|| at each fixed gamma. A minimizer is identified;
uniqueness is not claimed.

For n=12, the eigenvalues of B_C/gamma² are

    11/(24 gamma²) + (5/(2 gamma))lambda - lambda².

On the feasible range this is increasing in lambda, since
5/(2gamma)-2s >=30ell-2s >=8ell>0, using s<=11ell.
The maximum is thus the Perron value and decreases strictly as gamma
increases. The optimum of this variance bound is at gamma_star, with

    min ||B_sigma||/gamma² = 66ell²+30ell s-s².

The unavoidable Perron coefficient likewise satisfies

    f0/gamma²=55/(144gamma²)+5s/(3gamma),
    min f0/gamma²=55ell²+20ell s.

Outward strict Fourier enclosures certify the scalars in these formulas.
Numerically gamma_star is about 0.1222120805, the optimized variance
coefficient about 61.8939393, and the minimum Perron-leakage coefficient
about 48.2148582. These are engineering bounds for the declared processor,
not a derivation of an absolute physical normalization or a source for W.

In particular ||B|| alone is not a universal operational cost: for an
environment-only interaction I tensor Z, B can be positive while the
probe channel is exactly the identity. The variance optimization above
is not promoted to optimality over all interactions or algorithms.
The separate survival-loss calculation is what certifies real error here.

## Outstanding work

Structural tests and finite-N hierarchy checks are now implemented.
Primary sample-simulation/no-programming comparisons and the source
significance audit are being reconciled with the exact bounds. No physical
completion is claimed merely from a conditional simulation construction.
