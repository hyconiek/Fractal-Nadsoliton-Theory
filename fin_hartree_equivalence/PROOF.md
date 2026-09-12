# Identifying a microscopic interaction from projective self-dynamics

New goal after the completed quantum-source report. No PDF. The statements
below concern a declared instantaneous Hermitian two-body Hamiltonian and
the fixed normal-ordering basis, not all possible FIN microscopic theories.

Let S swap two n-dimensional factors, P_s=(I+S)/2, P_a=(I-S)/2,
D=sum_i |ii><ii| and Omega=sum_i |ii>. The earlier canonical lift is

    V=(|Omega><Omega|+S-2D)/2,
    L_V(rho)=Tr_2[V(I tensor rho)]=T(rho)=Pi Re rho.

Projective dynamics means F_V(rho)=i[L_V(rho),rho]. It does not fix the
phase of a state vector or identify a laboratory clock. Reciprocity below
means S W S=W for the compared microscopic interaction W.

## 1. Three genuinely different identification premises

**Full source.** If L_W(P)=T(P) for every pure projector P, then W=V.
Pure projectors span the Hermitian matrices, so the linear maps are equal
on their full domain; the partial-trace pairing identifies W uniquely.

**All mixed projector flows.** If F_W(rho)=F_V(rho) for every density
rho, then for reciprocal W,

    W=V+a S+b I,  a,b real.                                  (1)

For a general Hermitian W without reciprocity the corresponding expression
is V+aS+I tensor B with Hermitian B: the extra part is an environment-only
field, not a new two-body interaction. The second subsystem need not have
the same mean-field law in that nonreciprocal case.

Proof of (1): a real-linear Hermiticity-preserving map L satisfying
[L(rho),rho]=0 for all densities satisfies the same polynomial identity
for all Hermitian matrices. First extend on the trace-one affine space,
then use homogeneity and continuity at trace zero. Fix a basis. For every
simple-spectrum diagonal D, L(D) is diagonal, hence the whole diagonal
space is preserved. Polarizing [L(D+tX),D+tX]=0 with off-diagonal
Hermitian matrix units X shows that

    (L(D))_i-(L(D))_j=a_ij(d_i-d_j).

Comparison on triples of indices makes all a_ij equal to one real a
(there is only one pair for n=2). The same comparison for both real and
imaginary matrix units shows that L-a id maps every Hermitian matrix
to a diagonal one. For a matrix with connected nonzero off-diagonal
support, that diagonal matrix must be scalar. Linearity extends this to
all matrices. Thus L(X)=aX+f(X)I for a real linear functional f.
Reciprocity makes L self-adjoint for the Hilbert--Schmidt pairing, which
forces f(X)=b Tr X. The two-body representatives are S and I.
The standard commuting-map principle is known mathematics; this finite
proof and its FIN application do not establish global priority.

**Only pure projector flows.** If F_W(P)=F_V(P) only for all pure P,
the exact reciprocal classification is much larger:

    W=V+c P_s+B_a,   B_a=P_a B_a P_a Hermitian.               (2)

Proof: write Z=W-V. Reciprocity splits Z into symmetric and antisymmetric
blocks. E(psi)=<psi psi|Z|psi psi> on the unit sphere has derivative
4 Re<dpsi|L_Z(P)|psi>. Vanishing projective flow makes this derivative
zero on every tangent direction, so E is constant c. Homogeneity and
complex polynomial polarization imply P_s Z P_s=cP_s. Conversely an
antisymmetric block annihilates every psi tensor psi, while cP_s adds
only a common scalar to the pure-state Hamiltonian action.

At n=12, dim(antisymmetric)=66. The invisible real dimensions are therefore
0 for the full source, 2 for all mixed projector flows, and 1+66²=4357
for pure projector flows. Removing a global energy shift leaves 4356 in
the last case. Requiring the exact pure-vector phase removes the c freedom,
but not the 4356-dimensional antisymmetric block.

This does not mean the original full controller source has 4357 free
parameters. Its full T is a stronger premise. Nor may the supplied mixed
extension be treated as a consequence of pure-state dynamics alone.

## 2. The mixed-equivalent exchange family has one real exception

Put V_a=V+aS. Its bands are

    Omega: (n-1)/2+a; off-diagonal symmetric: 1/2+a;
    traceless diagonal: -1/2+a; antisymmetric: -1/2-a.

For n>=3 the solutions of [V_a,h tensor I+I tensor k]=0 are scalar h,k,
except at a=-1/2. At that value they are exactly diagonal h,k with
h+k=cI. The proof separates the four bands. Generic a has no band
coincidences. At a=0 the earlier canonical proof applies. At a=-n/4,
Omega meets the antisymmetric band, but testing traceless diagonal vectors
forces the imaginary parts to vanish and then forces h,k scalar. At
a=-1/2 the off-diagonal symmetric and antisymmetric bands coincide,
leaving the diagonal staggered fields. There are no other coincidences.

Thus a=-1/2 gives full-rank nonuniform product equilibria:

    C=diag(p_i),  E=C^-1/Tr(C^-1).

For example p=(1/2,1/3,1/6), E=diag(2,3,6)/11. This is an exact
counterexample to extending the OLD canonical 'only maximally mixed'
claim to every mixed-flow-equivalent lift. Both T(C) and T(E) are zero,
so it does not realize strict or another nonzero learned kernel.

For a connected graph of this exceptional interaction, a nonuniform
full-rank product must alternate C and its normalized inverse. An odd
cycle forces uniformity. On a bipartite graph the alternating product is
stationary when the local fields commute with the assigned factors.
Other exchange parameters retain the canonical full-rank obstruction.

## 3. The much larger pure-equivalent family still has an exact obstruction

**Rank-free identical-pair theorem.** For n>=3, every reciprocal W in (2)
and EVERY density C, including rank-deficient C,

    [W,C tensor C]=0  iff C=I/n.                            (3)

First prove that no nonzero proper subspace M of C^n has Sym² M invariant
under V_s=P_s V P_s. Suppose otherwise and choose psi in M. Subtracting
the known (1/2)psi tensor psi term in V_s(psi tensor psi) shows that the
diagonal matrix

    diag(q/2-psi_i²),  q=sum_i psi_i²,

has range contained in M. If any coordinate e_j belongs to M, taking
psi=e_j makes that diagonal matrix invertible, forcing M=C^n. Otherwise
no coordinate belongs to M, so all its diagonal coefficients must vanish
for every psi in M. Summing psi_i²=q/2 over i gives (n-2)q=0. Since n>=3,
q=0 and every psi_i=0, another contradiction. The argument is complex,
not restricted to real subspaces or a nonzero value of q.

If C has proper support M, stationarity would make M tensor M invariant
under W. Reciprocity then makes Sym² M invariant under its fixed
symmetric block V_s+cI, contradicting the lemma. Thus C is full rank.
Taking log(C tensor C) yields the collective field h tensor I+I tensor h,
h=log C. On the symmetric sector, its commutation with V_s forces h
scalar by the earlier isolated-Omega/diagonal/off-diagonal argument.
This proves (3). The n=2 exception is not removed.

The same conclusion holds for identical products C^tensorN in the
permutation-invariant all-to-all pair model, N>=2. For proper support M,
invariance of Sym^N M implies invariance of Sym² M under V_s: apply the
pair Hamiltonian to psi^tensorN and separate the components with one and
two factors in M-perp. The one-outside component is a symmetric product
of a vector with psi^tensor(N-2), and multiplication by a nonzero
psi^tensor(N-2) is injective. Hence both outside components of V_s on
Sym² M vanish. Full rank then reduces the logarithmic commutator to the
symmetric-sector two-body collective commutator by evaluating coherent
products. This extension assumes the stated permutation-invariant pair
sum; it is not silently applied to arbitrary driven graphs.

## 4. Unequal full-rank partners: a necessary matching condition

Let W be any reciprocal member of (2), with no extra Hamiltonian outside
that declared equivalence class, and suppose [W,C tensor E]=0 for full-rank
densities C,E. Then

    E=C^-1/Tr(C^-1),
    the off-diagonal support of C is a disjoint matching.    (4)

To prove this, put h=log C, k=log E. The symmetric block of their additive
field is the collective field of (h+k)/2. Its commutation with V_s forces
h+k scalar. Thus the nontrivial part is D_h=h tensor I-I tensor h,
which exchanges the symmetric and antisymmetric sectors. Set
L=P_s D_h P_a. Block commutation gives (V_s+cI)L=L W_a. Hermiticity
then implies [V_s,LL*]=0. The computational diagonal projector D is a
spectral sum of V_s for n>=3, so [D,LL*]=0 as well.

On a diagonal matrix diag(d), LL* acts as D_h² on the symmetric sector.
Its off-diagonal (i,j) entry has coefficient -2 h_ik h_jk multiplying
d_k whenever i,j,k are distinct. Every such coefficient must vanish.
Thus each column of the Hermitian h has at most one nonzero off-diagonal
entry. The matrix is block diagonal in singleton and two-label blocks.
Exponentiation preserves that partition, proving (4), even for complex C.

Consequently any full-rank density with T(C)=gamma W_strict, gamma!=0,
is excluded: the strict source has at least two nonzero edges incident
at each vertex. The same obstruction applies separately to the dense
canonical legacy cycle; no sign, role or arithmetic transfer is used.

This matching obstruction is genuinely weaker than uniformity. The exact
three-label witness exported by research.py has

    C=[[2,1,0],[1,2,0],[0,0,1]]/5,
    E=[[2,-1,0],[-1,2,0],[0,0,3]]/7,

and a rational reciprocal pure-equivalent Hamiltonian for which C tensor E
is stationary. Its symmetric block is exactly V_s. The nonzero real source
is one isolated edge, not a dense FIN kernel. The partner's edge has the
opposite sign. Rotating this witness while retaining the same fixed T is
not an admissible basis shortcut: normal ordering is basis dependent.

## 5. What remains invisible at finite N

For (1), the uniform all-to-all sum of swaps commutes with every permutation-
invariant density, so it changes no trajectory in that sector. For (2),
the antisymmetric perturbation annihilates the completely symmetric N-body
sector. Thus all lifts (2) give the same exact dynamics on that sector,
up to global phase, not merely the same leading pure Hartree equation.
An independent mixed product C^tensorN is permutation invariant but in
general is NOT supported in the completely symmetric sector. This
distinction is why mixed-state data can reveal the additional freedom.

None of these theorems fixes a microscopic source, tensor composition,
finite-speed memory controller, observable clock, or physical particle
statistics. The rank-free theorem concerns identical products; the
unequal-partner matching theorem explicitly requires full rank.
The prior canonical theorem with arbitrary local fields is not silently
promoted to the entire pure-equivalent class.

## 6. Rank-free strict-source obstruction and a uniform correlation floor

There is a stronger, directly FIN-specific certificate requiring neither
logarithms of the marginal states nor full rank. Let C and E be density
matrices individually stationary and self-consistent with the SAME strict W:

    T(C)=gamma_C W, [W,C]=0,
    T(E)=gamma_E W, [W,E]=0, gamma_C,gamma_E>0.

Connectivity forces Re C=I/n+gamma_C W, and similarly for E. The real
uniform Perron vector u and the real alternating mode v are simple strict
eigenvectors. Every imaginary antisymmetric part commuting with W
annihilates both vectors. Hence, for j=0,6,

    C u=c0 u, C v=c6 v, c_j=1/n+gamma_C lambda_j(W),
    E u=e0 u, E v=e6 v, e_j=1/n+gamma_E lambda_j(W).

In particular Delta=c0 e0-c6 e6>0 by positivity and lambda_0>lambda_6.
Put alpha=u tensor u, beta=v tensor v, kappa=(n-2)/(2n), n=12.
For EVERY interaction Wmic in (2),

    <beta|[Wmic,C tensor E]|alpha> = kappa Delta.             (5)

Indeed V alpha=alpha/2+kappa Omega and the same identity holds for beta;
the arbitrary antisymmetric block annihilates both, while the scalar
symmetric shift cancels. Thus NO product of two stationary strict source
states can be a microscopic equilibrium in this whole class, at any rank.
The opposite-chirality minimum-rank constructions are not exceptions.

This noncancellation also holds on any finite graph containing a nonzero
canonical-source edge, with an independent arbitrary pure-equivalent
antisymmetric block on each edge and arbitrary one-site Hamiltonians.
Take the matrix element from u^tensorN to the vector with v on the two
endpoints of a chosen edge and u elsewhere. Only that edge can connect
these two vectors. Its nonzero value is (5) times its coupling and the
positive factors c0 from the other sites. Local fields cannot change both
endpoints in one matrix element. This assertion concerns products of
stationary strict source marginals, not arbitrary unequal states.

Moreover, let R be ANY stationary two-body density for Wmic, with those
same marginals, and chi=R-C tensor E. Its one-body partial traces vanish.
Define d=alpha-beta and

    B0=kappa/2 (|d><Omega|+|Omega><d|),
    L=(|u><u|-|v><v|) tensor I + I tensor (|u><u|-|v><v|),
    B=B0-kappa L/2.

Stationarity implies Tr(chi B)=-kappa Delta. The subtraction is licensed
by the zero marginals of chi. The exact operator norm is
||B||=kappa sqrt((n-2)/2): on span{d,Omega-alpha-beta} it has these two
opposite eigenvalues, and on the one-u/one-v spectator sectors it has
eigenvalues +/-kappa/2. All remaining eigenvalues vanish. Therefore

    D_tr(R,C tensor E) >= Delta/sqrt(2(n-2)).                (6)

This is uniform over the ENTIRE arbitrary antisymmetric microscopic block,
over the scalar symmetric shift and over a nonzero overall time scale.
It is not merely failure of exact equality in one numerical example.

At gamma_C=gamma_E=1/20, outward strict Fourier bounds give
Delta>63/2500, while sqrt(20)<9/2. Hence

    D_tr(R,C tensor E)>7/1250 = 0.0056.

Quantum Pinsker further gives mutual information
I(A:B)>2*(7/1250)^2=49/781250 nats. This is a necessary dimensionless
correlation resource, not a measurement, an entanglement requirement,
or a proof that an admissible stationary R exists for every interaction.
The exact numerical bound is for the declared strict loading, not silently
for legacy, other loadings, altered state laws, additional local fields
outside (2), or externally added pair terms. The preceding graph no-product
statement allows local fields; this quantitative two-body floor has the
narrower Hamiltonian scope just stated.
