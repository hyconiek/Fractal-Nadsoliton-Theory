# RT: exact reference transport, equivalence and composition

Date: 2026-10-20. Self-audited finite-model mathematics. Existing positive
scaling, uniqueness and reversible compression results are reused. No
mathematical-priority, independent peer-review or physical-source claim.

## 1. Inherited source and a complete reference-equivalence test

For fixed positive pi and a finite symmetric reference K>=0 with positive
diagonal and connected support, the preceding hypothesis gives the unique
positive scaling and its reversible source law

    pi_x=z_x(Kz)_x, c^K_xy=K_xy z_x z_y,
    P^K_xy=c^K_xy/pi_x.                                 (1)

The variational objective is stationary history entropy rate plus
sum c_xy log K_xy. For the original 0/1 reference the latter is zero on
legal paths. K is a reference weight, **not** a jump covariance tensor.

Let K and R have the same support and positive diagonals. At any fixed
positive pi the following are equivalent:

1. Their selected step kernels coincide.
2. R_xy=b_x K_xy b_y for a positive vertex vector b.
3. On every allowed off-diagonal edge,

       R_xy^2/(R_xx R_yy)=K_xy^2/(K_xx K_yy).            (2)

Proof: equality of kernels with the same pi gives equality of their flows.
If z and u are their positive scalings, K_xy z_x z_y=R_xy u_x u_y,
so b_x=z_x/u_x proves 2. Conversely the reference objective changes by
the fixed constant 2 sum pi_x log b_x under 2, and its unique maximizing
law is unchanged. Diagonals determine b_x=sqrt(R_xx/K_xx); positivity
makes (2) precisely the remaining off-diagonal equations. Different
supports cannot select the same kernel, since all legal entries of the
maximizing flow are positive. Thus a support check and finitely many (2)
are a complete equivalence certificate, not just a sufficient test.

The invariants can also be read from the selected flow:

    c_xy^2/(c_xx c_yy)=K_xy^2/(K_xx K_yy).              (3)

This identifies reference classes conditional on the supplied model. It
does not source K from nature. Equality here concerns an entire discrete-
step law, not equivalence after arbitrary waiting-time changes. In
particular different holding probabilities can share a jump-order law.

## 2. A derived coarse reference, with actual hypotheses

Partition the states into nonempty fibres F_i and fix w_x>0 before fitting
any dynamics. Put W_i=sum_(x in F_i)w_x and

    V_xi=1_(x in F_i) sqrt(w_x/W_i), V^T V=I.

Require the **weighted block-invariance condition**

    K V=V Kdown, Kdown=V^T K V.                         (4)

Equivalently, for each i,j the weighted block sum

    B_ij(x)=sum_(y in F_j)K_xy sqrt(w_y/w_x)

is independent of x in F_i, and

    Kdown_ij=sqrt(W_i/W_j) B_ij.                        (5)

For symmetric K, W_i B_ij=W_j B_ji, so Kdown is symmetric. Its diagonal
is positive, and its support is connected when the fine support is.
With w=1, (4) is the usual constant block-row-sum condition. The formula
V^T K V always exists; its **dynamical** use needs (4), not just compression.

Now supply pi_x=p_i w_x/W_i with p_i>0. Let psi solve

    p_i=psi_i(Kdown psi)_i.

Set z=V psi. By (4), z_x(Kz)_x=(w_x/W_i)p_i=pi_x.
The uniqueness already proved makes this the actual fine scaling. Summing
its source rows over a target fibre gives

    sum_(y in F_j)P^K_xy=Kdown_ij psi_i psi_j/p_i,
                          x in F_i.                    (6)

Thus the **full projected history is Markov**, for every initial source
law, with exactly the kernel selected from the derived reference Kdown
and coarse stationary p. No stationary-rate replacement is used.

This transports the existing hypothesis; it does not impose another bare
maximum-entropy rule on the coarse carrier. The conditional-pi form and
reference masses are essential, not tunable corrections after failure.

## 3. A necessary criterion when fibres have only self transitions

The weighted block condition is sufficient in general, not claimed
necessary for arbitrary fibres with internal transitions. There is an
important complete subclass. Suppose K_xy=0 for distinct x,y in the same
fibre. Canonically normalize the diagonal by

    H_xy=K_xy/sqrt(K_xx K_yy), H_xx=1.                 (7)

This is a vertex gauge transformation, so it leaves the source law at pi
unchanged. For fixed w and pi_x=p_i w_x/W_i, the following are equivalent:

- H V=V(V^T H V).
- The selected source is strongly lumpable for one positive p.
- It is strongly lumpable for every positive p of that conditional form.
- Its stationary projected history is Markov for one such positive p.

To prove necessity, work with H. Since a fibre has no internal off-diagonal
edges, its source probability of staying in that fibre is P_xx=z_x^2/pi_x.
Strong lumpability makes this a positive constant r_i on F_i. Therefore
z_x=sqrt(w_x/W_i) psi_i, psi_i=sqrt(p_i r_i). In every target fibre j,
the lumped row formula is a positive constant times B^H_ij(x). Its
constancy forces every B^H_ij(x) to be independent of x, which is (4).
Sufficiency for all p is section 2.

The final equivalence uses the already paid reversible compression fact,
not a new Markov axiom. Let E be conditional expectation onto fibre
functions in L2(pi), and T the reversible fine step operator. Stationary
Markov output implies C_2=C_1^2, whereas

    E T^2 E-(E T E)^2=E T(I-E)T E=A^*A.

Hence A=0. The fibre subspace is invariant under T, giving strong
lumpability because pi is positive on every fine state. The converse
follows directly from the constant fibre transition rows.

The diagonal normalization in (7) is needed: multiplying an invariant
reference by a vertex gauge varying inside a fibre may destroy its raw
block sums while preserving its source law. Canonical H removes that
irrelevant freedom. Internal edges invalidate the necessity proof; an
explicit guard is given below. That diagonal-normalization test alone is
not extended to arbitrary internal edges. The next result uses a different
canonical gauge.

### A complete criterion with internal edges: balance before transport

Put W=sum_x w_x, pi0_x=w_x/W, and let z0 be the unique positive scaling
of K at pi0. Define the **balanced reference**

    R0_xy=K_xy z0_x z0_y/sqrt(w_x w_y).                (7b)

It obeys R0 sqrt(w)=sqrt(w)/W and is unchanged by any vertex gauge of K,
since z0 transforms by the reciprocal gauge. The base source is
P0_xy=W R0_xy sqrt(w_y/w_x).
The base law pi0 is an auxiliary normalization, not a physical preparation
or controller. Rescaling all w by one constant rescales R0 reciprocally
and leaves its reference class and selected laws unchanged; no clock or
absolute physical scale is fixed.

For an arbitrary partition, including internal edges, these are equivalent:

1. R0 satisfies weighted block invariance (4).
2. The source at pi0 is strongly lumpable.
3. The source at pi_x(p)=p_i w_x/W_i is strongly lumpable for at least
   one positive p.
4. It is strongly lumpable for every positive p of that form.
5. The stationary observed history is Markov for at least one such p.

The equivalence of 1 and 2 follows from the formula for P0. Condition 1
implies 4 by section 2 and gauge invariance, while 4 implies 3. For the
nontrivial converse, suppose 3 and let z be that source's scaling. Set
R_xy=K_xy z_x z_y/sqrt(w_x w_y), a fixed vertex gauge of K. Then

    sum_(y in F_j)R_xy sqrt(w_y/w_x)
       =(p_i/W_i) sum_(y in F_j)P^K_xy, x in F_i.

Strong lumpability makes every block sum independent of x, so R satisfies
(4). Section 2 gives lumpability for every p, and gauge invariance transfers
this back to K. Thus 3 implies 4 and in particular 2. The equivalence with
5 is the reversible compression argument already paid above. No new
reversibility axiom is being added: the selected source is reversible.

This is a complete conditional test, not a selected-p fit. It requires
only the base scaling at pi0 fixed by w, before any target coarse p is
chosen. When it passes, the gauge-invariant coarse reference
R0down=V^T R0 V is derived, and it is already balanced for coarse masses
W_i: R0down sqrt(W_i)=sqrt(W_i)/W.

When raw K itself satisfies (4), balancing and transport commute:

    B_(W_i)(V^T K V)=V^T B_w(K)V.                      (7c)

Indeed both references select the same coarse family by section 2 and
therefore are gauge equivalent by section 1. The right side is already
balanced, which fixes that gauge representative uniquely. Here B_w denotes
(7b); it is not the block-count matrix or the native HB generator.

The complete criterion does **not** assert that raw K, or diagonal-normalized
H, must be block-invariant. The internal-edge counterexample below only
rejects that incorrect shortcut. If the balanced test fails at pi0, no
choice of positive coarse p with that same conditional fine measure repairs
Markov closure. This is a structural result at a fixed reference, not a
universal no-go for every microscopic dynamics.

## 4. Nested reductions compose only with transported masses

For a second partition of the indices i into G_a, retain the reference
masses W_i. Put M_a=sum_(i in G_a)W_i and

    U_ia=1_(i in G_a)sqrt(W_i/M_a).

Then the product VU is exactly the direct embedding of the union fibres:

    (VU)_xa=1_(x in union_(i in G_a)F_i)sqrt(w_x/M_a).

If the first and second weighted block-invariance conditions hold,

    K V=V K1, K1 U=U K2,
    K2=U^T K1 U=(VU)^T K(VU).                          (8)

Thus direct and staged transported references coincide. If the supplied
fine law is pi_x=P_a w_x/M_a, then the intermediate law is
p_i=P_a W_i/M_a; using these **same** laws, direct and staged selected
source kernels and all projected path laws agree by section 2.

The second embedding generally is NOT uniform over the intermediate
indices. Resetting all W_i to one changes the conditional measure and
usually the reference class. This is a source of actual law changes,
not merely an error in matrix display. Formula (8) alone, without the
invariance and conditional-pi hypotheses, is just static algebra.

## 5. Actual N2/q12: labels -> difference -> inversion orbit

At the labelled level take the unchanged one-slot-change reference, one
unmarked self continuation and w=1. The q difference fibres all have
mass q. Direct neighbour counting gives B_dd=1 and B_dv=2 for v!=d,
so K1=2 ones-I. This derives the preceding source's multiplicity term;
it is not a new kernel choice.
The fine reference has constant degree 23 at q=12. Its balanced version is
therefore K/3312; both transported references inherit the same factor.
No new coarse normalization is selected to fit HB.

For q=12 the inherited symmetric circulant A7 gives p_d=p_(-d). Reduce
differences to inversion orbits r=min(d,12-d), r=0,...,6. Their sizes are

    ell=(1,2,2,2,2,2,1).

The conditional-pi condition holds. The intermediate masses are 12 per
difference and become 12 ell_r on the radial fibres. The second block
counts and derived symmetric reference are

    B2_rs=2 ell_s-1_(r=s),
    K2_rr=2 ell_r-1,
    K2_rs=2 sqrt(ell_r ell_s), r!=s.                    (9)

This reduction is strongly lumpable and composes with the labelled one.
The diagonal 3 in a two-element orbit includes real microscopic changes
d<->-d that are invisible in the reduced record. They are not simply
null attempts. Retaining them is part of the proof.

If the scaling psi_d is constant on an orbit (as follows from symmetry),
write its common value psi_r and the per-difference probability p_r.
The reduced probabilities and actual-change law are

    P2_rr=(2 ell_r-1)psi_r^2/p_r,
    P2_rs=2 ell_s psi_r psi_s/p_r, r!=s,
    J2_rs=ell_s psi_s / sum_(t!=r)ell_t psi_t.           (10)

Native HB instead has J2_HB,rs=ell_s p_s/sum_(t!=r)ell_t p_t. For two
different p_s,p_t the inherited strict odds inequality for psi survives
the multiplicity factors. Since there are seven orbits, a distinct current
orbit can be chosen. Thus the nonuniform hypothesis/HB branching
disagreement survives forgetting the orientation of the difference.
These are internal label orbits, not emergent physical positions or a
selector/orientation source.

For the exact non-A7 fixture a=(1,2,3,4,5,6,7,6,5,4,3,2), sum a=48,
D=4378 and p_d=a_d(96-a_d)/D. At current radial class r=2,

    J2_(2,0)=1/42,
    J2_HB_(2,0)=19/764,
    difference=-17/16044.                              (11)

The exact number is a fixture only. The actual A7 statement is qualitative
and conditional on its inherited static p, with no numerical A7 fit.

### The count-state reference for any finite number of labelled slots

There is also an exact all-finite-N consequence of the same hypothesis,
not a new microscopic law or a large-N sweep. For S={0,...,q-1}^N with
single-slot moves and one unmarked self step, partition by counts n.
The orbit size is m_n=N!/prod_a n_a!, and the supplied exchangeable Gibbs
law is uniform within each orbit. From any state with counts n there are
exactly n_a continuations to n'=n-e_a+e_b, a!=b. There are no internal
off-diagonal transitions in a count fibre. Thus

    Bcount_(n,n')=n_a,
    m_n n_a=m_n'(n_b+1),
    Kcount_(n,n')=sqrt(n_a(n_b+1)), Kcount_(n,n)=1.      (11b)

This is the transported symmetric reference, not a native HB rate or a
jump-covariance Gamma. With actual multinomial Gibbs count weights pi_n,
its selected count source is

    Pcount_(n,n')=Kcount_(n,n')psi_n psi_n'/pi_n,
    pi_n=psi_n(Kcount psi)_n.

Here the inherited all-pairs law is
pi_n proportional to [N!/prod n_a!] exp[(g/(2N))n^T A7 n]. The constant
self-pair convention remains included. This equilibrium input is reused,
not independently sourced by the transport theorem.

Full permutation invariance supplies the needed conditional-uniform law;
block counts above pay closure without enumerating large count lattices.
The square-root factors are finite orbit normalization, not quantum
creation operators, physical particles or a physical Hamiltonian.
At the auxiliary uniform-labelled base law, degree d=N(q-1)+1 gives
Kcount sqrt(m)=d sqrt(m), so the balanced count reference is
Kcount/[q^N d]. This normalization does not select nature's initial law.

For actual q12,N=2, the 144 labelled states can be reduced either through
12 differences or through 78 unordered count types, and then to the same
seven inversion orbits. The count-type masses are 1 or 2, not all equal.
Keeping these masses and the inherited p_d=p_(-d) makes both routes satisfy
section 4 and yield exactly (9)--(10). Thus a genuine two-route state-map
consistency square is paid. Directly maximizing bare entropy on count
states would instead drop (11b) and is not this derived law.

No conditional local-domain reflection, native mixing/trapping theorem,
finite-N time calibration or infinite-population limit is transferred to
the new hypothesis by (11b).
The finite-N formula is a counting corollary, not an A7/N/g search for a
larger system in which the rejected HB-source claim might become true.

## 6. Exact guards and counterexamples

### Nested unequal masses cannot simply be flattened

Take 50 fine states in blocks of sizes 9,16,25. Fine references are constant
on each pair of blocks, with block coefficients

    [[25/9,1,1],[1,2,1],[1,1,3]].

The first transported reference, second embedding and correct final
reference (merging the first two blocks) are

    K1=[[25,12,15],[12,32,20],[15,20,75]],
    U=[[3/5,0],[4/5,0],[0,1]],
    K2=[[41,25],[25,75]].

All invariance equations hold exactly. For final law (33/83,50/83), the
fine maximizing scaling is constant 1/sqrt(4150), making a completely
rational stationary source. Direct and staged source rows agree.

Flattening the intermediate masses would instead use equal 1/sqrt(2)
entries in the first column. Its compressed reference has diagonal 81/2,75
and squared off-diagonal 1225/2. The edge invariant is 49/243, not the
correct 25/123. By section 1 this is not a harmless vertex gauge and
selects a different law at the same final positive p. The attempted
intermediate uniform-conditional assumption also fails (actual masses are
9:16). No false dynamical reduction is licensed by a static compression.
This certificate distinguishes discrete-step laws; the final binary
first-change destination alone is trivial and is not a clock-free test.

### A compressed matrix can exist while the output has memory

Take the unweighted three-state path with self loops, uniform pi=1/3,
and fibres {0,1},{2}. Put t=(sqrt(5)-1)/2. The selected source is

    T=[[t,1-t,0],[1-t,2t-1,1-t],[0,1-t,t]].

Using t^2=1-t verifies stochasticity, pi and the scaling certificate.
The first fibre's rows to the second are 0 and 1-t, so strong lumpability
fails. Its stationary compressed kernels have

    tr(C2-C1^2)=3(2-3t)/4>0.                            (12)

The positivity follows from t<2/3. Hence the output is not Markov. This
is a finite countercontrol to automatic projection, not a physical memory
device or a new native q12 diagnostic law.
It also fails the balanced test at pi0=uniform, so the complete criterion
shows that changing only the positive coarse mass, keeping the within-fibre
uniform profile fixed, cannot make this source projection Markov.

### Internal edges prevent an unrestricted necessity claim

Let K=T be the symmetric stochastic matrix

    [[3/5,1/10,1/5,1/10],
     [1/10,1/2,3/10,1/10],
     [1/5,3/10,2/5,1/10],
     [1/10,1/10,1/10,7/10]].

At uniform pi=1/4 the selected scaling is z=1/2, so the selected source
is T itself. Fibres {0,1,2},{3} are strongly lumpable. Yet canonical H
has squared entries to state 3 equal to 1/42,1/35,1/28, so H is not
block-invariant with w=1. The first fibre has internal transitions, exactly
the excluded case of section 3. This weighted-reference guard is not an
unweighted A7 counterexample or a claim that block invariance is necessary
for every possible observation.

## 7. Result and remaining boundary

The new object is the conditional reference transport (K,w,partition)->Kdown,
with its mass law, composition rule, complete step-law gauge test and the
balanced all-positive-marginals closure criterion. It
pays a structural dependency of the authorized hypothesis without changing
that hypothesis to recover HB. It does not choose the initial reference,
carrier, partition, stationary weights or physical time from old FIN.
Strict F01/F03/F06 and all physical/ontology/selector/legacy-role gates
remain open. No natural geometry, conserved physical energy, calibration,
QM/GR or ToE closure is exported.
