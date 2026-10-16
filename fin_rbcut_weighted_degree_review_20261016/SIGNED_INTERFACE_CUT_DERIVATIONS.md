# New bounded proof support: signed incoming compensation and three paid cut classes

Date: 2026-10-16. Analytical consequences of the actual admitted model,
not a constructed flow or an all-N cut theorem. Frozen lambda0, exact g_c,
N in 12N, actual MT domain D_N, weights W and the unchanged phase predicate.

## 1. Exact residual, including zero-denominator cuts

Orient the native opposite-parity edges U upward in m, let H be the permitted
rare-departure support, and F=U\H. Let u_e be native conductances and
`s=A_U u=(N/2)W L0_D m`. For any union S of G_+ orbits, define

    B(S)=sum_(x in S) W_x # {upward H edges x->S^c},
    Hout(S)=sum_(H edges S->S^c) u_e,
    Fout(S)=sum_(F edges S->S^c) u_e,
    Uin(S)=sum_(U edges S^c->S) u_e,
    D(S)=Fout(S)-Uin(S).

Summing the exact native divergence gives **s(S)=Hout(S)+D(S)**. Put

    Xi_N(S)=[D(S)]_+/(N B(S)).

Use zero/zero=0 and positive/zero=infinity. Since Hout>=0 and each native
rate is at most p_max N,

    Xi_N(S) <= [s(S)]_+/(N B(S)) <= p_max+Xi_N(S).       (1)

For B=0, Hout=0, so the two ratios have the same zero/infinity value;
the statement includes capacity-free obstructions. Incoming **native**
crossings must remain in D(S), even if those edges are themselves forbidden.
A small unsigned Fourier-shell mass cannot replace their cutwise cancellation.

## 2. Bulk-exposed cuts

Let `w_F(S)=sum W_x` over distinct tails x in S with a forbidden upward edge
to S^c. These tails lie in the Fourier shell. The native total rate bound
gives Fout<=N w_F(S), without multiplying by the number of labels. Thus

    w_F(S)<=B(S)  =>  [s(S)]_+/(N B(S))<=1+p_max.       (2)

If both vanish, the source cannot be positive and the score is zero. This
settles this whole class, without an anti-concentration rate or a weight
normalization. A uniform inequality w_F<=C B on *all* cuts is not necessary:
closed whole-phase cuts can have B=0, w_F>0 and s=0 through incoming balance.
The independent finite fixture demonstrates this distinction, not actual
FIN feasibility. It is a reason **not** to make unsigned expansion the
new compulsory premise.

## 3. Odd-rotation-invariant cuts

An odd label rotation R leaves W,D_N and the phase predicate invariant,
while `s(Rx)=-s(x)`. Therefore every cut satisfying RS=S has s(S)=0 exactly.
This includes the entire phases A_N and B_N. Their zero scores are already
paid even where an unsigned interface quotient would be infinite.

## 4. Global parity prefixes: a quantitative native bound

For sufficiently large admissible N, each prefix S={m<=a} has score <=1/9.
Here is the additional payment beyond the old sign-only prefix statement.

Only the top occupied parity level can have upward edges leaving a prefix;
no native upward edge enters it. Outside E_N the previously proved phase
geometry supplies at least one same-phase k3 increment (Delta C,Delta S)
in {(+/-1,+/-1)} for the upward parity direction. Each such increment is
realized by **nine distinct** odd-to-even label pairs: each of the relevant
Fourier-axis values has three labels. All source labels are occupied in
the MT cylinder, and all neighbors are in D_N outside E_N. Hence there are
at least nine allowed upward exits at every such top-layer vertex.

At these vertices the total upward native rate is <=N and therefore at
most `(N/9) #H_up(x,S^c)`. At a top-layer vertex in E_N every internal native
edge is allowed; the same inequality follows edgewise from
`q0(x,y)<=p_max N<N/9`. Summing with W_x proves

    s(S)=Uout(S)<=(N/9)B(S).                          (3)

Empty outgoing sets give zero. The geometric threshold in the inherited
direction lemma is only part of the sufficiently-large-N qualification;
this is not a new practical MT N0. The proof is for full global prefixes,
not arbitrary fibre-dependent lower ideals or intersections of phases with
prefixes. Nine label pairs is not nine independent physical channels.

## 5. The genuinely unresolved class

Let R_N consist of G_+-orbit unions which are not global parity prefixes,
are not invariant under the odd rotation, and satisfy w_F(S)>B(S). Combining
(1)--(3), with the supremum of an empty class set to zero,

    sup_(S in R_N) Xi_N(S) <= Q_N
       <= max {1+p_max, p_max+sup_(S in R_N) Xi_N(S)}.  (4)

Together with native ellipticity, a uniform bound on the left residual is
equivalent to the original uniform vertex-cap existence problem. Three
classes of cuts are now paid; this is not a proof for the irregular,
interface-insulated class. Its exact source and incoming conductances cannot
be replaced by cubic drift, limiting Gaussian probabilities or an o(1)
absolute defect. A valid growing-N bad-cut family in this class would also
settle a genuine obstruction, with its stated N quantifiers.

The next proof attempt should establish **signed incoming compensation**
for this remaining class, or exhibit an actual bad-cut sequence. Another
abstract norm equivalence, an unsigned-shell inequality contradicted by
balanced cuts, or a small synthetic network would not pay this dependency.
