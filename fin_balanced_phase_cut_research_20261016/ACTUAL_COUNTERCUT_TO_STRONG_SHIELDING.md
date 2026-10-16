# An actual count-space counterexample to the stronger shielding lemma

Date: 2026-10-16. New analytical result under the admitted exact A7/MT
premises. **This disproves the proposed stronger proof route, not (SC),
uniform flow feasibility or universal IBA.** No new microscopic law is added.

## 1. Statement and unchanged objects

Use lambda0, exact g_c, the actual MT D_N, its native departure E_N,
weights W, native opposite-parity conductances u, phase A={R<=T(Z)} and
the permitted H / forbidden F edges. Write pi=W/W(D_N), a0>0 for the
admitted native ellipticity and B(S), Fout(S), Uin(S) as in the SC order.
For a cut S, Z_F(S) are forbidden upward departure tails in S having no
upward H exit from S. Set

    D(S)=Fout(S)-Uin(S),
    Xi(S)=[D(S)]_+/(N B(S)).

There is an actual sequence of G_+-orbit-union cuts S_N, N through multiples
of 12000, in the previously unresolved irregular class, such that

    Xi(S_N)=0,
    [N W(Z_F(S_N))-Uin(S_N)]/(N B(S_N)) -> +infinity.   (1)

Thus the sufficient shielded-**weight** inequality proposed by SC is false
even on its intended remaining cut class. The original signed residual
is benign on these cuts. It must not be replaced by the stronger numerator.

## 2. A genuine forbidden upward tail of polynomial stationary weight

Put T0=T(0)>0 and k=floor(sqrt(N T0/8)). From the balanced uniform type n*
make the two legal upward moves

    y0=n*+k(e0-e2),
    y1=y0-e3+e0,
    y2=y1-e1+e0.

Their (C,S) pairs are (2k,0), (2k+1,1), (2k+2,0). Their m values are
0, 2/N, 4/N. All counts are positive for large N. We have

    R(y0)=8k²/N<=T0,
    R(y2)=8(k+1)²/N>T0>=T(Z(y2)).

Hence at least one of the two upward edges is an A-to-B crossing. Each
state and all its one-copy neighbors have |p-u|=O(N^-1/2), F=O(1/N),
and belong to the central MT component for sufficiently large N. None is
in E_N. The first crossing is therefore forbidden. Denote its A tail y_N.
This construction does not decide a transcendental tie by floating point:
phase A includes ties, and the strict final inequality is exact.

At lambda0 the admitted coercivity gives F>=0 on D_N. The elementary type
bounds give, with q=12 and the exact all-pairs convention,

    12^N (N+1)^-12 exp(-NF(n/N)) <= W(n)
        <=12^N exp(-NF(n/N)).

There are at most (N+1)^12 types. Since NF(y_N/N)<=C, and the balanced
type has F=0,

    pi(y_N)>=exp(-C)(N+1)^-24.                        (2)

For completeness, the lower type bound follows by writing a multinomial
coefficient as successive binomial coefficients. At parameter k/M the
binomial k term is a mode (check adjacent probability ratios), hence has
probability at least 1/(M+1). Multiply these bounds, using M+1<=N+1 and
the entropy chain identity. The upper bound follows since any type's own
multinomial probability is at most one. No unbound local-limit estimate
is hidden in this polynomial lower bound.

These deliberately crude powers are sufficient. They require no local
limit rate or claim that any small N already lies in the MT domain.

## 3. A deep interior singleton orbit with exponentially small weight

Fix theta=1/1000 and N in 12000N. Define

    x_N,j=N(1+theta sigma_j)/12.

These are integer counts. G_+ fixes x_N, while an odd label rotation sends
it to the distinct negative-parity type. Its C=S=0, so it and all one-copy
neighbors are in A for large N. On the parity line,

    f_theta=F(p_theta)=I(theta)-theta²/2,
    theta^4/12 <= f_theta <= theta^4/[12(1-theta²)] < h.

Thus x_N and its neighbors are strictly inside the MT cylinder/sublevel;
full local-sublevel connectivity puts them in D_N. Their fixed energy
margin puts x_N and its neighbors outside E_N eventually. Its normalized
weight satisfies

    pi(x_N)<=(N+1)^12 exp(-N f_theta).                 (3)

No forbidden edge is incident to x_N. Set S_N=A_N\{x_N}.

## 4. Exact flow cancellation and the larger shielding surplus

Odd-rotation symmetry gives s(A)=0 exactly, so

    Uin(A)=Fout(A)+Hout(A).                           (4)

Every H edge crossing the whole phase has an endpoint in E_N. Charging
conductances to that endpoint and using the native total-rate cap gives

    Hout(A)<=N W(E_N), B(A)<=W(E_N)/a0.               (5)

Split the forbidden upward tails in A into exposed E_F(A) and shielded
Z_F(A). Every exposed tail contributes its weight to B(A), hence
W(E_F(A))<=W(E_N)/a0. Each forbidden tail has at least nine safe upward
native edges, each of rate >=a0 N. Consequently its forbidden upward
outflux is at most (1-9a0)N times its weight. Equations (4)--(5) imply

    N W(Z_F(A))-Uin(A)
       >=9a0 N W(Z_F(A))-(1+a0^-1)N W(E_N).          (6)

The central tail y_N from section 2 has no neighbors in E_N, so it is
shielded for the whole-phase cut. Therefore W(Z_F(A))>=W(y_N).

Removing x_N changes no forbidden edge or forbidden-tail shielding: x_N
and all its neighbors are well away from the Fourier interface. Exactly

    Fout(S_N)=Fout(A), Z_F(S_N)=Z_F(A),
    Uin(S_N)=Uin(A)+W(x_N) sum_(up x_N->A) q0,
    D(S_N)=-Hout(A)-W(x_N) sum_(up x_N->A)q0 <=0.     (7)

All upward neighbors of x_N are in A. Native ellipticity and detailed
balance bound the newly exposed incoming-to-x_N boundary by W(x_N)/a0.
Thus B(S_N)>0 and

    B(S_N)<=[W(E_N)+W(x_N)]/a0.                     (8)

Combining (6)--(8), the stronger numerator is at least

    N[9a0 W(y_N)-(1+a0^-1)W(E_N)-W(x_N)].            (9)

After dividing by W(D_N), (2) is polynomially small, whereas the actual
departure mass and (3) are exponentially small. Equations (8)--(9) prove
the divergence in (1). Equation (7) proves Xi=0 exactly.

## 5. This is within the intended unresolved class

S_N is a G_+-orbit union but is not odd-rotation invariant. It contains
n* and excludes the m=0 type n*+(k+1)(e0-e2), which lies in B; it therefore
is not a global parity prefix. Its forbidden-tail weight is at least
W(y_N), while its B(S_N) is exponentially smaller relative to normalization,
so w_F(S_N)>B(S_N) eventually. These are actual irregular interface-insulated
cuts, not a generic graph or a finite-N numerical anecdote.

The replacement W(Z_F) erased a fixed positive part of the outgoing/incoming
cancellation. A new proof should retain **actual signed currents** and the
exact zero source of a whole phase. The valid local nine-way lemma does
not repair the stronger inequality, and its failure is not an adverse
result for the original flow problem.
