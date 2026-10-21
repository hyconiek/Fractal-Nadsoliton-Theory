# A finite-rank canonical state map, without deleting the continuum

Date: 2026-10-21. Same-author proof; not formalization or independent review.
All constants below are uniform in the compact parameter interval and fixed
finite M4,M6. The model and invariant-sector scope are in [CONTRACT.md](CONTRACT.md).

## 1. The actual Hamiltonian and phase-space split

Write eta=c g+v, psi=z f+u and their canonical momenta
eta_t=r g+v_m, psi_t=p f+u_m. All four continuum coordinates remain;
v,v_m are odd and orthogonal to g, u,u_m even and orthogonal to f.
These are bounded symplectically orthogonal projections of the energy
space. In particular the kinetic energy, not only field position, is split.

Put ell=3 and w=a0+2. Subtract the static kink energy from the actual
Hamiltonian, integrate its static linear term by parts, and obtain exactly

    E=E2+E3+E4+E6,
    E2=(p^2+w z^2+r^2+ell c^2)/2+E2_cont,
    E3=integral [2K eta^3+6K eta psi^2]dx,
    E4=integral [eta^4/2+3eta^2 psi^2+g4 psi^4/4]dx,
    E6=(g6/6)integral psi^6 dx.                            (1)

The continuum quadratic form contains H on v and L on u, and both
continuum momenta. On this sector H>=3, L>=5/4; the quadratic form is
equivalent to the H1 x L2 norm, since both potentials are bounded. The
nonlinear energies are smooth on X by the one-dimensional Sobolev algebra.
No spatially integrated constant background term is individually retained.

Let F=K f^2. Direct beta integrals give

    d=<F,g>=(3pi/64)sqrt(3/2), d^2=27pi^2/8192,
    integral f^4=18/35, integral g^4=9/35,
    integral f^2 g^2=6/35, 2integral K g^3=4d.

For example integral sech^3(x)tanh^4(x)dx=pi/16. Therefore the
pure-bound cubic and quartic parts are

    E3b=d(6z^2 c+4c^3),
    E4b=(9g4/70)z^4+(18/35)z^2 c^2+(9/70)c^4.             (2)

Every term of E3m=E3-E3b or E4m=E4-E4b contains a continuum coordinate.
This is a polynomial expansion, not a physical finite-mode truncation.

## 2. Explicit cubic homological equation and signs

Use {A,B}=A_z B_p-A_p B_z+A_c B_r-A_r B_c, and put

    D A={A,E2}=p A_z-w z A_p+r A_c-ell c A_r,
    Delta=ell-4w in[-4,-2].

Define

    G3=d[A z^2 r+B z p c+C p^2 r+E c^2 r+F0 r^3],
    A=6(2w-ell)/(ell Delta), B=12/Delta,
    C=12/(ell Delta), E=-4/ell, F0=-8/(3ell^2).             (3)

Then D G3=E3b. This follows by applying D to the five displayed
monomials: D(z^2r)=2zpr-ell z^2c,
D(zpc)=p^2c-wz^2c+zpr, D(p^2r)=-2wzpr-ell p^2c,
D(c^2r)=2cr^2-ell c^3, D(r^3)=-3ell cr^2.
These identities also fix the sign of the transformation.

Let Phi3 be the auxiliary time-one Hamiltonian flow of G3. OLD physical
coordinates equal Phi3(NEW coordinates). Hence E becomes E composed
with Phi3, whose cubic term is E3+{E2,G3}=E3m. The negative sign
{E2,G3}=-D G3 is essential; using the opposite flow would double E3b.

All denominators in (3) are separated from zero on the full frozen
interval. X_G3 is a quadratic finite-dimensional vector field, acts only
on z,p,c,r, and has zero continuum components. A small ball, uniformly chosen, has
the auxiliary flow for parameters -1<=tau<=1. It is invertible and
canonical there. Extend it by identity on the continuum. The result is
a finite-rank canonical map of X, not a new evolution equation.

For clarity, the last assertion follows directly: preserve the finite
canonical two-form by Hamilton's ODE and leave the symplectically
orthogonal continuum factor unchanged. Smooth dependence and the ODE
bound |X_G3|<=C|mode|^2 give Phi3=Id+X_G3+O(||mode||^3).
Smooth profiles f,g show the same map preserves each strong Sobolev
phase space used by the actual PDE. The auxiliary parameter tau is not t.

Taylor expansion along this finite-rank flow gives

    E composed Phi3 = E2+E3m+Q4b+Q4m+R5,
    Q4b=E4b+{E3b,G3}/2,
    Q4m=E4m+{E3m,G3}.                                    (4)

Every term of Q4m still contains a continuum coordinate: G3 differentiates
only bound variables. On a small energy ball, |R5|<=C||state||_X^5 and
||dR5||_(X*)<=C||state||_X^4. The signed E6 is included in this remainder,
not removed from the exact Hamiltonian. These estimates follow from
Taylor's integral remainder for smooth energy functionals on X; all
Poisson brackets here involve only finite-rank derivatives, so no
unbounded continuum commutator/domain assertion is concealed.

## 3. Constructive quartic map on the entire parameter interval

Let V4 be the19-dimensional space of real homogeneous quartic polynomials
with even total degree separately in (z,p) and (c,r). There are five
pure-z/p monomials, five pure-c/r monomials and nine with degree2 each.
Both Q4b and D preserve this space. Let Pi be independent free-angle
averaging in the two oscillators. With

    Jz=p^2+w z^2, Jc=r^2+ell c^2,

the range of Pi is span(Jz^2,Jz Jc,Jc^2). For a monomial z^(2i)p^(2j),
its average is

    (2i)!(2j)!/[4^(i+j)i!j!(i+j)!] w^(-i) Jz^(i+j),       (5)

and odd powers average to zero. Apply the same formula to c,r with ell.
This defines the matrix Pi exactly, without quadrature or orbit fitting.

The real matrix D+Pi is invertible on V4. Indeed the pure-mode blocks
have frequencies0, +/-2sqrt(w), +/-4sqrt(w) and their ell analogues.
The mixed block has frequencies0, +/-2sqrt(w), +/-2sqrt(ell),
+/-2(sqrt(w)+sqrt(ell)), +/-2(sqrt(w)-sqrt(ell)). Thus

    det(D+Pi)=2^20 w^3 ell^3(w-ell)^2>0.                  (6)

The zero frequencies acquire eigenvalue1 from Pi. This proves uniform
invertibility on the frozen interval, not merely at three numerical points.
It uses the precise even-pair space: no exclusion of an unlisted resonance
is being assumed for a larger polynomial class.

Define constructively

    Z4=Pi Q4b,
    G4=(D+Pi)^(-1)(Q4b-Z4).                               (7)

Since Pi D=D Pi=0 and Pi^2=Pi, (7) has Pi G4=0 and
D G4=Q4b-Z4. It is the unique zero-average solution in V4, with rational
parameter coefficients multiplying the explicit source constants g4,d^2.
No optimizer, symmetry-breaking source or physical selector is inferred
from this coordinate gauge. Its auxiliary Hamiltonian flow Phi4 is local
and canonical by the same argument as Phi3, with a cubic vector field.

Writing Z4=Az Jz^2+Bz Jz Jc+Cz Jc^2 gives the explicit coefficients

    Az=27g4/(560w^2)+9d^2(8w-3ell)/(4w^2 ell Delta),
    Bz=9/(70w ell)+72d^2/(ell^2 Delta),
    Cz=27/(560ell^2)-15d^2/ell^3.                         (8)

One can check (8) simply by differentiating (2)--(3), then using
Pi(z^4)=3Jz^2/(8w^2), Pi(z^2p^2)=Jz^2/(8w),
Pi(z^2c^2)=Jz Jc/(4w ell), and their c,r analogues.
These are source-model coherent energy coefficients, not measured
constants or an averaged damping law.

## 4. Full state map and the interaction that must remain

Set OLD=Phi3 composed Phi4(NEW), with both maps identity on the continuum.
Equation (4) and the quartic homological equation now prove

    E_new=E2+E3m+Z4+Q4m+O_X(||state||^5).                 (9)

The derivative remainder has the bound specified after (4). Because
G4 has degree4, its interaction with E3m begins at degree5. In particular
the following actual cubic term is still present in E3m:

    6 z^2 <P_c F,v>.                                    (10)

Here P_c is the odd H continuum projection; it does not refer to a new
field. The quadratic drive in (10) has carrier2sqrt(w)>2, inside the
continuous spectrum of sqrt(H). No bounded L2 inverse at that carrier is
used to eliminate it. The other continuum terms remain as well. The
submanifold with all continuum data zero is not invariant under (9).
Therefore (9) is NOT a closed four-dimensional dynamical model.

For strong actual solutions, composing their state with the inverse maps
gives the canonically conjugate Hamiltonian equation, by the chain rule
and symplectic identity. This is an identity on their existing interval
of existence. It does not extend that interval or supply a uniform slow
envelope bound. The transformed canonical momenta are not generally the
time derivatives of the transformed positions; Hamilton's equation uses
the complete transformed energy, including mixed interactions.

## 5. The preparation is transformed, not refitted

Let the actual old Cauchy state be (z,p,c,r)=(epsilon,0,0,0) and all
continuum coordinates zero. Both G3 and G4 are odd under reversal of all
momenta. For G4 this follows either from (7), since D reverses this parity,
or uniqueness on the zero-average space. Both maps preserve the fixed
zero-momentum set exactly and commute with (z,p)->(-z,-p).
Consequently the actual new initial state satisfies

    Z0=epsilon+O(epsilon^3), P0=R0=0 exactly,
    C0=-d A epsilon^2+O(epsilon^4),                       (11)

with zero continuum data. The absence of an epsilon^3 term in C0 follows
from the psi-sign symmetry. This is a coordinate slip, not added radiation.
Setting C0 to zero by hand usually changes the physical preparation.

For a check using the actual retarded response, its leading shape is

    c_old(t)=-3d epsilon^2[(1-cos(sqrt(ell)t))/ell
           +(cos(2sqrt(w)t)-cos(sqrt(ell)t))/Delta].        (12)

The coefficient of the free shape cosine is -d A epsilon^2. Reconstruct
old c from the new free oscillators using

    c_old=C+d[A Z^2+C P^2]+higher terms,
    A+C w=-6/ell, A-C w=-6/Delta.

This recovers (12) with the correct initial zero field. Here C in the
square brackets is the coefficient C from (3), not the new position C.

An important exact exception occurs at w=ell/2=3/2 (a0=-1/2): A=0.
X_G3 vanishes at every raw prepared point (epsilon,0,0,0), so Phi3 fixes
it exactly. G4 is even in the shape pair, so its flow preserves c=r=0.
Thus C0=R0=0 exactly there, not just to second order. Nevertheless

    c_old(t)=-d epsilon^2[1-cos(2sqrt(w)t)]+higher terms

is nonzero. The absence of the free second-order shape cosine is not
absence of forced shape response, an invariant zero-shape theorem or
absence of higher-order excitation. Nonzero overlap d alone cannot imply
a nonzero free shape transient for every parameter in this interval.

## 6. Paid dependency delta and precise remaining gap

CM pays the explicit local canonical state map, its quartic coherent
normal form and its preparation correction in the actual field space.
It does NOT prove that the transformed actual PDE amplitude satisfies
a uniform C2 slow-envelope estimate or a long-window dispersive bound.
Those are still needed to insert the actual dynamics into OM on times
epsilon^-2. In particular the mixed cubic terms in (9) cannot be bounded
over such an interval merely by the finite-order energy jet.

No Markov closure, independent noise, nonlinear decay coefficient,
epsilon^-2 full-PDE theorem or all-time scattering is claimed. The
existing shorter NM reduction and conditional RP flux theorem retain
their own paid windows. F01 is not closed by a canonical re-description
of an action whose choice, coefficients, space and preparation are inputs.
