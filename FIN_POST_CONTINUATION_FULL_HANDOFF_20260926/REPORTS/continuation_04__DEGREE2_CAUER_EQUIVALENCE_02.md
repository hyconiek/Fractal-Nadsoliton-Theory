# TREE-DEGREE2-EQUIVALENCE-02 — positive-storage Cauer uniqueness

Status: **CONSTRUCTIVE CONDITIONAL THEOREM FOR A MAXIMAL DEGREE-2 PATH**.

Consider a two-terminal path with positive edge conductances g_0,...,g_m and
internal shunt storages c_1,...,c_m >=0.  Ground the right boundary and let
Y(z)=Lambda_11(z) be the left input admittance.

For positive storages the ladder recursion is

Y_0 = series(g_0, z c_1 + Y_1),
series(g,Y)=gY/(g+Y),

ending with z c_m + g_m.  Therefore

g_0 = lim_(z->infty) Y_0(z).

Define the inverse Möbius step

T_0 = g_0 Y_0/(g_0-Y_0)=z c_1+Y_1.

Then

c_1 = lim_(z->infty) T_0/z,
g_1 = lim_(z->infty) [T_0-z c_1].

Repeating the same transform reconstructs all c_i and g_i.  Thus with labelled
terminals and all g_i,c_i>0 the degree-2 path realization is unique from the
full dynamic boundary response; if terminal labels are forgotten, reversal is
the remaining path symmetry.

If c_i=0, that vertex carries no dynamic state.  Its adjacent conductances
collapse exactly to their series equivalent g_eq=g_L g_R/(g_L+g_R).  Hence the
minimal response-equivalence class suppresses precisely zero-storage degree-2
vertices.  Positive-storage degree-2 subdivisions are dynamically observable.

An exact SymPy replay reconstructs the rational fixture
(g)=(2,3,5,7), (c)=(11,13,17) coefficient by coefficient and verifies the
zero-storage series collapse exactly.

This exact uniqueness does not imply uniform numerical stability: near pole
collisions and weak residues remain ill-conditioned as quantified separately by
the Hankel/Vandermonde determinant.
