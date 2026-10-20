# MH-L: exact six cosets and the rational phase projection

Date: 2026-10-20. sigma is +1 on even labels and -1 on odd labels.
Let l be the even-label total, m=2l/N-1. Set U=sqrt(N)Y. Its entries are
n_j-l/6 on the even block and n_j-(N-l)/6 on the odd block. Both block
sums vanish. Work in their ten-dimensional Euclidean direct sum E.

## 1. Full slice lattice

Let Lambda=A5_even direct-sum A5_odd, where each A5 is the integer zero-
sum lattice in six coordinates. Let rho=(5/6,-1/6,-1/6,-1/6,-1/6,-1/6)
in the ordered block and w=(rho,-rho). Since N is divisible by 12,

    U in Lambda+r w, r=l mod6, r=0,...,5.             (1)

To prove equality, subtract r from the first even count and add r to the
first odd count after subtracting their common integer means. The residual
is integer with zero sums in both blocks; the reverse construction gives
every point of this coset before positivity/domain restriction. Also
6w belongs to Lambda, so precisely these six cosets suffice.

Bases e_i-e_5 have Gram determinant 6 per block. Thus covol(Lambda)=6,
and the scaled slice covolume is 6h^10=6N^-5, agreeing with PC and MT.
This is not an assumed lattice law or a rounded parity level.

## 2. Phase coordinates and their eight-dimensional kernel

Let c_j=cos(pi j/2), s_j=sin(pi j/2). These are integer vectors:
c is supported on even labels with alternating signs, s on odd labels.
Both have squared norm 6 and are orthogonal. Write C(U)=c.U, S(U)=s.U.
On Lambda their values are independent even integers; on (1),

    C(U)=r+2u, S(U)=-r+2v, u,v in Z.                 (2)

Choose b_e=e_0-e_2 and b_o=e_1-e_3, so C(b_e)=S(b_o)=2 and the other
pairings vanish. The kernel Lambda0=Lambda intersect ker C intersect ker S
is exactly FOUR A2 lattices: within each parity block the three c-positive
and three c-negative entries (or s-positive/s-negative) each sum to zero.
It has dimension 8 and covolume (sqrt3)^4=9.

The orthonormal phase coordinates are x=(C(Y)/sqrt6,S(Y)/sqrt6). Their
projected lattice has area (2/3)h^2. The product 9*(2/3)=6 recovers the
full slice covolume, so no multiplicity is lost by projection.

At fixed (u,v), subtract the orthogonal projection from rw+u b_e+v b_o.
The remaining kernel shift is

    eta=rw+u b_e+v b_o
          -(r+2u)c/6-(-r+2v)s/6.                    (3)

Increasing u by 3 changes eta by 3b_e-c in Lambda0; increasing v by 3
changes it by 3b_o-s in Lambda0. Thus there are at most 6*3*3=54 kernel
residue shifts. Within each pair of u,v residue classes the phase grid
has area 6h^2; its weight is 1/9 of that cell area. Summing all nine
classes gives exactly the projected area (2/3)h^2.

This finite residue reduction is what permits a two-dimensional sharp-
disk estimate after an eight-dimensional smooth quadrature. It does not
assume that the whole ten-dimensional lattice is an orthogonal product.

## 3. Actual phase and ball membership

The original integer Fourier sums satisfy C(n)=sqrt(N)C(Y), and similarly
for S, since the parity means have zero c,s projections. Hence

    R_N=2[C(Y)^2+S(Y)^2]=12||x||^2,
    phase A iff ||x||^2<=T(z)/12.                    (4)

Inclusive ties remain inclusive. With p_A(z)>=1/2 and 0<gamma3<1,
T(z)>1, whereas 12L^2=3/25<1. The frozen ball and its fixed-step
neighborhood therefore lie strictly inside phase A eventually. Their
actual MT membership and distance from E follow from the native central
energy bound, not from replacing D_N by an artificial cylinder.
