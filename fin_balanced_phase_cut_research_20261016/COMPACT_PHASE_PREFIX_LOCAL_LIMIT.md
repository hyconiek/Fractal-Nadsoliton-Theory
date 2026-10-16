# New actual cut family: compact-critical phase-prefix scores

Date: 2026-10-16. Analytical local-limit extension of the admitted MT
Stirling/lattice/tail proof; not an all-cut theorem, finite-N numeric rate
or a new process limit. Fix lambda0 and exact g_c.

## 1. A slice local limit, not just weak convergence

Write rho0(z) for the normalized exp(-z^4/12) density, and

    q_A(z)=rho0(z)[1/2+exp(-z²)/8],
    q_B(z)=rho0(z)[1/2-exp(-z²)/8].

For every fixed K, uniformly over actual parity levels z_l=N^(1/4)m_l
with |z_l|<=K and either phase P,

    pi_D(P intersect {Z_N=z_l})
       =2N^(-3/4)[q_P(z_l)+o(1)].                    (1)

Here is the additional lattice/tail payment. At a fixed even-label count,
the integer tangent lattice is the product of the two six-label zero-sum
lattices. Bases e_i-e_5 in each parity have Gram matrices I_5+J_5,
determinant 6 each. Their product has Euclidean covolume **6**. In the
orthonormal fast coordinates Y=sqrt(N)E_fast^T w the slice cell volume is
`6N^-5`. The soft mesh is `Delta z=2N^-3/4`. Their product is the actual
full MT cell volume `12N^-23/4`; no missing count multiplicity or factor
of twelve is introduced.

Uniform interior Stirling on a compact scaled (z,Y) box and the exact
entropy expansion give a weight proportional to

    N^(-11/2) exp(-z^4/12-Y^T H_fast Y/2)[1+o(1)].

For bounded z the remainder is uniform, even as the parity-level lattice
translation changes. The slice meshes have fixed shape and shrink; shifted
Riemann sums on a compact box converge uniformly in the translation and z.
The whole compact scaled box lies in the true D_N eventually by MT's
central connectivity, not a replacement cylinder.

For the expanding tails, F>=c_m m^4+c_w||w||² and the positive occupation
floor give the uniform interior Stirling upper bound with the same power
N^-11/2. The scaled slice-cell oscillation of exp(-c_w||Y||²) is uniformly
bounded, since |Y|=O(sqrt(N)), cell diameter O(N^-1/2). Thus the tail sums,
after multiplication by the slice volume, are uniformly bounded by Gaussian
tail integrals. The phase indicator has limiting boundary R=T(z), of zero
Gaussian measure uniformly on |z|<=K: the smooth positive T varies over a
compact positive interval. The same argument makes the shrinking Fourier-
shell contribution within a slice o(N^-3/4), uniformly there.

Normalizing by the full MT partition asymptotic proves (1). This derivation
uses MT's full Stirling/coercivity structure; weak convergence alone would
not prove (1), and no quantitative remainder or practical N0 is claimed.

## 2. Uniform cumulative native source

The paid native drift remainder gives stationary L1 convergence of
N^(3/4)d0 to b0(Z_N)=-Z_N^3/3. The joint phase/soft limit and actual moment
UI imply, uniformly on compact z intervals,

    N^(3/4) E_pi[1_P 1_{Z_N<=z} d0] -> J_P(z),
    J_P(z)=integral_(-infinity)^z b0(u)q_P(u)du.       (2)

For uniformity, split b0 into its positive and negative parts. Each yields
a finite positive cumulative measure converging to an atomless limit, with
its total mass paid by moment UI. A finite partition controlling the limit
mass between adjacent endpoints upgrades pointwise convergence to uniform
cumulative convergence. The L1 drift remainder is uniform in the cutoff.
No unbounded function is transferred by unspecified TV error.

The q_P are positive and even, b0 is odd and positive on the negative
half-line. Its total integral vanishes, so J_P(z)>0 for every finite z.
In particular its minimum on a fixed compact interval is positive.

## 3. The actual weighted boundary and its normalization

For S_P,l=P intersect {Z_N<=z_l}, only its top parity level contributes
ordinary within-phase upward boundary edges. Outside the Fourier shell
and departure E_N, every one of the 36 upward replacements remains in P
and in D_N. The count-degree error on the top slice is bounded by 36 times
its shell/departure mass. Cross-phase H boundary from lower slices is
incident to E_N; ellipticity bounds its normalized weighted degree by
pi(E_N)/a0. Both errors are o(N^-3/4), uniformly on fixed |z_l|<=K.
Consequently

    B(S_P,l)/W(D_N)
       =72N^(-3/4)[q_P(z_l)+o(1)].                   (3)

Since s=(N/2)W d0, (2)--(3) yield the genuinely microscopic score limit

    [s(S_P,l)]_+/(N B(S_P,l))
       -> J_P(z_l)/(144 q_P(z_l))
        =g_P(z_l)/288,                              (4)

uniformly on every fixed compact critical interval. The positive part
is eventually unnecessary there by positivity of J_P. For the whole
unpartitioned parity prefix q=rho0, J=rho0 and the limit is 1/144.
The earlier 1/9 was a valid coarse bound, not this sharp limiting value.

The inherited FP inequalities g_A<=2 and g_B<=8/3 bound the limiting
profiles, but (4) does **not** extend the finite-N compact uniformity to
all z_N growing with N. Far soft tails, departure-layer cuts and irregular
long-depth defects remain separate estimates. This theorem rules out
compact-critical phase prefixes as diverging-score countercuts; it does
not construct a diagnostic flow or refute IBA.
