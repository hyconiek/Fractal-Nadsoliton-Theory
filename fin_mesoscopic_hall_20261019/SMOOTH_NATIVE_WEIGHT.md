# MH-W: smooth exact-count interpolation with uniform native payments

Date: 2026-10-20. Exact at actual count states; not a new Gibbs or kinetic law.

## 1. An analysis interpolation, including multiplicity

On a parity slice write n_j^0=N(1+sigma_j m)/12 and n_j=n_j^0+sqrt(N)Y_j,
Y in E. Extend the ACTUAL weight by

    W_N(m,Y)=const_N exp[g_c n^T A7 n/(2N)]
                                      /product_j Gamma(n_j+1),
    a_N(m,Y)=W_N(m,Y)/W_N(m,0).                      (1)

The fixed count-independent const_N is the original normalization, e.g.
N! or N!/12^N according to the source weight convention. It cancels from
(1) and is never changed between marginals. Gamma equals the factorial
on integer counts; noninteger arguments are analysis coordinates only.

Choose R_N=N^(1/20), and a fixed radial smooth cutoff equal to one on
||Y||<=R_N and zero beyond 2R_N. Use it to extend a_N smoothly by zero.
Its entire support and all fixed-step neighbors lie inside the TRUE D_N\E_N
eventually: m=O(N^-1/4), w=Y/sqrt(N), and the native energy/root-path
bounds give F<=C(m^4+||w||^2+N^-1)<h_MT/2 along a legal balancing path.
The positive cylinder is respected. This does not replace the actual D_N.

## 2. Uniform smoothness and Gaussian domination

The needed Gamma derivative estimates can be proved directly from Euler's
product. If psi is the log-Gamma derivative, then for x>0 and r>=1,

    psi^(r)(x+1)=(-1)^(r+1) r! sum_(k>=1)(x+k)^(-r-1),
    1/(x+1)<=psi'(x+1)<=1/x,
    |psi^(r)(x+1)|<=(r-1)! x^(-r).                   (2)

Uniform convergence of the differentiated series and the decreasing-
function integral test justify these statements. In particular psi'=1/x+O(x^-2).
No unproved differentiated Stirling expansion is assumed.

Both first fast derivatives of log W and of its m derivative vanish at
Y=0: every Gamma coefficient is constant within a parity block, and
g_c A7 p_m=m sigma is orthogonal to E. On E the exact fast Hessian is

    D_Y^2 log a_N=-N diag[psi'(n_j+1)]+g_c A7.       (3)

The admitted strict transverse spectral gap and (2) show that this Hessian
is uniformly negative definite on the cutoff support for large N. One may
first fix a sufficiently small m-neighborhood and then include all |z|<=1
slices in it. Both upper and lower Hessian bounds are finite and uniform.
Thus on the support a_N<=exp(-c||Y||^2), and on any fixed ball
a_N>=exp(-C||Y||^2).

Explicitly, let gamma5=1-Lambda5/Lambda6>0 and choose |m|<=gamma5/8.
For large N, 2qR_N/sqrt(N)+q/N<=gamma5/8. The lower bound in (2) then
gives N psi'(n_j+1)>=q/(1+gamma5/4)>=q(1-gamma5/4).
Subtracting g_c A7 on E leaves precision at least 3q gamma5/4.
Also n_j>=N/(2q), so the absolute Hessian norm is at most 3q. These
explicit inequalities justify both Gaussian bounds without a fitted Hessian.

All higher mixed Y/m derivatives follow from (2): two or more fast
derivatives have bounded coefficients since their N power is N^(1-r/2)
for r fast derivatives, independent of the number of fixed m derivatives.
The zero fast gradient at the centre then bounds log a_N and its m
derivative by C||Y||^2, and their fast gradients by C||Y||. Products with
the Gaussian envelope show that a_N AND partial_m a_N have uniformly
bounded Schwartz seminorms of every fixed order. Derivatives of the cutoff
are supported in its Gaussian tail and preserve those bounds.

## 3. No actual rare states are silently dropped

The admitted MT fast coercivity, positive population floor and uniform
interior Stirling bounds give, on the actual compact-z slice,

    W_N(n)/W_N(m,0)<=C exp(-c||Y||^2).

To justify the denominator without assuming a continuous Stirling theorem,
choose the actual within-parity balanced integer type. Its scaled fast
coordinate is O(h). The exact smooth Hessian and zero fast gradient already
proved imply its weight divided by W_N(m,0) is exp(O(h^2)). The inherited
INTEGER Stirling/tail comparison can therefore use this real count state
as the reference and gives the displayed bound with the same uniform constants.

The known cell-oscillation comparison bounds Gaussian lattice tails uniformly
in the six translations. After multiplying by 6h^10, weight omitted beyond
R_N is bounded by a polynomial times exp(-cN^(1/10)), hence O(N^-K)
for EVERY fixed K. The true E_N slice has the stronger inherited exponential
barrier bound. Even a factor N from the global native source cap preserves
these negligible errors. These estimates charge the states omitted by
the smooth interpolation; they are not set to zero in the actual marginals.

## 4. Smooth phase integrals and sharp lattice phase sums

Let b_N(m)=integral_(||Y||<=L) a_N(m,Y)dY and

    A_N(m,z)=integral_(12||projection_k3 Y||^2<=T(z)) a_N(m,Y)dY.

Both have uniformly bounded m derivatives; A_N also has a uniformly bounded
z derivative since T(z) ranges smoothly in a positive compact interval.
Their denominators A_N are bounded below by a positive constant: the fixed
L ball lies strictly in phase A and its density has the lower bound above.
Their upper bounds are uniform by Gaussian domination.

For the actual normalized lattice sums, the ten-dimensional ball lemma gives

    6h^10 W(ball)/W_N(m,0)=b_N(m)+O(h^(20/11)).       (4)

For phase A, use the exact 8+2 decomposition. On each of its 54 kernel
residue shifts, eight-dimensional smooth Poisson summation has error
O(h^K)(1+||x||)^(-M) for arbitrary fixed K,M. The leading term is the
kernel integral, a uniformly Schwartz function of the two phase coordinates.
The projected disk radius is sqrt(T(z)/12). Apply the d=2 sharp-disk lemma
on each of the nine phase subgrids, retaining their 1/9 cell weights. Thus

    6h^10 W(phase A slice)/W_N(m,0)
                          =A_N(m,z)+O(h^(4/3)).     (5)

Inclusive phase ties and all lattice shifts are covered by the sharp-ball
envelopes. Actual D/E and cutoff errors are smaller than either displayed
power. This supplies a RATE, not only the older weak/slice local limit.
