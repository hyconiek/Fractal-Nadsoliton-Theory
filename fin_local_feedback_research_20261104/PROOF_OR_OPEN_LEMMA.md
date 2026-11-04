# LF: absorbing localized continuum feedback on critical time

Date: 2026-11-04. Base: baa489be73538429540ac50fd236cc9a2ebc8175.
Conditional analytical reduction for the unchanged supplied RP model.
Self-audited; not formal verification or independent review.

## 1. Exact spatial split after QH

Use physical eta=cg+v, psi=zf+u, with whole actual continuum fields u,v.
QH paid every bound-only term of the even source, together with BT.
Split its remaining QH10 source as R_C=R_L+R_X, where

    R_X=Pc,e[-12Kvu-6v^2u-g4 u^3-g6 u^5],               (LF1)
    R_L=Pc,e[-12K(cg u+zf v)
       -6(c^2g^2u+2cz gf v+2cg vu+zf v^2)
       -g4(3z^2 f^2u+3zf u^2)
       -g6 sum_(j=1)^4 binom(5,j)(zf)^(5-j)u^j].       (LF2)

The subscripts denote localized and spatially extended sources. The
projection is retained on both; its finite-rank subtraction is bounded
on every weight used here. Extended does not mean a new incoming wave.
Define u_L,u_X to be the zero-data responses to these exact sources,
so u=u_D+u_Q+u_L+u_X. The response u_X depends on the whole u,v through
LF1, not merely on u_X itself. This is not an autonomous reduced PDE.

Let E(t) be the norm of the actual physical continuous phase with weight
<x>^-3: the H1 position and L2 velocity norms of both u and v, combined
in a Euclidean sum. VD's weighted equivalence gives, for every finite t,

    F(t):=integral_0^t E(s)^2 ds
       <=A epsilon^2+B epsilon^4 t+D integral_0^t |q_psi(s)|^2ds. (LF3)

A,B,D are finite source-dependent constants, independent of epsilon,t.
The physical/QP phase correction is O(epsilon^2) globally and its
integrated cost is already included in B epsilon^4 t. Bound modes are
not included in E. q_psi remains VD's actual QP compact observation.

## 2. Weighted localization paid by existing energy, not new moments

Every LF2 monomial has a fixed exponentially localized f or g factor,
at least one continuum factor, and at least one other small factor.
The actual global bound controls c,z and u,v in H1 by C epsilon.
For a fixed localized h,

    ||<x>^6 h u||_2 <= ||<x>^9 h||_infinity ||<x>^-3 u||_2.

For additional field factors use their global L-infinity bounds from
one-dimensional H1 embedding. Thus quadratic terms cost C epsilon E,
cubic terms C epsilon^2 E, and quintic terms C epsilon^4 E. For epsilon<=1,

    ||R_L(t)||_L2_6 <= C epsilon E(t).                 (LF4)

Pc,e is bounded on L2_6 and does not spoil this estimate. No positive
moment of u or v, high Sobolev norm, or small shape hierarchy is assumed.
LE, followed by the bounded compact observation after the Darboux map,
gives for q_L=q(Wcal* u_L)

    |q_L(t)|<=C epsilon integral_0^t k(t-s)E(s)ds,
    k(r)=(1+r)^-1/2.                                  (LF5)

We use only the position component for this trace; LE also controls the
physical response velocity. QP's momentum is never identified with its
position derivative.

## 3. Why the direct finite-window absorption is insufficient

Unweighted Young gives ||q_L||_L2(0,T)<=2C epsilon sqrt(T)||E||_L2(0,T).
At T=S epsilon^-2 its coefficient is 2C sqrt(S), not o(1) as epsilon
shrinks. An absorption justified only by 'epsilon is small' would be
circular. Short slow-time intervals alone would also require careful
control of accumulated history. The following global weighted argument
keeps that history exactly.

## 4. Exponential time weights close this feedback, not the extended source

For mu>=1 set alpha=mu epsilon^2 and define

    ||h||_alpha,T^2=integral_0^T exp(-2alpha t)|h(t)|^2dt.

This is an analytical norm, not physical damping, a changed time variable,
a fitted coefficient, or a modification of the evolution. Multiplying LF5
by exp(-alpha t), extending inputs by zero outside [0,T], and applying
Young yields

    ||q_L||_alpha,T
       <=C epsilon ||exp(-alpha r)k(r)||_L1(0,infinity) ||E||_alpha,T
       <=C sqrt(pi/mu)||E||_alpha,T.                  (LF6)

Indeed k(r)<=r^-1/2 for r>0 and its exponential majorant has integral
sqrt(pi/alpha). The integrable singularity is just a majorant; the
original kernel is bounded at zero. Write the squared coefficient C_L/mu.

The cumulative inequality LF3 also transfers to this weighted norm.
With beta=2alpha, integration by parts gives

    integral_0^T exp(-beta t)dF(t)
      =exp(-beta T)F(T)+beta integral_0^T exp(-beta t)F(t)dt.

Inserting LF3, the constant term contributes at most A epsilon^2,
the linear term B epsilon^4(1-exp(-beta T))/beta, and the trace term
exactly its weighted integral. Therefore

    ||E||_alpha,T^2
       <=[A+B/(2mu)]epsilon^2+D||q_psi||_alpha,T^2.    (LF7)

No pointwise differential energy inequality is assumed in this step.
The nonnegative cumulative inequality is enough, including its endpoint.

BT/QH and the inverse-QP correction give the exact observation identity

    q_psi=p+q_L+q_X,
    p=q_D+q_Q+delta_q,
    |p(t)|<=C[epsilon^2+epsilon^3 sqrt(1+t)].           (LF8)

Here q_X=q(Wcal* u_X), and the sign of delta_q is fixed by the exact
inverse chart. Its absolute bound is Cepsilon^2. The component estimates
used in LF8 hold at every finite t, even when their bounds become crude.
Integrating their square against exp(-2mu epsilon^2 t) gives, for mu>=1,
epsilon<=1 and any finite T,

    ||p||_alpha,T^2<=C epsilon^2.                     (LF9)

The three contributions scale as epsilon^4/beta, epsilon^6/beta and
epsilon^6/beta^2. The constant is independent of T and epsilon.

By the square-of-three-summands inequality, LF6--LF9 imply

    ||E||_alpha,T^2
       <=C epsilon^2+(3D C_L/mu)||E||_alpha,T^2
                                      +3D||q_X||_alpha,T^2.

Choose ONCE mu>=max(1,6D C_L). This depends only on fixed source bounds,
not on the unknown solution size or time horizon. Absorb the middle term:

    ||E||_alpha,T^2<=C epsilon^2+6D||q_X||_alpha,T^2.   (LF10)

Returning to unweighted time T=S epsilon^-2 yields

    integral_0^T E(t)^2dt
       <=C exp(2mu S)[epsilon^2+integral_0^T |q_X(t)|^2dt]. (LF11)

In particular a proved O_S(epsilon^2) bound on the extended-source trace
would close the local continuum energy budget on each fixed critical
window. LF11 itself is an actual reduction, not that missing bound.
It is not uniform as S tends to infinity. The constant is not optimized
and has not been evaluated for quantitative physical predictions.

## 5. Why the remaining weighted forcing estimate cannot be assumed

Global energy smallness alone does not put Kvu in L2_6 uniformly.
Choose a nonzero smooth h supported in (-1/4,1/4), and R>2. At a fixed
time let

    v_R=epsilon[h(x-R)-h(x+R)],
    u_R=epsilon[h(x-R)+h(x+R)].

With even h these have the required odd/even parities and global H1
size Cepsilon independent of R. Their product is epsilon^2 times the
difference of the two disjoint squared packets. Hence

    ||<x>^6 K v_R u_R||_2 >=c epsilon^2 R^6.

Subtract the projections on g and f to make both fields continuous.
The projection coefficients are exponentially small in R, since f,g
are exponentially localized. The correction to this product's weighted
norm is exponentially small times a polynomial in R, and cannot cancel
the lower bound. This is a norm-class countercontrol, NOT the actual
prepared trajectory or authorization to inject a radiation packet.
It identifies the missing premise in applying LE directly to LF1.

## 6. Next exact obligation and scope

The next target is q_X for LF1. Apply a source-specific version of the
quadratic far-field transform to the continuous fields, or another
justified spacetime argument, keeping commutators with Pc, bound-mode
forcing, noncanonical phase corrections, and cubic histories. QP's
whole-state formula cannot simply be relabelled as a closed u,v formula.
A generic unweighted cubic norm is still too weak by CT; retained cubic
structure must be used. Modal action errors in all three CD channels
remain a separate problem even after the continuum estimate is closed.

This reduces the open forcing to four extended continuum-only monomials
and removes localized feedback as a separate obstruction on each fixed
critical window. It proves no full critical damping, nonlinear instability,
all-time scattering, photon law, natural action/clock/geometry, F01,
strict/legacy bridge or selector closure.
