# ML: mixed local/global estimates and the remaining cubic trace gate

Date: 2026-11-04. Base: c13621867efa3421cacf2749d4b7b6bd4854e044.
Conditional analytical result for the unchanged supplied RP model.
Self-audited; not formal verification or independent review.

## 1. The missing estimate in CN

CN9 gives explicit localized quadratic half-wave forcing R_2loc built
from the whole physical continuous phase h=(h_phi,h_psi). Its previously
paid bound was only O(epsilon^2) in positive weighted H1. The next target
was a stronger bound involving actual LOCAL continuum energy E(t):

    ||R_2loc(t)||_(H1_N pair)<=C_N epsilon E(t).         (ML1)

E is LF's physical continuous H1-position/L2-velocity norm with weight
<x>^-3; bound modes are excluded. SP gives the global phase size
Cepsilon. Below we prove ML1 by transferring a polynomial weight through
the explicit bilinear symbols and localized linear kernels. We do not
assume the desired integrated energy estimate.

Throughout, H1_s means ||<x>^s h||H1. Integer positive-moment norms used
in QP are equivalent to these positive-weight norms for fixed integer s.
Let W(x)=(1+x^2)^2=<x>^4. The use of exponent four is a proof convenience;
H1_-4 is controlled by H1_-3 via the bounded smooth multiplier <x>^-1.

## 2. An exact polynomial identity for the bilinear operator

For B_b(f,g) with Fourier symbol b(p,q), QP paid bounded derivatives of
b and the elementary H1 x H1 -> H1 estimate. With the convention
fhat(p)=integral exp(-ipx)f(x)dx, integration by parts in the independent
variable p gives

    B_b(x f,g)=x B_b(f,g)-i B_(partial_p b)(f,g).

Repeated application yields, for any polynomial P,

    B_b(P f,g)=sum_(j=0)^deg(P)
                  P^(j)(x)/j! (-i)^j B_(partial_p^j b)(f,g). (ML2)

The derivative is taken with q FIXED in the independent (p,q)
representation, not with output frequency fixed. The q-version holds
for multiplication on the other input. There are no boundary terms on
a Schwartz core. For P=W every W^(j)/W and its first x derivative are
bounded. Put f=W^-1 h and multiply ML2 by W^-1. The finite sum of QP's
H1 bilinear estimates proves

    ||B_b(h,g)||H1_-4<=C ||h||H1_-4 ||g||H1,            (ML3)

and the same bound with inputs interchanged. Symbol derivatives through
order four suffice for this estimate. The derivative of each bounded
coefficient W^(j)/W is included in the H1 product estimate.

For the actual global H1 inputs, approximate by Schwartz functions in
both global and weighted norms. ML2 identifies the weighted limit with
the original global-H1 bilinear output. Thus ML3 does not define a
different operator or assume positive moments of the unknown fields.

Likewise the ordinary product obeys

    ||h g||H1_-4=||(W^-1 h)g||H1
                          <=C||h||H1_-4||g||H1.       (ML4)

Multiplication by K and its bounded derivative preserves H1_-4.

## 3. Localized linear kernels accept a negative-weight input

The required linear operators are

    T_alpha=[Omega^alpha,K], alpha=+1 or -1,
    T_V=Omega^-1 V, V=6sech^2(x).

Their Fourier kernels, already displayed in QP, are

    t_alpha(xi,p)=(1/i)
       [(Omega(xi)^alpha-Omega(p)^alpha)/(xi-p)] S2(xi-p),
    t_V(xi,p)=6 Omega(xi)^-1 S2(xi-p),                  (ML5)

where S2(k)=pi k/sinh(pi k/2), continuously extended at zero. The divided
difference is an integral of the derivative of the smooth massive symbol
along the segment between p and xi. All its mixed xi,p derivatives are
uniformly bounded on the compact mass band. S2 and every derivative are
Schwartz. Consequently, for any fixed a,b,M,

    |partial_xi^a partial_p^b t(xi,p)|
                                      <=C_(a,b,M)<xi-p>^-M. (ML6)

This statement holds for both kernels; the apparent divided-difference
pole is removable. In particular the moment proof is not limited to
output derivatives only: INPUT-frequency derivatives have the same bound.

For h=W f, integrate the Fourier derivatives of fhat by parts. Since
Fourier(W f)=(1-2partial_p^2+partial_p^4)fhat, the resulting operator on
f has kernel

    t_W=t-2partial_p^2 t+partial_p^4 t.               (ML7)

Its arbitrary fixed output derivatives still satisfy ML6. Young's
inequality, with <xi>/<p><=C<xi-p>, gives the H1 bound after any required
output moment. Therefore for every fixed integer N>=0,

    ||T_alpha h||H1_N+||T_V h||H1_N
                                          <=C_N||h||H1_-4. (ML8)

A Schwartz-core argument and density establish the extension; the actual
inputs also belong to global H1, where it agrees with the previous maps.
This is the extra input-weight payment missing from CN. It uses the
explicit localized kernels, not a generic assertion that every bounded
operator preserves negative-to-positive weights.

## 4. Apply the estimates to every term of CN9

The physical-to-free-half-wave map gives

    sum_i,sigma ||h_i^sigma||H1 <=Cepsilon,
    sum_i,sigma ||h_i^sigma||H1_-4 <=C E(t).           (ML9)

For the velocity contribution use Omega^-1:L2_-3 -> H1_-3. Its proof
follows from the weighted multiplier bounds used in QP/VD: bounded
symbol derivatives give positive integer moments, duality gives negative
weighted L2, and the same argument for partial_x Omega^-1 supplies the
H1 derivative. Derivatives of the polynomial weight are bounded relative
to that weight. No physical velocity is replaced by a QP phase momentum.

Now inspect CN9 term by term:

* Q is a sum of K B_b(h_phi,h_psi). ML3 and ML9 give
  ||Q||H1_-4<=Cepsilon E. Thus ML8 applied to Omega^-1 V Q_pos gives
  C_N epsilon E in H1_N.
* [Omega^-1,K](h_phi h_psi) has the same bound by ML4 and ML8.
* [Omega,K]B_b(h_phi,h_psi) has the same bound by ML3 and ML8.
* F_V,i=sigma i Omega_i^-1 V r_i has H1_N norm C_N E by ML8.
  QP's one-localized-input bilinear estimate bounds each B_b(F_V,h)
  or B_b(h,F_V) by C_N epsilon E. Multiplication by K preserves H1_N.

This proves ML1 for every branch with the stated finite source constants.
Reality, parity and all branch sums are unchanged. CN's conjugated phase
projection, phase conversion and Darboux maps preserve these estimates.
In particular P is not commuted with a free multiplier.

## 5. Absorb both localized feedback histories at once

CN11 and LF8 together give the exact compact observation decomposition

    q_psi=p_new+q_L+q_2loc+q_ge3,
    p_new=q_D+q_Q+delta_q+q_corr+q_3loc.               (ML10)

Every quantity is the previously defined actual source response or chart
correction. Their sources depend on the whole actual fields. The already
paid estimates imply

    |p_new(t)|<=C[epsilon^2+epsilon^3 sqrt(1+t)].       (ML11)

LF4 controls q_L's forcing by Cepsilon E. ML1 controls q_2loc's forcing
by the same quantity. LE and bounded observation give for their SUM

    |q_L(t)+q_2loc(t)|
                <=Cepsilon integral_0^t(1+t-s)^-1/2 E(s)ds. (ML12)

Set alpha=mu epsilon^2 and use ||h||_alpha,T^2=integral_0^T
exp(-2alpha t)|h(t)|^2dt. LF's Young estimate and its cumulative virial
conversion, now applied to ML10, give

    ||q_L+q_2loc||_alpha,T^2 <=(C_fb/mu)||E||_alpha,T^2,
    ||p_new||_alpha,T^2<=Cepsilon^2,  mu>=1,
    ||E||_alpha,T^2<=Cepsilon^2+D||q_psi||_alpha,T^2.   (ML13)

Choose once mu>=max(1,6D C_fb). The square-of-three-summands inequality
absorbs the feedback coefficient 3D C_fb/mu<=1/2. Hence

    ||E||_alpha,T^2<=Cepsilon^2+6D||q_ge3||_alpha,T^2.  (ML14)

At T=S epsilon^-2, returning to unweighted time proves

    integral_0^T E(t)^2dt
            <=C_S[epsilon^2+integral_0^T|q_ge3(t)|^2dt], (ML15)

where C_S can grow like exp(CS). The same right-hand form bounds the
whole integral of |q_psi|^2 by ML10--ML14. No all-S uniformity or physical
damping rate is inferred from this analytical time weight.

All localized quadratic feedback is now accounted for. The remaining
trace is CN9's explicit structured cubic-and-higher response, including
its quadratic cascades and actual projected input forcing.

## 6. An actual whole-continuum partial bound

CN10 and CT's free observation estimate give, for this exact q_ge3,

    integral_0^T |q_ge3|^2dt
                  <=Cepsilon^6 T^2 sqrt(1+T)<=C_S epsilon.

Together with ML15 this yields a new WHOLE physical local-continuum bound

    integral_0^(S epsilon^-2) E(t)^2dt<=C_S epsilon,
    integral_0^(S epsilon^-2) |q_psi(t)|^2dt<=C_S epsilon. (ML16)

Unlike a controlled source component, E includes both actual continuous
fields and their physical velocities. It still excludes the bound modes.
This improves the O_S(1) integrated local bound from global energy alone.
It is weaker than the required O_S(epsilon^2), and does not close the
critical modal damping theorem. CT's sharp forcing-class countercontrol
prevents claiming the extra power from generic cubic size alone; it is
not a counterexample for this specific nonlinear prepared trajectory.

## 7. Next source-specific gate and unchanged scope

Analyze CN9's full R_ge3 for a better integrated observation estimate.
Retain direct -6v^2u-g4u^3-g6u^5, the bilinear terms with F_G, and the
terms with the projected actual F_C. Seek source-specific cancellation,
frequency decomposition or spacetime control, and bound the remainder.
Do not declare the cascades negligible, assume moments of radiation, or
replace actual continuum fields by their controlled source pieces.

The three-channel modal action error remains a separate gate. No full
damping, scattering, nonlinear instability, physical source selection,
clock/geometry, strict/legacy or empirical FIN validation is claimed.
Historical source proofs and their hashes are preserved.
