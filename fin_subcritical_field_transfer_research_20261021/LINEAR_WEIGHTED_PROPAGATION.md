# LE: weighted local energy propagation at the actual forcing regularity

Date: 2026-10-21. Conditional analytical operator proof, same author.
No appeal to a generic dispersive theorem or numerical kernel sample.

## 1. Operator, projections and statement

Let A_m=-partial_x^2+m^2-6sech^2(x), with2<=m<=sqrt(23/4).
Let Pc remove both L2 bound profiles f,g from its continuous subspace.
They are the same normalized Darboux profiles for all m. Define

    U_m(t)=( cos(t sqrt(A_m)), sin(t sqrt(A_m))/sqrt(A_m);
             -sqrt(A_m)sin(t sqrt(A_m)), cos(t sqrt(A_m)) ) Pc,

as the2x2 phase propagator on (position,momentum). The projection is
applied to both coordinates. Every square root is on the continuous
spectrum, where A_m>=m^2; no zero-mode inverse is introduced.
For real sigma use

    X_sigma(a,b)=||<x>^sigma a||H1+||<x>^sigma b||2.

For all t>=0,

    ||U_m(t)(a,b)||X_-6<=C(1+t)^(-1/2)||(a,b)||X_6.       (1)

For ODD a,b this improves to C(1+t)^(-3/2). Shape and translation modes
have NOT been claimed to decay; they have been projected separately.

For t>=1 there is the stronger rank-two phase expansion

    U_m(t)(a,b)
      =r0(x)/sqrt(2pi m t)
           ( m cos(theta) mu_a+sin(theta)mu_b;
            -m^2 sin(theta)mu_a+m cos(theta)mu_b )
          +R_m(t)(a,b),                                 (2)
    r0=(3tanh^2-1)/2, theta=mt+pi/4,
    mu_h=integral r0 h, ||R_m(t)||_(X_6->X_-6)<=C t^-3/2.

All constants are uniform in the stated compact mass interval. The input
has exactly one spatial derivative on position and zero on momentum,
plus the weights. No N16 derivative premise is hidden.

## 2. Global energy and the explicit spectral kernel

Use the actual Darboux completeness already proved in RP. Its Jost state
is e^(ikx)(3tanh^2 x-1-k^2-3ik tanh x)/[(1-ik)(2-ik)].
The continuous measure is dk/(2pi), frequency nu_m(k)=sqrt(k^2+m^2).
The product of two Jost states is

    e^(ik(x-y)) sum_(r=0)^4 a_r(x,y) k^r/[(1+k^2)(4+k^2)],

with bounded smooth coefficients and bounded spatial derivatives of
every fixed order. At k=0 it equals r0(x)r0(y).

On Pc the A_m quadratic norm is uniformly equivalent to H1, since
A_m>=m^2 on Pc and the potential is bounded. Consequently U_m(t) is
uniformly bounded in X_0, and Pc is bounded in weighted energy spaces
by its finite-rank exponentially localized subtraction formula.
Also mu_f=mu_g=0, so inserting Pc does not alter the threshold moments.

## 3. Low frequencies: an actual weighted operator expansion

Choose a smooth even spectral cutoff chi(k), equal1 for |k|<=8 and
equal0 for |k|>=9. For each phase block and the one output derivative
needed for its H1 norm, change

    k=u sqrt(2m+u^2), nu_m=m+u^2.

The low-frequency amplitude B_(x,y)(u) is smooth compactly supported.
Its first3 u derivatives have L2 norm <=C<x-y>^3, uniformly in x,y,m.
All additional tanh coefficient derivatives are bounded. Apply the
Gaussian-regularized identity

    integral B(u)e^(it u^2)du
      =e^(i pi/4)sqrt(pi/t)/(2pi)
                  integral Bhat(xi)e^(-i xi^2/(4t))dxi.

Subtract the constant term. The bound
|exp(-i xi^2/(4t))-1|<=xi^2/(4t), followed by Cauchy-Schwarz/Plancherel,
bounds the error by C t^-3/2(||B''||2+||B'''||2).
The negative-frequency phase is treated the same way.
The values B(0) give exactly (2)'s matrix, with the correct frequency
factors in all four blocks. For the derivative of the position block,
differentiate its actual spectral kernel before this argument.

After multiplying output and input kernels by <x>^-6 and <y>^-6,
the remainder is Hilbert-Schmidt: <x-y>^3<=C<x>^3<y>^3 and the resulting
squared weight is integrable in both variables. This proves the low-
frequency remainder X_6->X_-6 bound. Only the weighted L2 parts of the
input norm were needed in this low-frequency step.

## 4. High frequencies without a hidden loss of derivatives

The remaining spectral phase is t[v k+/-nu_m(k)] with v=(x-y)/t.
In the middle cone |v|<=1/2 and |k|>=8,

    |partial_k(vk+/-nu_m)|>=1/4.

Indeed8/sqrt(64+23/4)>3/4. After one output derivative for position, the
largest high-frequency amplitude order among the four phase blocks is1
(the cos position derivative or the -sqrt(A)sin block). Coefficients
and every k derivative satisfy the corresponding classical symbol bounds,
uniformly in x,y. Repeated integration by parts with
(partial_k phase)^-1 partial_k lowers the amplitude order by1 at each
step. Eight steps give a uniformly L1 amplitude and kernel O(t^-8)
throughout that middle cone.

Here is the cutoff/domain payment: first insert a smooth compact frequency
cutoff theta(delta k). Terms differentiated on its transition region
are O(delta^(8-M-1)) for initial order M<=1, hence vanish as delta->0.
The undifferentiated-cutoff term converges by L1 domination after eight
steps. This defines the ordinary smooth kernel off the propagation front
and agrees with the spectral distribution. There is no discarded
wavefront delta or derivative in the middle cone.

To control the rest, use even spatial cutoffs supported in |x|<=t/4
and equal1 on |x|<=t/8, and the same for y. If either input or output
is outside the inner region, the weight gives O(t^-6), using the global
energy bound on the high-frequency spectral propagator. For the input
H1 cutoff its derivative costs only C/t, and the weighted H1 norm
controls it. The output cutoff derivative has the same bound.
On the remaining middle piece the just-proved O(t^-8) kernel and its
needed position derivative give the weighted Hilbert-Schmidt bound.
Thus high frequencies are O(t^-6) X_6->X_-6. No input H2 or H16 norm
was smuggled into a pointwise Fourier estimate.

Combine with section3. For0<=t<=1 use global energy and bounded weights.
This proves (1)--(2). Odd input has mu_a=mu_b=0, proving its improvement.

## 5. Scope and sharpness of the leading channel

The matrix in (2) has determinant m^2 because cos^2+sin^2=1.
For even smooth input with at least one nonzero threshold moment, its
leading vector cannot vanish simultaneously in position and momentum.
Uniformly on the mass band its phase-space norm is bounded below by a
positive constant times the moment vector norm. Therefore a generic
even continuous input really has local energy order t^-1/2, not t^-3/2.
The odd gain comes from parity, not deletion of a physical mode.

This is a weighted LOCAL energy theorem. It does not imply decay of
global radiating energy, pointwise derivative traces, flux, stochastic
noise, a physical clock or a natural source for the supplied operator.

