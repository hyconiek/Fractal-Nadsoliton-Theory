# MH-P: a translation-uniform weighted ball discrepancy

Date: 2026-10-20. An elementary specialization of known Poisson-smoothing
methods. The general method and its exponents are not FIN discoveries.

## 1. Statement with the hypotheses actually needed

Let Lambda be a fixed full-rank lattice in R^d, d>=2. Let a>=0 range over
a family with uniformly bounded Schwartz seminorms of every fixed order.
Let radii r stay in a fixed compact subset of (0,infinity). Uniformly in
ALL lattice translations b and those amplitudes,

    |covol(Lambda) h^d sum_(y in h(Lambda+b),||y||<=r) a(y)
                       -integral_(||y||<=r) a(y)dy|
         <=C h^(2d/(d+1)), 0<h<1.                    (1)

The closed-ball convention is included. No probabilistic CLT is used.

## 2. Fourier decay for the ball, with uniform constants

For a sphere of radius in the stated compact range, a smooth surface
amplitude has Fourier transform O((1+||xi||)^(-(d-1)/2)). Here is a direct
bound sufficient for this proof. Rotate xi to the vertical axis and split
the sphere into two pole charts and a region with nonzero tangential phase
gradient. Away from the poles, repeated tangential integration by parts
gives any inverse power. At either pole the graph phase has gradient
comparable to its chart radius and an invertible definite Hessian, uniformly
in r and rotation. The cap of radius t^-1/2 has volume O(t^(-(d-1)/2)).
On smooth dyadic annuli of radius rho>=t^-1/2, M integrations by parts
with gradient/its squared norm cost at most C t^-M rho^-2M. Their volume
is O(rho^(d-1)); choose 2M>d-1 and sum. The same bound follows on both
caps. All constants require only finitely many amplitude derivatives.

Applying the divergence theorem in direction xi/||xi|| to a ball integral
gains one inverse ||xi|| and a surface term. Iterating on the remaining
interior derivative a finite number of times gives

    |Fourier(1_Br)(xi)|<=C(1+||xi||)^(-(d+1)/2).      (2)

The same argument would allow smooth weights; below we use (2) and the
Schwartz transform of a. Uniformity under rotations is paid by the uniform
seminorm bounds, not a special lattice orientation.

## 3. Positive upper/lower envelopes and Poisson summation

Take a nonnegative C-infinity mollifier rho_delta of integral one, supported
in the radius-delta ball. Put

    q_- = 1_(B_(r-delta))*rho_delta,
    q_+ = 1_(B_(r+delta))*rho_delta.

For delta small relative to the fixed minimum radius,
q_-<=1_(closed B_r)<=q_+ pointwise. Their integral error against a is O(delta)
because the difference is confined to a bounded-width annulus. Nonnegativity
of a makes a q_-,a q_+ genuine pointwise envelopes, not heuristic smoothing.

Let p=(d+1)/2. From (2), rapid decay of the mollifier transform and convolution
with the uniformly Schwartz Fourier transform of a, for arbitrary fixed M,

    |Fourier(a q_+/-)(xi)|
      <=C_M[(1+||xi||)^(-p)(1+delta||xi||)^(-M)
                                  +(1+||xi||)^(-M)]. (3)

To check the convolution bound, split its integration variable zeta at
||zeta||=||xi||/2. The near region retains both xi factors in (3); in the
far region the arbitrarily fast tail of Fourier(a) gives the second term,
using the uniform L1 norm of q. Constants are independent of delta.

Periodizing a Schwartz function on Lambda and expanding its Fourier series
proves Poisson summation. At nonzero dual lattice points xi=k/h the translation
contributes only a modulus-one phase. Absolute summation of (3), for delta<h,
is bounded by

    C[h^d delta^(-(d-1)/2)+h^M].                      (4)

Indeed the dual-point sum with p<d up to ||k||~h/delta is O((h/delta)^(d-p));
its tail has the same bound for M sufficiently large. The fixed dual lattice
has O(R^d) points in radius R, proved by disjoint fundamental-cell volume
comparison. There is no assumption about a generic or averaged translation.

Balance the smoothing loss with (4) by delta=h^(2d/(d+1)). Bracketing the
sharp sum by its two smooth envelopes proves (1). All constants depend on
the fixed lattice, radii and amplitude seminorm bounds only.

## 4. The two rates needed in FIN

For d=10 the normalized ball error is O(h^(20/11))=O(N^(-10/11)), which
is o(h). It is uniform even when two radii differ by h/8. For d=2 the
normalized disk error is O(h^(4/3))=O(N^(-2/3)), also o(h). The latter
must be applied after the exact phase/kernel residue decomposition,
not to a fictitious fully curved ten-dimensional cylinder.
