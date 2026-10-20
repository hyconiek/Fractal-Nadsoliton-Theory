# Causal and spectral audit on the actual proof bytes

Date: 2026-10-20. Same-author full-input audit, not an independent review.

## 1. Darboux normalization and source transform

The inherited proof is sufficient: with D_j=partial_x+j tanh x and
T=D_2^*D_1^*, T^*T=(H_f+1)(H_f+4). Thus T times the inverse positive
square root is an isometry. Solving D_1D_2 h=0 gives precisely the L2
span of sech^2 x and sech x tanh x. This pays the missing two modes and
the continuous range, rather than inferring completeness from a Jost ODE.
The denominator (1-ik)(2-ik) differs from the positive Fourier normalization
only by a unit phase, so dk/(2pi) is retained.

An additional direct source check uses y=tanh x and the inherited transform
recurrence F_6=(k^2+16)F_4/20. In the conjugated Jost numerator,

    integral tanh sech^4 exp(-ikx)=-ik F_4/4,
    integral tanh^3 sech^4 exp(-ikx)=-ik F_4/4+ik F_6/6,
    integral tanh^2 sech^4 exp(-ikx)=F_4-F_6.

Their combination is ik(k^2+4)F_4/8. Multiplication by (9/4)/Dbar
reproduces the source M with its sign and complex conjugation. Both
spectral roots give Gamma=|M(q)|^2/(2q). No factor of two is missing.

## 2. Zero-data causality and the pole signs

In each oscillator solve z''+nu^2 z=-exp(-iOmega t) by the zero-data
Duhamel integral. Its exact value is the source formula (5). The combined
numerator vanishes at nu=Omega; its first time derivative does too.
Splitting with nu^2-Omega^2-i0 is therefore legitimate and changes neither
data nor forcing. The exp(+inu t) coefficient is the SMOOTH
1/[2nu(nu+Omega)]; the exp(-inu t) coefficient keeps both poles.

The partial-fraction prescriptions are opposite at the two roots:

    1/(k^2-q^2-i0)
       =[1/(k-q-i0)-1/(k+q+i0)]/(2q).

The outgoing Green function i exp(iq|x|)/(2q) has derivative jump -1,
as required by -partial_x^2-q^2. With forcing -S the positive-side
amplitude is therefore -i M(q)/(2q). Reversing that sign would give the
wrong Green jump. The source sign agrees with this separate test.

## 3. Explicit uniform cutoff geometry

Take pole cutoff radius delta=1/16 at +/-q, q in[1,sqrt3]. On the positive
pole support k>=15/16 and

    k/sqrt(k^2+4)>=15/sqrt(1249)>5/12.

The last strict inequality follows by squaring: 32400>31225. For every
v=x/t in[1/6,1/3], the transient phase kv-nu(k) consequently has derivative
less than -1/12 there. On the negative pole support its derivative is
greater than 1/6+5/12=7/12. Hence neither transient pole contributes a
selected residue; the static exp(ikx) pole at +q does. The stationary
points have |k_s|<=1/sqrt2<3/4, leaving a gap from the pole supports.
These constants are uniform in the actual compact parameter interval.

For smooth complementary amplitudes, the one-dimensional stationary-phase
argument in the source proof is valid. Near the unique critical point its
second derivative 4/(k^2+4)^(3/2) is bounded away from zero on a fixed
compact neighborhood. At infinity v stays below one, so the first phase
derivative stays separated from zero. The source M and every needed
polynomial multiple decay rapidly; cutoff complements preserve those
uniform seminorm bounds. No unbounded-frequency or threshold endpoint
has been hidden in the error.

## 4. Differentiation and all non-radiative pieces

Partial t/x derivatives are applied to the EXACT split integrals, not to
an error stated only along a ray. Their factors Omega, nu and k preserve
the same amplitude classes and the residue sign tests. The selected residue
then differentiates to the advertised w_t,w_x. The Jost correction and
its x derivative are uniformly O(exp(-2x)); the UNSPLIT oscillator gives
their full contribution O((1+t)exp(-t/3)) on the band, also for u_t.

The DC static part is rapidly decaying on x~t; its transient has the same
smooth stationary-phase bound. Odd shape-mode value/derivatives decay as
exp(-x), with uniformly bounded coefficients since Omega^2-3>=2. The
translation projection is zero by the actual odd source, not artificial
pinning. At k=0, nu and M have smooth removable values on the full real
line. Thus the claimed O(t^-1/2) error for u,u_t,u_x is paid uniformly.

No specific contradiction or missing analytical inequality was found in
these source steps. This disposition remains same-author, conditional on
the supplied canonical model and inherited normalized spectral proof.
