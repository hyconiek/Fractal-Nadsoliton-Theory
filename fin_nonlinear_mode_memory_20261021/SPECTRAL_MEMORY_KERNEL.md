# The explicit memory: bound shape, continuum tail and passive response

Date: 2026-10-21. Reuse the accepted normalized H Darboux resolution;
the new object is its ACTUAL nonlinear feedback profile F=K f^2.

## 1. Profile, moments and the nondecaying piece

Here F=(3/4)tanh x sech^4 x is odd. It is exactly orthogonal to the even
translation mode. The normalized odd shape is
g=sqrt(3/2)sech x tanh x, Hg=3g. Direct beta/polynomial integrals give

    ||F||^2=2/35,
    d=<F,g>=(3pi/64)sqrt(3/2), d^2=27pi^2/8192,
    integral f^4=18/35, integral f^6=24/77.                (1)

One further exact check is H^(-1)F=(1/8)sech^2 x tanh x on the odd sector,
so the static susceptibility is <F,H^(-1)F>=1/70. This is a gapped odd
inverse, not an inverse on the translation mode.

For example integral sech^5 x tanh^2 x dx=pi/16. Translation and shape
are BOTH accounted for; removing only the zero mode would be insufficient.

The causal kernel is

    k(t)=(d^2/sqrt(3))sin(sqrt(3)t)+k_c(t).               (2)

The shape was not an independent initial wave. The actual quadratic
forcing excites it through this nonzero overlap even from eta=eta_t=0.

## 2. Continuum transform in the inherited convention

Use the source Jost functions with continuum measure dk/(2pi). For this
normalization write M_F(k)=integral conjugate(psi_k)F. The already paid
Fourier recurrence gives, with F4(k)=pi k(k^2+4)/[6sinh(pi k/2)],

    M_F(k)=3ik(k^2+4)F4(k)/[32(1+ik)(2+ik)],
    |M_F(k)|^2=pi^2 k^4(k^2+4)^3/
                    [4096(1+k^2)sinh^2(pi k/2)].         (3)

The removable value obeys M_F(k)/(ik)->1/4. The source's old second-
harmonic S=3F has M_S=3M_F; no new Gamma digits are being calculated.
At Omega^2=4+q^2 the both-root normalization gives exactly

    Gamma_old=9 Im<F,(H-Omega^2-i0)^(-1)F>.               (4)

At the frozen Omega^2 in[5,7], the shape contribution d^2/(3-Omega^2)
is real. It does not change the imaginary part used by the old P_c Gamma.

Thus the relation to the admitted old spectral density retains its factor
nine, retarded sign and both roots, rather than silently redefining Gamma.

With nu(k)=sqrt(k^2+4),

    k_c(t)=integral_R |M_F(k)|^2 sin(t nu(k))/nu(k) dk/(2pi). (5)

All integrands are smooth at k=0 and decay exponentially at infinity.
Parseval gives its spectral mass ||P_c F||^2=2/35-d^2>0.

## 3. A paid long-time tail of the KERNEL, not of nonlinear solutions

Define A(k)=|M_F(k)|^2/[2pi nu(k)]. Near zero,

    A(k)=k^2/(64pi)+O(k^4).

Make the exact odd change k=u sqrt(u^2+4), so nu=2+u^2. Then
B(u)=A(k(u))k'(u) is even Schwartz and

    B(0)=0, B''(0)=1/(4pi).

Fourier inversion with a Gaussian regularizer gives

    integral B(u)e^(it u^2)du
      =e^(i pi/4)sqrt(pi/t)/(2pi)
                        integral Bhat(xi)e^(-i xi^2/(4t)) dxi.

Since integral Bhat=0 and integral xi^2 Bhat=-2pi B''(0), expand the last
exponential once. Its error is <=xi^4/(32t^2); the required Fourier moment
is finite because B is Schwartz. Taking the imaginary part yields

    k_c(t)=sin(2t+3pi/4)/(16sqrt(pi) t^(3/2))
                                      +O(t^(-5/2)).      (6)

This proves |k_c(t)|<=C(1+t)^(-3/2) for t>=1; near zero
|k_c(t)|<=t ||P_cF||^2. Hence k_c is absolutely integrable. Its discrete
shape counterpart is not; the full k is NOT in L1 and does not decay.
Moreover integral_0^infinity t|k_c(t)|dt diverges: the nonzero leading
oscillatory tail has divergent t^(-1/2) absolute integral, while the weighted
remainder is integrable. An unproved finite first-moment estimate cannot
be used to justify instantaneous slowly modulated friction.
Equation(6) is a linear response-kernel theorem. It is not an all-time
asymptotic theorem for the full nonlinear fields or mode amplitude.

## 4. Passivity without a fabricated bath

For Re s>0,

    ktilde(s)=<F,(H+s^2)^(-1)F>.

The limit ktilde(s)->1/70 as positive real s->0 is an Abel/static response.
It is not an ordinary integral of the full k, whose shape oscillation
prevents convergence of that integral.

The velocity-response transfer is s ktilde(s). At each spectral value
lambda>=3, for s=sigma+i tau, sigma>0,

    Re[s/(s^2+lambda)]
      =sigma(|s|^2+lambda)/|s^2+lambda|^2>0.              (7)

Equivalently solve v_tt+Hv=F r(t), zero data, as a characterization of
the linear response operator. Its work satisfies
integral_0^T r<F,v_t>dt=(||v_t||^2+<v,Hv>)/2>=0.
The actual source uses r=-6z^2; no independent drive is appended to FIN.
This is passivity/storage of that positive linear channel, not positivity
of every signed nonlinear sector or a source of conserved physical energy.

## 5. Exact finite-mode closure and noise are not supplied

Any finite-dimensional constant-coefficient linear realization has an
impulse response which is a finite sum of polynomial-times-exponentials
(and possible distributions at t=0). Subtracting the explicit shape
oscillation would leave such a response for k_c. Its nonzero oscillatory
t^(-3/2) tail in (6) is impossible for such a finite sum. Thus the actual
response cannot be represented EXACTLY by finitely many linear auxiliary
oscillators/constant-coefficient Markov states. A finite-time approximation
would need its own paid error. This is a scoped realization statement,
not a theorem excluding every nonlinear hidden source or asymptotic closure.

For arbitrary initial scalar-wave data the elimination also contains its
homogeneous solution projected against F. Under the ACTUAL RP data this
term is zero. With no supplied random initial measure, no independent
Gaussian forcing or FDT noise is derived by calling that missing term a
bath. A positive Gamma or passive kernel does not choose a stochastic law.
