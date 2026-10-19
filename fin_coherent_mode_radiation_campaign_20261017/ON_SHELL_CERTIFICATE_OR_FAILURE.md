# RP-003 — exact strictly positive scalar continuum spectral density

**Disposition: CERTIFIED_CONDITIONAL_ON_SHELL_COUPLING** for the predeclared F_RP family, with no retuning following the calculation. It is a spectral/source theorem, NOT a nonlinear radiation or spatial energy-transport theorem.

## Vertex and source normalization

The existing source-compatible vertex is V_int=phi^2 sum_j y_j |Psi_j|^2. For the neutral data phi(0)=K, phi_t(0)=0, Psi(0)=epsilon f, Psi_t(0)=0, the leading linear mode is Psi^(1)=epsilon f cos(omega t). The exact quadratic Phi equation obtained from -V_phi has forcing

    -2 K sum_j y_j |Psi_j^(1)|^2 = -epsilon^2 S(x)[1+cos(2 omega t)],
    S(x)=K f^T Y f=(9/4) tanh x sech^4 x.

This factor and the DC term are compulsory. No free continuum wave was added. In particular the normalized 2omega cosine coefficient is **-epsilon^2 S**, not -2epsilon^2 S. U is quartic, no phi6.

The scalar linear operator is L_phi=-d_x^2+4-6 sech^2 x, with eigenvalues 0 (even translational K') and 3 (odd shape state) and continuum [4,infinity). The bound-mode harmonic satisfies Omega=2omega, Omega^2=4(a+2) in [5,7], strictly >4, away from threshold by at least 1. Let q=sqrt(Omega^2-4)=2sqrt(a+1) in [1,sqrt(3)]. The static S is odd: it is orthogonal to even translation, but not automatically to the odd shape mode. Define P_c by removing BOTH scalar bound states.

## Exact on-shell transform

Use scattering states normalized to unit outgoing amplitude at +infinity (the reflectionless l=2 Poeschel–Teller problem), with the full q in R continuum normalization dq/(2pi):

    psi_q(x) = exp(i q x) [3 tanh^2 x -3 i q tanh x -(q^2+1)] /
                          [(1-i q)(2-i q)],   q real.

Direct differentiation gives (L_phi-(q^2+4))psi_q=0. The scattering map is unitary with this normalization: psi_q asymptotes to unit-modulus plane waves at both ends, with no reflection, and dE_c=|psi_q><psi_q| dq/(2pi). Bound states are orthogonal to psi_q; hence <psi_q,P_c S>=<psi_q,S>.

Let F4(q)=int_R sech^4 x exp(-i q x) dx= pi q(q^2+4)/[6 sinh(pi q/2)] (positive for q>0). Integration by parts and the Fourier recurrence F6/F4=(q^2+16)/20 give, with C=9/4,

    M(q)=int_R conj(psi_q(x)) S(x) dx
        = C * i q(q^2+4)F4(q) / [8(1+i q)(2+i q)].

For q>0 none of the factors vanish. This is more than P_c S!=0: it is a nonzero **at the actual frequency** on-shell transition amplitude.

Apply Sokhotski–Plemelj / limiting absorption for exponentially localized S in the explicit reflectionless operator, and the stated continuum spectral resolution. The absorption sign is positive for (L_phi-Omega^2-i eta)^(-1), eta>0. Both scattering roots +q and -q contribute:

    Gamma=lim_(eta -> 0+) Im <P_c S,(L_phi-Omega^2-i eta)^(-1)P_c S>
         = |M(q)|^2/(2q)
         = 81 q(q^2+4) F4(q)^2 / [2048 (1+q^2)]  > 0.

The strict positivity holds for **every** a in [-3/4,-1/4], not a point-picked resonance. Compactness and continuity imply uniform positive infimum over this closed interval. For the diagnostic center a=-1/2, q=sqrt(2) and Gamma ~0.106376851619 in normalized model units. These numbers do not select the FIN source coefficients or measure experimentally emitted power.

## What the certificate does not yield

Gamma is an infinite-time, retarded spectral response coefficient for the *linearized forced scalar channel*. It does not by itself say that the actual nonlinear neutral Cauchy solution emits a positive integrated canonical outward energy flux through +/-R on a uniformly controlled, declared time scale. The actual source is switched on at t=0, includes its DC/static and shape-mode terms, and generates an initial slip; the scalar and transverse errors must be compared on the same trajectory. RP-004 is not earned by this file.
