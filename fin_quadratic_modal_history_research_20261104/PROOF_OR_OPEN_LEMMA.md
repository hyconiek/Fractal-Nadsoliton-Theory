# QH: the actual quadratic modal history on critical time

Date: 2026-11-04. Base: 370dba36f68919521016cb645e19dab166baed49.
Conditional analytical component result in the unchanged RP model.
Not formal verification, independent review, or a full damping theorem.

## 1. Dependency delta and exact amplitudes

BT paid the bound-only cubic/quintic response. Its BT2 residual retained
the quadratic source -12cz B, B=Pc,e(Kgf), together with every source
term containing continuum fields. Here we pay that quadratic response
using ACTUAL modal histories. This extends the frozen CD response, not
its spectral-density precision. SP's energy bound and LE's weighted
phase expansion are the analytical inputs.

Keep physical eta=cg+v, psi=zf+u, c=<g,eta>, z=<f,psi>, with the same
odd/even sector, preparation and coefficients. Let lambda=sqrt(3),
omega=sqrt(w), m=sqrt(4+w), w in [5/4,7/4]. The projected equations are

    c''+lambda^2 c=F_c,
    F_c=<g,-6K(eta^2+psi^2)-6eta psi^2-2eta^3>,
    z''+omega^2 z=F_z,
    F_z=<f,-12K eta psi-6eta^2 psi-g4 psi^3-g6 psi^5>.   (QH1)

SP gives global energy size C epsilon for the actual small solution.
In one dimension ||ab||_2<=||a||_infinity ||b||_2 and H1 embeds in
L-infinity. Hence |F_c|+|F_z|<=C epsilon^2 for all finite times,
uniformly in the parameter band and fixed finite signed g4,g6 bounds.
The same energy bound controls c,c',z,z' by C epsilon.

Define exact complex interaction amplitudes, not prescribed envelopes:

    a(t)=exp(i lambda t)(c+i c'/lambda),
    b(t)=exp(i omega t)(z+i z'/omega).

Their derivatives are exactly

    a'=i exp(i lambda t)F_c/lambda,
    b'=i exp(i omega t)F_z/omega.                       (QH2)

Thus |a|+|b|<=C epsilon and |a'|+|b'|<=C epsilon^2.
No second derivative, adiabatic assumption, or c=O(epsilon^2) hierarchy
has been used. Physical preparation gives a(0)=0, b(0)=epsilon.
These are physical oscillator coordinates, not CM or QP momenta.

Put a_+=a, a_-=conj(a), and similarly b_tau. Then EXACTLY

    cz=(1/4) sum_(sigma,tau=+/-1)
          A_sigma,tau(t) exp(-i nu_sigma,tau t),
    A_sigma,tau=a_sigma b_tau,
    nu_sigma,tau=sigma lambda+tau omega.               (QH3)

Each A starts at zero, is O(epsilon^2), and has derivative O(epsilon^3).
These envelopes may contain nonlinear oscillations. No spectral support
restriction on them is inferred from their derivative bound.

## 2. A uniform gap for these quadratic carriers only

The positive sum and difference obey

    (lambda+omega)^2-m^2=-1+2sqrt(3w)>=sqrt(15)-1>2,
    m^2-(lambda-omega)^2=1+2sqrt(3w)>=1+sqrt(15)>4.     (QH4)

Also lambda+omega+m<6 on the whole band. Factoring the squared gaps
therefore gives |nu_sigma,tau +/- m|>1/3 for every branch. This is
separation from the THRESHOLD, not absence of continuum resonance:
lambda+omega lies in the continuous spectrum and its outgoing response
is not a global L2 state. CD's mixed radiative channel remains present.

BT's cubic carrier 2lambda-omega still reaches m at w=4/3. QH4 neither
contradicts nor removes it. We use different quadratic carriers here.

## 3. A convergent local primitive and its tail

Work in LE's X_-6=H1_-6 x L2_-6 phase norm. For the fixed exponentially
localized continuous profile B define

    K(t)=U_m(t)(0,B)=(S_L(t)B, cos(t sqrt(L))B),
    P_nu=integral_0^infinity exp(i nu t) K(t) dt,
    E_nu(t)=integral_t^infinity exp(i nu s) K(s) ds.    (QH5)

Both integrals are improper limits in the weighted LOCAL phase space.
LE's expansion gives K(t) as a linear combination of
t^-1/2 exp(+/- i m t) times fixed weighted phase vectors, plus an
X_-6 remainder O(t^-3/2). The short-time kernel is bounded.
For any real delta with |delta|>1/3, integration by parts gives

    |integral_t^infinity s^-1/2 exp(i delta s)ds|
                         <=2 |delta|^-1 t^-1/2, t>=1.

The absolute integral of the remainder is O(t^-1/2). Consequently,
uniformly for all four carriers and the entire w band,

    ||P_nu||_X_-6<=C,
    ||E_nu(t)||_X_-6<=C(1+t)^-1/2.                    (QH6)

This proves the primitive and its tail directly; no unproved uniform
global resolvent bound is imported. Negative carriers follow by complex
conjugation. The momentum component uses LE's cosine expansion, not
differentiation of an unspecified error term. The continuous projection
in B is essential: retained bound modes would invalidate QH5.

## 4. Exact retarded integration by parts, with the initial transient

For a C1 scalar envelope A and carrier nu the zero-data phase response is

    Y_nu,A(t)=integral_0^t K(t-s)exp(-i nu s)A(s)ds.

Since the finite primitive equals P_nu-E_nu(t), integration by parts
on the finite interval gives the exact LOCAL identity

    Y_nu,A(t)=exp(-i nu t)[P_nu A(t)-E_nu(t)A(0)
                         -integral_0^t E_nu(r)A'(t-r)dr]. (QH7)

At t=0, E_nu(0)=P_nu, so it has exactly zero initial phase. The formula
does not replace the actual Cauchy state by a steady outgoing wave.
Its last two terms retain switch-on and all envelope history. P_nu is
an intermediate local vector; no global energy is assigned to it.

For the actual envelopes in QH3, A(0)=0. Equations QH2 and QH6 yield

    ||Y_nu,A(t)||_X_-6
                  <=C[epsilon^2+epsilon^3 sqrt(1+t)].  (QH8)

The second term follows from integral_0^t(1+r)^-1/2 dr<=2sqrt(1+t).
Let u_Q be the exact zero-data response to -12cz B. Then
(u_Q,u_Q,t)=-3 sum Y_nu,A, with all four branches retained, and for
fixed S, T=S epsilon^-2,

    sup_(0<=t<=T)||(u_Q,u_Q,t)||_X_-6<=C_S epsilon^2,
    integral_0^T ||(u_Q,u_Q,t)||_X_-6^2dt<=C_S epsilon^2. (QH9)

Indeed the integral is bounded by C[epsilon^4 T+
epsilon^6(T+T^2)]. The amplitude has only the already proved O(epsilon)
bound. The improvement over absolute Duhamel comes from the oscillatory
primitive and its decaying tail, not from an assumed small shape mode.

For comparison, keeping only ||K(t)||<=C(1+t)^-1/2 and |cz|<=Cepsilon^2
would give O(epsilon) local size and O(1) integrated square at T.
Those crude bounds do not close the gate. QH7 supplies the missing gain.

## 5. The exact remaining source contains continuum fields in every term

Use BT's u_D without modification, and define u_C=u-u_D-u_Q. It has
zero physical Cauchy data and the exact source

    R_C=Pc,e[-12K(cg u+zf v+vu)
       -6(c^2 g^2 u+2cz gf v+2cg vu+z f v^2+v^2u)
       -g4(3z^2 f^2u+3z f u^2+u^3)
       -g6 sum_(j=1)^5 binom(5,j)(zf)^(5-j)u^j].       (QH10)

Here u,v are the WHOLE actual physical continuum fields, not u_C or
their leading approximations. This is an exact source decomposition,
not a closed equation with omitted feedback. Every monomial in R_C
contains at least one continuum factor. The coefficient 2cz in the
mixed cubic term is retained, and signed g4,g6 are unchanged.

Let q_D,q_Q,q_C be the compact observation after VD's normalized
Darboux map of the respective physical positions. BT's inverse-QP
correction delta_q has pointwise size Cepsilon^2. Consequently

    q_psi=q_D+q_Q+q_C+delta_q,
    integral_0^T |q_psi|^2dt
                       <=4 integral_0^T |q_C|^2dt+C_S epsilon^2,
    integral_0^T E_rad^2dt
                       <=C_S epsilon^2+C integral_0^T |q_C|^2dt. (QH11)

Weighted Darboux boundedness transfers QH9 to q_Q. Position corrections
only are used; QP momentum is never equated to a position derivative.

## 6. The next gate, stated without a circular bootstrap

The remaining trace is now driven entirely by continuum-bearing source
terms. Their presence does NOT make their response small automatically.
For example Kvu is not localized by K, and global H1 smallness does not
pay a positive spatial moment of vu. The gapped energy estimate alone
cannot be inserted into LE as a weighted forcing bound. Terms with f or
g do localize but depend on unknown actual local continuum norms.

The next step is to combine this exact source reduction with QP's
quadratic far-field transformation and VD's local energy identity, paying
the chart corrections and a noncircular estimate for q_C. A bound which
simply assumes integral |q_C|^2=O(epsilon^2), or has an unabsorbable
coefficient depending on epsilon^-1, does not close it. Retain all three
CD action channels; their actual accumulated error is a separate gate.

No full critical damping, scattering, nonlinear instability, physical
radiation law, natural action/clock/geometry, F01, strict/legacy or
selector closure follows. No physical premise or new preparation is added.
