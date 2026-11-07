# PA: weighted propagation and causal-response audit

Date: 2026-11-07. Parent: 7d8864a21b35836aaaf6a40534458dc21726a71d.
Same-author analytical audit. ACCEPT_AT_STATED_DEPENDENCY_BOUNDARY.
Independent review and formal verification remain open.

## 1. Scope

This audit reconstructs LE's weighted propagator expansion, QH's mixed
carrier gap and local primitive, and their use in RW4--RW6. It inherits
the exact Darboux spectral resolution, SP smallness and MA's canonical
envelope bounds. It does not certify CM's construction or CD's entire
spectral-density calculation. No modal damping conclusion is used as an
input. All constants are uniform on the original compact mass band;
constants involving a critical window may depend on fixed S.

## 2. Threshold matrix and forcing regularity

The continuous spectral kernel uses measure dk/(2pi), frequency
nu_m(k)=sqrt(m^2+k^2), and a Jost product whose value at k=0 is
r0(x)r0(y), where r0=(3tanh^2(x)-1)/2. With a fixed even low cutoff,
put k=u sqrt(2m+u^2). Then nu_m=m+u^2 and dk/du at zero is sqrt(2m).
The stationary-phase constant is therefore sqrt(m/(2pi t)). Combining
the cosine, sine/nu and -nu sine blocks gives, for t>=1,

    U_m(t)(a,b) = r0/sqrt(2pi m t)
        (m cos(theta) mu_a+sin(theta) mu_b;
         -m^2 sin(theta) mu_a+m cos(theta) mu_b)+R_m(t)(a,b),
    theta=mt+pi/4, mu_h=integral r0 h.                 (PA1)

The momentum factor and phase shift are fixed by the actual propagator.
They are not obtained by differentiating an unquantified asymptotic error.
The inner matrix has determinant m^2 and is uniformly invertible on the
mass band. Odd input has both moments zero. Even input need not.

For each low-frequency phase block, including one position derivative,
the transformed amplitude has three u derivatives bounded in L2 by
C<x-y>^3. The Gaussian Fourier identity bounds the remainder by
Ct^-3/2(||B''||2+||B'''||2). Explicitly, the Fourier integral of
xi^2|Bhat(xi)| is bounded by C times those two norms, by Cauchy-Schwarz
with (1+xi^2)^-1. The weights <x>^-6<y>^-6 make the resulting kernel
Hilbert-Schmidt. Differentiating the weighted output adds a bounded
weight derivative. Thus no hidden input derivative beyond X_6 is needed
in this low-frequency step.

For high frequencies |k|>=8 and |x-y|<=t/2, the phase derivative is
bounded below by 1/4 since 8/sqrt(64+23/4)>3/4. The largest symbol order
is one, including the momentum-from-position block. Eight integrations
by parts reduce it to an integrable symbol. A temporary compact frequency
cutoff has transition errors O(delta^(8-M-1)), M<=1, tending to zero.
This proves the smooth middle-cone kernel estimate without discarding a
wavefront distribution.

To cover the rest, split input and output by smooth spatial cutoffs at
scale t. An exterior input costs t^-6 through X_6; an exterior output
costs t^-6 through X_-6. The global energy bound of the continuous
high-frequency propagator pays each such term. Cutoff derivatives cost
O(1/t) and are controlled at H1 regularity. This separates the wavefront
from the middle-cone calculation. Small times use global energy directly.

Consequently LE's X_6 to X_-6 decay and remainder follow from the
inherited spectral resolution: t^-1/2 in general and t^-3/2 for odd
input. It is local weighted decay, not decay of global radiation energy.

## 3. A nonzero threshold coefficient for the actual mixed profile

For f=sqrt(3/4)sech^2 x, g=sqrt(3/2)sech x tanh x and
B=P_c(Kfg), the projected threshold moment is exactly

    mu_B = 3pi/(64sqrt(2)) != 0.                      (PA2)

Indeed integral r0 f=0 and integral r0 g=0, so projection does not alter
this moment. Setting z=tanh x reduces the unprojected integral to

    (3/(2sqrt(2))) integral_(-1)^1
                       z^2(3z^2-1)sqrt(1-z^2)/2 dz.

The two elementary moments are pi/8 and pi/16, yielding PA2. For initial
phase (0,B), PA1 is therefore a genuine nonzero leading phase vector.
This rules out deleting the even threshold term for this profile.
The exact coefficient supplies a reproducible normalization control;
it is not a new nonlinear damping coefficient.

## 4. Uniform quadratic carrier gaps and their limits

For lambda=sqrt(3), omega=sqrt(w), m=sqrt(4+w), the mixed carriers are
nu=+/-lambda+/-omega. The exact squared gaps are

    (lambda+omega)^2-m^2=-1+2sqrt(3w)>=sqrt(15)-1>2,
    m^2-(lambda-omega)^2=1+2sqrt(3w)>=1+sqrt(15)>4.

Since lambda+omega+m<6 uniformly, all distances |nu +/- m| exceed 1/3.
This is a gap from threshold, not a gap from the continuous spectrum:
the sum carrier is open, and its outgoing stationary response is not a
global L2 vector. The difference carrier is closed. Negative carriers
are retained by conjugation. The cubic carrier 2lambda-omega touching
threshold at w=4/3 is outside this quadratic inventory and is not removed
by this argument. No parameter value is excluded or fitted.

## 5. Local primitive, susceptibility sign and causal endpoint

For K(t)=U_m(t)(0,B), PA1 writes the large-time kernel as a sum of
fixed phase vectors times t^-1/2 exp(+/-imt), plus an X_-6 remainder
O(t^-3/2). Multiplication by exp(i nu t) and the preceding gap give

    |integral_t^infinity s^-1/2 exp(i delta s)ds|
       <=2|delta|^-1 t^-1/2, |delta|>1/3.

Only the explicit leading term is integrated by parts. The remainder
is integrated absolutely; its time derivative is not assumed. Together
with the bounded short-time kernel this proves in the local phase space

    P_nu=integral_0^infinity exp(i nu s)K(s)ds,
    E_nu(t)=integral_t^infinity exp(i nu s)K(s)ds,
    ||P_nu||+ (1+t)^1/2||E_nu(t)|| <= C.               (PA3)

At a threshold carrier the leading integral would instead contain a
nonoscillatory s^-1/2 component, and PA2 shows why it would fail. The
actual gap is essential. Continuous projection is also essential:
bound-mode terms do not have LE's decaying remainder.

With an Abel factor exp(-a t), a>0, the scalar sine/nu transform is
1/[Omega^2-(nu+ia)^2], and the cosine transform is
-i(nu+ia)/[Omega^2-(nu+ia)^2]. Thus the causal position susceptibility
in RW has the stated +i0 convention. Passage to the local improper
limit is justified by PA3. It is not a claim of convergence in global
energy norm or a new calculation of CD's densities.

For every C1 envelope A, finite-time integration by parts gives exactly

    integral_0^t K(t-s)exp(-i nu s)A(s)ds
      =exp(-i nu t)[P_nu A(t)-E_nu(t)A(0)
                       -integral_0^t E_nu(r)A'(t-r)dr]. (PA4)

The endpoint E_nu(0)=P_nu enforces zero initial phase. QH's physical
mixed envelope has A(0)=0, but RW uses canonical envelopes and includes
self channels with nonzero preparation. Its initial term must be kept.
An exact scalar example in the checker tests both this term and the
history sign; omitting the transient gives a nonzero initial response.

## 6. RW work budget and audited dependency boundary

For odd Fa,Fb, parity makes the propagator absolutely integrable locally;
its tail is O(t^-1/2), including the DC carrier. For even G, PA3 gives
the same bound by modulation. These are precisely RW's required profiles.

MA supplies sup|A|<=C epsilon^2 and ||A'||L2(0,T)<=C_S epsilon^3.
The initial transient has L2 norm at most
C epsilon^2 sqrt(log(2+T)). The convolution history has L2 norm at most
C sqrt(1+T)||A'||L2 by Young. At T=S epsilon^-2 this is RW5.
The actual work weight has L2 norm O_S(epsilon), so the resulting work
error is O_S(epsilon^3 sqrt(1+log(1/epsilon))). Every bound holds uniformly
in the terminal time; no constant envelope or pointwise slow modulation
was assumed. The absolute remainder at times near zero is paid by the
bounded local kernel, not a singular asymptotic approximation.

No blocking defect was found in these inspected propagation and causal
transfer steps. The result remains conditional on the exact spectral
resolution and the stated global/envelope bounds. It is a same-author
audit, not independent review of the complete action law.

Next audit CM's actual finite-rank canonical map, Hamiltonian remainder
regularity, physical reconstruction and preparation shifts. Those are
central inputs to MA's work identities and were not established by the
propagation calculation. All physical-source and empirical gates remain
unchanged; no all-time scattering or waveform theorem is inferred.
