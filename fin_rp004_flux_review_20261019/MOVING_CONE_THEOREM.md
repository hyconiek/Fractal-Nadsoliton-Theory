# Retarded limiting amplitude on the original moving cone

Date: 2026-10-19. **New analytical work in this review**, not a claim made
by the incoming A report. Conditional on the inherited, normalized Darboux
resolution in the RP intake. This is not a mechanized or peer-reviewed proof.

## 1. Statement and uniformity

Write `L=-d_x^2+4-6sech^2 x`, `S=(9/4)tanh x sech^4 x` and

    Omega=2sqrt(a+2), q=sqrt(Omega^2-4), a in[-3/4,-1/4].

Thus `q in[1,sqrt(3)]`. Let **u be the zero-data solution**

    u_tt+Lu=-S(1+cos(Omega t)).                         (1)

Let `M(k)=<psi_k,S>` in the inherited `dk/(2pi)` normalization, with

    psi_k(x)=exp(ikx)[3tanh^2 x-3ik tanh x-(k^2+1)]/D(k),
    D(k)=(1-ik)(2-ik),
    M(k)=(9/4)ik(k^2+4)F4(k)/[8(1+ik)(2+ik)],
    F4(k)=pi k(k^2+4)/[6sinh(pi k/2)].

All removable values at zero are understood continuously. The functions M,
and every polynomial multiple needed below, are Schwartz on the real axis.
Define `C(a)=-i M(q)/(2q)` and `w=Re[C(a)exp(iqx-iOmega t)]`.
There is a finite constant C, independent of a, t and x, such that for t>=1
and `t/6<=x<=t/3`,

    |u-w|+|u_t-w_t|+|u_x-w_x| <= C t^(-1/2).          (2)

The negative band follows by the exact oddness of u. The derivatives in (2)
are partial derivatives at fixed x or t. They are **not** derivatives along
the ray, nor derivatives of an asymptotic statement with x/t held fixed.

## 2. Two elementary oscillatory-integral estimates

These estimates are proved here to expose the hypotheses actually used.
Set `nu(k)=sqrt(k^2+4)` and `v=x/t in[1/6,1/3]`.

**Smooth amplitudes.** For a uniformly Schwartz family A(k),

    integral exp(it[kv +/- nu(k)]) A(k) dk = O(t^-1/2). (3)

The stationary point is `k_s=mp 2v/sqrt(1-v^2)`, in a fixed compact interval
with `|k_s|<=1/sqrt(2)`. On a fixed neighborhood its second derivative has
absolute value `nu''=4/nu^3` bounded below. On `|k-k_s|<=t^-1/2` the integral
is bounded by its length times sup|A|. On the remaining compact neighborhood
integrate once by parts with `(it phi')^-1 d_k`. Since
`|phi'|>=c|k-k_s|`, the endpoint terms are O(t^-1/2), the A' term is
O(t^-1 log t), and the derivative of 1/phi' contributes O(t^-1/2).
Outside that neighborhood phi' is bounded away from zero, including at
infinity because v<=1/3<1; Schwartz bounds pay the integrations and endpoints.
This proves (3), with uniform constants. The same argument without a
stationary neighborhood gives O(t^-1) for a smooth nonstationary amplitude.

**A nonstationary simple pole.** For a compact smooth amplitude A near k0,
suppose phi'(k0) has fixed sign and is bounded away from zero on the support.
With integration normalized by 1/(2pi),

    integral exp(it phi(k)) A(k)/(k-k0-i0) dk/(2pi)
      = i A(k0) exp(it phi(k0)) 1_{phi'(k0)>0}+O(t^-1),
    integral exp(it phi(k)) A(k)/(k-k0+i0) dk/(2pi)
      =-i A(k0) exp(it phi(k0)) 1_{phi'(k0)<0}+O(t^-1). (4)

Here `i0` is a distributional boundary value, not a finite numerical
regularization. To prove (4), set
`r=(phi(k)-phi(k0))/phi'(k0)` on a small fixed support. This is an
orientation-preserving diffeomorphism with derivative 1 at k0. After a
smooth cutoff equal to one at zero, subtract `A(k0)/(r mp i0)`; the
remaining amplitude is smooth and compact, hence O(t^-1). The identities

    integral exp(i b r)/(r-i0) dr/(2pi)= i 1_{b>0},
    integral exp(i b r)/(r+i0) dr/(2pi)=-i 1_{b<0}

follow first for `r mp i eta`, eta>0, by closing the contour in the upper
or lower half-plane, then taking eta down to zero. Replacing the full-line
singular distribution by the cutoff changes the answer by O(1/|b|):
outside a neighborhood of zero the derivative of the discarded 1/r tail
is integrable, and its endpoint at infinity vanishes. These steps also
prove uniformity of (4) for the compact parameter ranges used here.

## 3. Exact causal response before splitting poles

For clarity, solve the **complex** harmonic problem with forcing
`-S exp(-iOmega t)` and zero data. In each spectral oscillator its multiplier
is exactly

    z_nu(t)=-[exp(-iOmega t)-cos(nu t)
                         +i(Omega/nu)sin(nu t)]/(nu^2-Omega^2). (5)

The numerator vanishes at nu=Omega, also after its first t derivative;
(5) is a continuous finite-time multiplier and has zero value/derivative
at t=0. Its real part is the harmonic contribution to u. Consequently the
denominator in the split below can consistently be read as
`nu^2-Omega^2-i0`: multiplying the numerator by its added delta term gives
zero. This prescription changes **no Cauchy data or forcing**.

Replace psi_k by exp(ikx) for the moment. The continuum part of (5) is

    - exp(-iOmega t) integral exp(ikx) M(k)/(k^2-q^2-i0) dk/(2pi)
    + integral exp(i[kx+nu t]) M(k)/[2nu(nu+Omega)] dk/(2pi)
    + integral exp(i[kx-nu t])
                M(k)(nu+Omega)/[2nu(k^2-q^2-i0)] dk/(2pi).      (6)

Equation (6) is an identity of distributions combined into the smooth
finite-time multiplier (5), not an assumption of an eternal outgoing wave.
The middle term is O(t^-1/2) by (3).

Use the exact distribution identity

    1/(k^2-q^2-i0)=(1/(2q))
                  [1/(k-q-i0)-1/(k+q+i0)].                   (7)

Choose fixed-radius smooth cutoffs at +q and -q, disjoint from the stationary
points. This is possible uniformly: q>=1, while |k_s|<=1/sqrt(2).
For the first integral in (6), the phase is kv. At +q its derivative v>0
selects the `-i0` residue in (4), while the `+i0` pole at -q contributes
no leading term. The smooth complement is O(t^-1). This gives precisely

    -i M(q)/(2q) exp(iqx-iOmega t).                           (8)

For the last integral in (6), the phase is `kv-nu(k)`. At +q,

    phi'=v-q/Omega <=1/3-1/sqrt(5)<0;

the `-i0` residue is therefore absent. At -q, `phi'=v+q/Omega>0`,
so the `+i0` residue is absent as well. Both pole neighborhoods are
O(t^-1) by (4). Their complement has uniformly Schwartz amplitudes and
is O(t^-1/2) by (3). The strict margin `1/sqrt(5)-1/3` is important:
this proof does not assert a sharp front for every frequency.

Applying partial t or x differentiation directly to (6) multiplies the
amplitudes by `-iOmega`, `+/-i nu`, or `ik`. These are still in the same
uniform classes. At the selected pole they give exactly the derivatives
of (8); at the unselected poles the sign tests are unchanged. Thus (8)
and its first derivatives have the stated uniform errors.

## 4. Pay the Jost correction, DC, discrete mode and threshold

For x>=t/6>0 the exact Jost correction is

    psi_k(x)-exp(ikx)
     =exp(ikx)[3(tanh^2 x-1)-3ik(tanh x-1)]/D(k).              (9)

Its value and first x derivative are bounded by `C exp(-2x)` uniformly
in real k: the denominator grows quadratically and the differentiated
numerator has degree at most two. Use the original, unsplit oscillator
formula or its Duhamel integral here. Since nu>=2,
`|z_nu|+|partial_t z_nu|<=C(1+t)` uniformly in nu and Omega. M is integrable,
so the correction in all three quantities in (2) is at most
`C(1+t)exp(-t/3)`. There is no multiplication of unrelated pole bounds.

The DC component is `-(1-cos(t sqrt L))L^-1 S`. Its plane continuum
static multiplier `-M(k)/(k^2+4)` is Schwartz, so its Fourier transform
and first x derivative decay faster than any fixed inverse power at
x~t; its time derivative is zero. Its cosine/sine transients and their
first derivatives satisfy (3). Formula (9) pays their Jost corrections.

The odd shape eigenfunction is `g=sqrt(3/2)sech x tanh x`, eigenvalue 3.
Both the DC and harmonic coefficients and their first time derivatives
are uniformly bounded, since `Omega^2-3>=2`. The mode and its first x
derivative are O(exp(-x)); their contribution in the band is exponentially
small. It is **retained**, not projected out. The even translation mode
has zero source coefficient by actual odd parity. The continuum threshold
k=0 causes no omitted endpoint: integration is on the whole real line,
nu(k) is smooth there, and M(k) has a smooth removable value. The stationary
points in (3) are not threshold endpoints.

Combining (3)--(9) proves (2) including the switch-on, both spectral roots,
DC and shape mode, with all constants uniform on the stated a interval.

## 5. A genuinely spatial energy lower bound

Let `B_T={T/6<|x|<T/3}` and use the actual scalar quadratic energy there:

    E_B(T)=1/2 integral_B [u_t^2+u_x^2+(4-6sech^2 x)u^2] dx.

For `w=|C|cos(qx-Omega T+arg C)`, its constant-potential density is

    e_w=Omega^2|C|^2/2
           -q^2|C|^2 cos(2qx-2Omega T+2arg C)/2.

The oscillatory integral on either band is uniformly O(1), because q>=1.
The pointwise error (2), integrated over total length T/3, changes energy
by O(sqrt T+1); the potential correction is exponentially small. Oddness
doubles the positive-side band. Since `|C|^2=Gamma/(2q)`,

    E_B(T)/T = Omega^2 Gamma/(12q)+O(T^-1/2).                 (10)

This proves the incoming report's previously conditional band calculation.
Put `gamma0=125/73728` and `c0=5 gamma0/21=625/1548288`. Since Omega^2>=5,
q<=sqrt(3)<7/4 and the inherited Gamma>gamma0, the coefficient in (10)
exceeds c0. For T>=6, the exterior density is positive at |x|>T/6,
so `E_ext(T;T/6)>=E_B(T)>=c0 T-C sqrt T-C`. This is **A-GATE**, now paid.

Finally the actual forced energy equation is
`e_t+(-u_tu_x)_x=-S(1+cos Omega t)u_t`. Integration on the exterior,
with radius held fixed during the time integral, gives

    F_u(R,T)=E_ext(T;R)
                  +integral_0^T integral_|x|>R S(1+cos Omega t)u_t dx dt.

The work term has absolute value at most `C T^2 exp(-4R)` by the global
energy estimate and the exponential S tail. Thus on the original ray

    F_u(T/6,T)>=c0 T-C sqrt T-C-C T^2 exp(-2T/3).             (11)

There exists a **uniform finite**, not numerically calibrated T0 such that
`F_u(T/6,T)>=(c0/4)T` for every T>=T0. Radius T/6 is the final radius of
each experiment, not a moving integration boundary R(t). No boundary
velocity term is omitted. This is a conditional classical reference theorem,
not yet the nonlinear theorem or a physical radiation claim.
