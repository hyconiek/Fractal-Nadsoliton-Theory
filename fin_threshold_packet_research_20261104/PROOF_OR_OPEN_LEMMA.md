# TP: finite-energy threshold packets on the dispersive time scale

Date: 2026-11-04. Base: fb9464d6a91039052100bf9ec64b0136bdda0bdd.
Exact linear comparison family and controlled packet asymptotics for the
actual potential. NOT a new preparation for the nonlinear RP solution.
No formal verification, independent review or full cubic trace theorem.

## 1. Packets with the correct potential, parity and continuous projection

Use VD's normalized Darboux map Wcal, with Wcal*Wcal=1,
Wcal Wcal*=Pc and L Wcal=Wcal(-partial_x^2+m^2), m^2=4+w.
Let

    chi(y)=pi^(-1/4)exp(-y^2/2),
    h_R(x)=R^(-1/2)chi(x/R), R>=1,
    U_R(t)=Wcal exp(-it sqrt(-partial_x^2+m^2)) h_R.    (TP1)

This is a complex positive-frequency LINEAR solution for the exact
potential. Its real part is a real even solution with position Wcal h_R
and zero physical velocity at t=0. These are comparison data only: the
original nonlinear preparation has zero continuous initial field and is
not changed or approximated by inserting this packet.

Unitarity/intertwining give exactly

    ||U_R(0)||2=1,  Pc U_R(t)=U_R(t),
    ||sqrt(L)U_R(0)||2^2=m^2+1/(2R^2),
    ||(L-m^2)U_R(0)||2=R^-2 ||chi''||2.              (TP2)

The phase energy is finite. In distorted spectral variables k=q/R the
packet has Gaussian density proportional to exp(-q^2). For any delta>0,
the residual also gives the elementary spectral-tail bound

    ||1_(|k|>delta) U_R(0)||2^2
                        <=(R delta)^-4 ||chi''||2^2. (TP3)

This is a spectral projection estimate, not multiplication by a physical
spatial cutoff. No bound component is introduced by cutting off r0.

## 2. Frozen threshold coherence has a finite time scale

By |exp(-ita)-exp(-itb)|<=t|a-b| and
sqrt(m^2+k^2)-m<=k^2/(2m),

    ||U_R(t)-exp(-imt)U_R(0)||2
                          <=t R^-2 ||chi''||2/(2m).  (TP4)

Thus t=o(R^2) is sufficient for a vanishing global frozen-carrier error.
At the critical matching R=epsilon^-1, t=S epsilon^-2 is of order R^2;
small epsilon alone does not justify a frozen envelope throughout it.

This limitation is genuine for TP1, not merely a loose upper bound.
Define its normalized dephased autocorrelation at t=tau R^2:

    C_R(tau)=<U_R(0), exp(imt)U_R(t)>
       =pi^-1/2 integral exp(-q^2)
          exp(-i tau R^2[sqrt(m^2+q^2/R^2)-m])dq.

For fixed tau, dominated convergence yields

    C_R(tau) -> [1+i tau/(2m)]^(-1/2).                (TP5)

The principal square root is fixed by value one at tau=0. The limiting
modulus is less than one for nonzero tau, and the frozen-field distance
squared tends to 2-2 Re[1+i tau/(2m)]^-1/2>0. This is dispersive change
in a unitary linear solution, not loss of total energy or nonlinear damping.

## 3. A controlled local threshold envelope through t=O(R^2)

Write tau=t/R^2 in [0,S] and define

    A(tau)=pi^(-1/4)[1+i tau/m]^(-1/2).                (TP6)

For the exact-potential packet, uniformly in tau and the fixed mass band,

    U_R(t,x)=R^-1/2 exp(-imt)r0(x)A(tau)+e_R(t,x),
    ||e_R(t)||H1_-6<=C_S R^-5/2,
    ||U_R,t(t)+im R^-1/2 exp(-imt)r0 A(tau)||L2_-6
                                                   <=C_S R^-5/2. (TP7)

Here r0=(3tanh^2-1)/2. The remainder is local weighted, not global:
the displayed non-L2 leading profile is never subtracted in global energy.
The momentum estimate is derived spectrally, not by differentiating an
unspecified error bound.

Details of the uniform estimate: in the rescaled Fourier integral,

    0<=q^2/(2m)-R^2[sqrt(m^2+q^2/R^2)-m]
                                      <=q^4/(8m^3 R^2).

The resulting phase error is bounded by C_S q^4/R^2 against a Gaussian.
The exact Darboux multiplier on exp(ikx) is

    a(x,k)=[3K^2-1-k^2-3ikK]/sqrt((1+k^2)(4+k^2)).

Its even part differs from r0 by O(k^2), uniformly in x. Its odd part
is O(k); paired with an even Fourier amplitude, its integral includes
sin(qx/R), so |sin(qx/R)|<=|qx|/R gives O(|x|/R^2).
The leading free envelope is

    Phi_tau(y)=pi^(-1/4)(1+i tau/m)^-1/2
                      exp[-y^2/(2(1+i tau/m))].

It is even with uniformly bounded derivatives on [0,S], hence
|Phi_tau(x/R)-Phi_tau(0)|<=C_S x^2/R^2. Differentiating the actual
integral once in x has the same polynomial bound: derivatives of K are
bounded, while each differentiated Fourier factor gives q/R. Gaussian
moments control every resulting term. The phase-error multiplier and
its x derivative are bounded by C(1+|q|) times that phase error for R>=1.
Thus the pointwise error and its derivative are bounded by
C_S R^-5/2 <x>^2, which is integrable with the stated H1_-6 weight.
For velocity, insert -i sqrt(m^2+k^2) into the original integral; its
difference from -im costs O(k^2), proving the last line of TP7.

## 4. The actual compact observation of this comparison family

The Darboux observation in VD becomes exactly q of the free packet.
Let c_q=integral_-1^1 <x>^-1 dx=2arsinh(1). Then

    q(Wcal*U_R(t))
       =c_q R^-1/2 exp(-imt)A(t/R^2)+O_S(R^-5/2).     (TP8)

For a fixed S>0, squaring and changing time variable gives

    integral_0^(S R^2)|q(Wcal*U_R)|^2dt
       =c_q^2 R [m/sqrt(pi)]arsinh(S/m)+O_S(R^-1).    (TP9)

For the REAL packet Re U_R, the leading coefficient is one half:

    integral_0^(S R^2)|q(Wcal*Re U_R)|^2dt
       =(c_q^2/2)R [m/sqrt(pi)]arsinh(S/m)+O_S(R^-1). (TP10)

The oscillatory square term is bounded by integration by parts in tau;
its phase frequency is 2mR^2 and A(tau)^2 has bounded derivative on the
fixed interval. The amplitude error contributes O(R^-1) after integration.

Scaling the real packet by epsilon and choosing R=epsilon^-1 gives
an integrated squared trace of order epsilon, with a positive explicit
leading constant and O_S(epsilon^3) error. This is a sharp comparison
obstruction to obtaining O(epsilon^2) from parity, continuous projection,
small global energy and near-threshold concentration ALONE. It is NOT a
counterexample to the original prepared nonlinear trajectory, which has
no such initial radiation. No independent incoming packet is added there.

## 5. A finite-energy input interface for DT's quadratic cascade

Let P_R=Re U_R be this comparison even packet. TP7 gives, in the same
negative-weight H1 norm,

    K P_R(t)^2=R^-1 K r0^2
       [A(tau)^2 exp(-2imt)+2|A(tau)|^2
                          +conj(A(tau))^2 exp(2imt)]/4
                                      +O_(H1_-6,S)(R^-3). (TP11)

For the product estimate, the exact packet and its x derivative are
uniformly O_S(R^-1/2) by their Gaussian Fourier integrals. Multiply TP7's
weighted error by these bounded factors; the leading profile has the same
bounds. This yields R^-3 without assuming a global non-L2 remainder.

Multiplication by the source factor -6 identifies DT's static coefficient
-3R^-1|A(tau)|^2 and second-harmonic coefficient
-(3/2)R^-1 A(tau)^2. Their slow envelopes are now explicit and vary on
R^2 time. DT's static and outgoing operator responses are therefore the
correct local reference profiles to investigate, but not yet a proved
retarded approximation on that whole interval.

The restriction is essential: TP11's error is controlled with NEGATIVE
spatial weight. LE's Duhamel estimate requires a positive-weight input.
One cannot feed TP11 into LE and declare the accumulated cascade error
small. The far spatial part of the finite packet and the outgoing harmonic
history must be estimated in a norm that the propagator actually accepts.

## 6. Next exact obligation and unchanged actual status

Derive a retarded quadratic-cascade estimate for these finite packets,
with spatial tails and the slowly varying harmonic amplitudes retained.
Then prove an appropriate link to the actual generated continuum; TP1
is a test family, not that link. The width-time relation and the outgoing
intermediate harmonic in DT must remain explicit.

ML's actual whole local-continuum O(epsilon) bound is unchanged. The
O(epsilon^2) cubic trace and all modal action errors remain open. The
current result concerns the actual LINEAR operator with comparison data,
not a new nonlinear preparation, empirical law or natural FIN source.
