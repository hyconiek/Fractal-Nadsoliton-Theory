# RW - retarded quadratic work transfers to the three-channel polynomial

Date 2026-11-07. Parent 22be29e26926cefbe17be3e4cf58cd7d7727a0b9.
Self-audited conditional analysis using CM, CD, QH, TP2 and MA.

## 1. Actual canonical histories drive specified comparison responses

Use the same actual canonical variables (Z,P,C,R), actions I,J and
interaction amplitudes a,b as in MA. Define zero-data continuous responses

    v_Q,tt+H v_Q=-6 P_o(Fa Z^2+Fb C^2),
    u_Q,tt+L u_Q=-12 P_e(G ZC),                       (RW1)

where Fa=K f^2, Fb=K g^2 and G=K f g. All coefficients are evaluated on
the actual canonical trajectory. These profiles are exactly the quadratic
continuum forces from MA3. Higher terms of the transformed Hamiltonian
are retained in the difference from the whole continuum in section 5.
The initial comparison fields and velocities are zero.

Define W_I^Q,W_J^Q by MA9's work integrals with v_Q,u_Q in place of v,u.
We prove, uniformly for 0<=t<=T=S epsilon^-2,

    W_I^Q(t)=-integral_0^t (k_a I^2+k_x I J)ds+e_I^Q(t),
    W_J^Q(t)=-integral_0^t (k_b J^2+k_x I J)ds+e_J^Q(t),
    |e_I^Q(t)|+|e_J^Q(t)|
                  <=C_S epsilon^3 sqrt(1+log(1/epsilon)). (RW2)

The coefficients are precisely CD's k_a,k_b,k_x. RW2 concerns the causal
quadratic responses RW1. The replacement of the whole actual continuum
by these responses is the explicit remaining term in section 5.

## 2. Carrier decomposition with actual envelope control

By definition of the canonical interaction amplitudes,

    C=(a exp(-i lambda t)+conjugate(a)exp(i lambda t))/2,
    R=lambda(a exp(-i lambda t)-conjugate(a)exp(i lambda t))/(2i),
    Z=(b exp(-i omega t)+conjugate(b)exp(i omega t))/2,
    P=omega(b exp(-i omega t)-conjugate(b)exp(i omega t))/(2i).

These identities are algebraic and use canonical momenta. Consequently
Z^2 has carriers 0,+/-2omega with envelopes |b|^2/2,b^2/4,conjugate(b)^2/4.
C^2 has the analogous lambda carriers. ZC has all four carriers
+/-lambda+/-omega, each with the corresponding product envelope divided
by four. Each envelope A satisfies

    ||A||Linfinity<=C epsilon^2, ||A'||L2(0,T)<=C_S epsilon^3. (RW3)

This is MA14, including conjugate products. Initial values are retained.
For example the Z^2 envelopes have nonzero values at time zero.

For a projected profile B and one of its carriers nu, set

    K_B(t)=U_A(t)(0,B),
    P_(B,nu)=integral_0^infinity exp(i nu t)K_B(t)dt,
    E_(B,nu)(t)=integral_t^infinity exp(i nu s)K_B(s)ds.

All vectors are measured in the local weighted phase space X_-6. For the
odd Fa,Fb profiles the LE decay t^-3/2 makes the primitive absolutely
convergent and gives a tail at most C(1+t)^-1/2, including the DC carrier.
For even G the four carriers are uniformly separated from the threshold
by QH4. QH5--QH6 give the same tail bound. Thus

    ||P_(B,nu)||X_-6<=C,
    ||E_(B,nu)(t)||X_-6<=C(1+t)^-1/2.                (RW4)

The open on-shell responses are used only in this local space. No global
energy norm is assigned to a stationary outgoing response.

## 3. The initial transient and envelope-history error are small in work

QH7 gives the exact identity for each branch

    Y_(B,nu,A)(t)=exp(-i nu t)[P_(B,nu)A(t)
       -E_(B,nu)(t)A(0)-integral_0^t E_(B,nu)(r)A'(t-r)dr].

Let Y_ad retain the first term. By RW3--RW4 and Young's convolution bound,

    ||Y-Y_ad||L2(0,T;X_-6)
       <=C epsilon^2 sqrt(log(2+T))+C sqrt(1+T)||A'||L2
       <=C_S epsilon^2 sqrt(1+log(1/epsilon)).        (RW5)

The initial transient is included by the first term of this bound.
Summing the finitely many branches preserves the estimate for both fields.
The modal product weights PZ,PC,RC,RZ have L2 time norm at most C_S epsilon.
Pairing RW5 with the exponentially localized profiles in the work integrals
therefore costs at most

    C_S epsilon^3 sqrt(1+log(1/epsilon))              (RW6)

uniformly in the terminal time. In particular this is o(epsilon^2).

## 4. Averaging the local response retains every channel

The position component of P_(B,nu) is the causal scalar susceptibility
(A-(nu+i0)^2)^-1 B. At positive open frequency its imaginary pairing has
the sign and normalization used in CD. Negative branches are conjugate.
The DC and mixed difference susceptibilities are real.

Insert the algebraic expressions for Z,P,C,R into the work of Y_ad.
Each term is a quartic envelope times exp(i ell t), with ell a combination
of omega and lambda. The zero-frequency terms are exactly

    -(36 d_a/omega^2)I^2-(72 d_x/(omega lambda))I J,
    -(36 d_b/lambda^2)J^2-(72 d_x/(omega lambda))I J.  (RW7)

This is CD's work calculation with the complex envelopes evaluated at the
current time. Constant phases can be absorbed into the two oscillator
angles, so the zero-frequency coefficient depends only on |a|^2,|b|^2.
The closed difference response and the real susceptibility parts have
zero zero-frequency action-work coefficient. Cross terms between the
Fa and Fb responses are retained among the nonzero frequencies.

The possible nonzero frequencies are among

    +/-2omega, +/-2lambda, +/-4omega, +/-4lambda,
    +/-2(omega+lambda), +/-2(omega-lambda).

Their magnitudes have a uniform positive lower bound on the parameter
band. For example 2(lambda-omega)>=2(sqrt(3)-sqrt(7/4))>1/2.
For each quartic envelope C_4,

    ||C_4||Linfinity<=C epsilon^4,
    ||C_4'||L1(0,T)
       <=C epsilon^3 sqrt(T)(||a'||L2+||b'||L2)
       <=C_S epsilon^4.

One integration by parts therefore bounds every nonzero-frequency work
integral by C_S epsilon^4, uniformly in its upper endpoint. Combining
this with RW6 proves RW2. No time derivative of a reference susceptibility
is used, since its mass parameter and carrier are fixed.

## 5. Exact remaining whole-continuum work

Define delta_v=v-v_Q and delta_u=u-u_Q. They have zero initial continuum
phase, since CM leaves the physical continuum unchanged. To display the
complete forcing without a chart approximation, reconstruct the exact
physical fields eta=c g+v and psi=z f+u from the canonical state and put

    N_phi=-6K(eta^2+psi^2)-6eta psi^2-2eta^3,
    N_psi=-12K eta psi-6eta^2 psi-g4 psi^3-g6 psi^5.

Then the remainder equations are exactly

    delta_v,tt+H delta_v=P_o[N_phi+6(Fa Z^2+Fb C^2)],
    delta_u,tt+L delta_u=P_e[N_psi+12G ZC].           (RW8)

The physical modes c,z and the canonical C,Z are distinguished in RW8.
It contains every higher-order, continuum-bearing and reconstruction term.
Define the remaining work

    D_I(t)=-(12/omega)integral_0^t P[Z<Fa,delta_v>+C<G,delta_u>]ds,
    D_J(t)=-(12/lambda)integral_0^t R[C<Fb,delta_v>+Z<G,delta_u>]ds. (RW9)

Combining MA9 and RW2 yields the actual identities

    I(t)-I(0)+integral_0^t(k_a I^2+k_x I J)ds
          =D_I(t)+O_S(epsilon^3 sqrt(1+log(1/epsilon))),
    J(t)-J(0)+integral_0^t(k_b J^2+k_x I J)ds
          =D_J(t)+O_S(epsilon^3 sqrt(1+log(1/epsilon))). (RW10)

The current absolute estimate gives |D_I|+|D_J|<=C_S epsilon^2.
Indeed TP2 bounds the whole continuum in L2 time local norm by C epsilon.
The leading local response has size O(epsilon^2), hence L2 time norm
O(epsilon), and RW5 is smaller. Thus delta has L2 time local norm
O(epsilon). The work weight has L2 norm O(epsilon), giving only O(epsilon^2).
This exact borderline estimate is the unpaid step in admitting damping.

Next estimate RW9 using the full source RW8 and the oscillatory modal
weights, seeking sup_(t<=T)(|D_I|+|D_J|)=o(epsilon^2). A stronger standalone
field approximation is sufficient but may be unnecessary. Control of the
actual work pairing is the target. Critical damping, independent proof
review and empirical validation remain separate obligations.
