# Explicit strong-norm and bootstrap payments

Date: 2026-10-20. Supplement to the existing proof audit; no model/ray repair.
Use ||y||_H2^2=||y||_2^2+||y'||_2^2+||y''||_2^2 and the corresponding H1.

## 1. A uniform linear constant and actual point traces

For A=L_phi on odd functions or A=L_a on even functions, A>=gamma=5/4
and its real multiplication potential satisfies |V|<=6. If g=||Ay||_2,

    ||y||_2<=4g/5,
    ||y'||_2^2<=116g^2/25,
    ||y''||_2<=29g/5,
    ||y||_H2^2<=973g^2/25<49g^2.

For velocity v, ||v||_2^2<=4||A^(1/2)v||_2^2/5 and
||v||_H1^2<=33||A^(1/2)v||_2^2/5<9||A^(1/2)v||_2^2.
The graph energy group and ||A^(1/2)F||_2<=sqrt(6)||F||_H1 therefore give
for ZERO data the safe explicit bound

    ||(y,y_t)(t)||_(H2 x H1)<=21 integral_0^t ||F(s)||_H1 ds. (1)

For smooth compact f, |f(R)|^2=2integral_(-infinity)^R f f' is at most
||f||_H1^2; density extends the uniform trace bound to H1. This applies
to y_t and y_x, not to an unsupported derivative trace of energy data.
The same inequality yields ||fg||_H1<=3||f||_H1||g||_H1. Multiplication
by K=tanh x has H1 operator norm at most 2 since |K|,|K'|<=1.

## 2. Fixed reference constants, without bootstrap dependence

With f=sqrt(3/4)sech^2 x, direct y=tanh x integration gives

    ||f||_2^2=1, ||f'||_2^2=4/5, ||f''||_2^2=16/7.

Since omega^2<=7/4, the prepared v=f cos(omega t) has
||(v,v_t)||^2<=1013/140<9, so choose C_v=3. The polynomial integral for
S=(9/4)y(1-y^2)^2 gives ||S||_H1^2=270/77<100. Equation (1) then allows
||(u,u_t)||<=420(1+t), so choose C_u=420. No numerical quadrature is needed.

## 3. Explicit order of constants for signed finite g4,g6

Write tau=1+t and X=||(r,r_t)||, Y=||(w,w_t)||. Bootstrap
X<=A0 epsilon^3 tau^2 and Y<=B0 epsilon^4 tau^3. Provided
A0 epsilon^2 tau^2<=1 and B0 epsilon^2 tau^2<=1,

    ||psi||_H2<=4epsilon, ||eta||_H2<=421epsilon^2 tau.

The exact source terms, product constant 3 and K multiplier constant 2
give, with epsilon^2 tau<=1,

    ||F_r||_H1<=C1 epsilon^3 tau,
    C1=38405304+576M4+82944M6,
    ||F_w||_H1<=C2[epsilon X+X^2+epsilon^4 tau^2],
    C2=1400000000.                                   (2)

The contributions to the constant part of C1 are 121248 from 12K eta psi
and 38284056 from 6eta^2 psi. For C2 the three pure-source coefficients
sum to 1349876718; the other coefficients 216 and 36 are smaller than C2.
All absolute signs of g4,g6 have been retained via M4,M6.

Choose, in this order,

    Cr=1+21C1, Cw=1+21C2(Cr+2),
    A0=4Cr, B0=4Cw,
    K0=1/[16(1+A0+B0+Cr^2)].                         (3)

Under epsilon^2 tau_T^3<=K0, the first part of (1),(2) improves X to
Cr epsilon^3 tau^2, independently of A0,B0. Then Cr^2 epsilon^2 tau^2<=1
pays the X^2 term, and (1),(2) improve Y to Cw epsilon^4 tau^3.
Both improve their bootstrap by a factor at least four. All smallness
requirements follow from the SAME K0, since tau>=1. Local polynomial
H2 x H1 wave continuation therefore closes the actual zero-data remainder
problem through T. This is not a circular simultaneous choice of constants.

These very large constants are deliberately sufficient, not numerical
performance advice or a claim of a useful epsilon_*.
