# SQ-002: stationary reflected critical process theorem

Date: 2026-10-20. Conditional on actual MT geometry, Gibbs moments and
spectra, with exact SQRT characteristics in
[MICROSCOPIC_CHARACTERISTICS.md](MICROSCOPIC_CHARACTERISTICS.md).
Not mechanized or independently peer reviewed.

## 1. Pay the stationary drift and bracket over the whole horizon

Write a=chi lambda, f(z)=q(a z-z^3/3), and let beta_N^D,gamma_N^D be
the reflected critical drift/bracket. Put e_N=beta_N^D-f(Z_N).
MT supplies, for every fixed power r,

    sup_N E_piD(|Z_N|^r+||sqrt(N)w||^r)<infinity.

It also supplies Z_N(0)=>rho_lambda, with density proportional to
exp(a z^2/2-z^4/12). Let epsilon_N^partial=pi_D(B_D). The actual boundary
weight is an UPPER bound: at most polynomial times exp(-Nb). This uses
the corrected MT direction and its positive local partition lower bound.
In particular N^r epsilon_N^partial->0 for every fixed r.

Equations SQ-001(1),(4),(6) give pointwise on the fixed MT domain

    |e_N|<=C[N^(-1/4)+N^(-1/2)(|Z_N|^3+|Z_N|^5)
                       +N^(3/4)||w||^2]
                  +2R N^(3/4)1_B_D.                (1)

The first term is the true finite-N LOO error. Therefore for p=1,2,4,

    E_piD|e_N|^p
      <=C_p[N^(-p/4)+N^(3p/4)epsilon_N^partial]->0,
    sup_N E_piD|beta_N^D|^p<infinity.               (2)

The O(N^(-p/2)) line remainder was bounded by the displayed slower power.
All ten fast coordinates are paid by their actual high Gibbs moments;
no fast-process theorem, NC/KR drift matching or constant full gap is used.
Stationarity converts (2) into bounds over [0,S]. In particular

    E sup_(t<=S)|integral_0^t e_N ds|^2
       <=S^2 E_piD|e_N|^2 ->0.                      (3)

Likewise SQ-001(1),(7) and the actual boundary correction give

    E sup_(t<=S)|integral_0^t gamma_N^D ds-2q t|
      <=S E_piD|gamma_N^D-2q|
      <=C_S[N^(-1/4)+epsilon_N^partial]->0.          (4)

This is bracket convergence, not an inference from its stationary mean.

## 2. Maximum bounds and tightness, not just one-time moments

For the full count filtration the exact compensated martingale satisfies

    Z_N(t)=Z_N(0)+integral_0^t beta_N^D ds+M_N(t),
    <M_N>_t=integral_0^t gamma_N^D ds<=4Rt.          (5)

A martingale jump is a Z jump, bounded by j_N=2N^(-3/4). A fourth-power
jump expansion, first stopped at |M_N|<=L, has compensating integrand

    sum_j rate_j[6 M_-^2 j^2+4 M_- j^3+j^4].

Since sum rate_j j^2<=4R, |j|<=j_N, and the stopped second moment is at
most 4Rt, its absolute upper expectation integrates to at most

    4R integral_0^S [6(4Rt)+4j_N sqrt(4Rt)+j_N^2]dt.

Fatou removes the stop; Doob's L4 inequality then gives
sup_N E sup_(t<=S)|M_N(t)|^4<infinity. Finite-jump localization justifies
the expansion and zero mean of its martingale term. This proves UI of
M_N(t) and M_N(t)^2, not merely bounded expected brackets.

Using (2),(5) and Holder,

    E sup_(t<=S)|Z_N(t)|^4
      <=C[E|Z_N(0)|^4+S^3 E integral_0^S|beta_N^D|^4 ds
                              +E sup|M_N|^4]<=C_S.   (6)

This is genuine maximum containment. For any stopping time tau<=S and
deterministic 0<=theta<=1, a pathwise Holder bound and stationarity yield

    E|integral_tau^(tau+theta) beta_N^D ds|^2
       <=theta E integral_0^(S+1)|beta_N^D|^2 ds<=C_S theta,
    E|M_N(tau+theta)-M_N(tau)|^2<=4R theta.          (7)

In particular one never replaces the law at a stopping time by pi_D.
The stopping-time tightness criterion, (6),(7) and vanishing jumps show
C-tightness in D([0,S]),J1. The same estimates give joint tightness with M_N.

## 3. Joint martingale identification and uniqueness

Along a convergent subsequence limits are continuous. On bounded path
sets the cubic drift integral is a continuous functional; containment (6)
localizes, and (3) removes the residual. Hence

    M(t)=Z(t)-Z(0)-integral_0^t f(Z(s))ds.

For 0<=s<t test the finite-N martingale and square-minus-bracket identities
against any bounded continuous cylinder of (Z_N,M_N,Z_N(0)) at times <=s.
Joint convergence at continuous paths and the uniform fourth bounds pass
these expectations. Equation (4) identifies M^2-2q t as a martingale.
A monotone-class extension gives both identities in the joint natural
past filtration, including Z(0). No limiting fast-path filtration is claimed.

Thus M is continuous with <M>_t=2q t. The continuous-martingale Brownian
characterization gives M=sqrt(2q)B in this same filtration. In particular
future Brownian increments are independent of the initial quartic variable.

Local Lipschitz drift gives a unique local solution. For its generator
G=q partial_z^2+f partial_z,

    G(z^2)=2q a z^2-(2q/3)z^4+2q
             <=2q+(3q/2)(max(a,0))^2.

Stopped second-moment estimates exclude finite-time explosion. Local
pathwise uniqueness then gives uniqueness in law of the global SDE, so
every subsequence has the same limit. The final theorem is

    Z_N => Z in D([0,S]), J1,
    dZ=q(chi lambda Z-Z^3/3)ds+sqrt(2q)dB,
    Z(0)~rho_lambda, B future-independent of Z(0).     (8)

The flux identity f rho_lambda-q rho_lambda'=0, with quartic tails, checks
the invariant law; stationarity also passes from the finite stationary
processes. At lambda=0 its fourth moment is 3, consistent with (8).
The limit has the law of native stationary critical HB observed at q times
its model time. This is a derived critical-amplitude universality statement,
NOT equality of finite-N laws or physical clock calibration.

No rate of weak convergence, labelled history, posterior/information limit,
global equilibrium, source of the coefficients or physical Brownian field
is supplied by this theorem.
