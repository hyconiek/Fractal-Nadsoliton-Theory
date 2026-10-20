# CPL-001: an explicit process using the existing critical coordinates

Date: 2026-10-20. Diagnostic law, not a finite-N FIN limit theorem.

## 1. The original state and its invertible charts

The target state is (z,Y) in R times the ten-dimensional MT fast space.
Its measure is rho(z)dz times N(0,H_fast^(-1)), with
H_fast=E_fast^T(12I-g_c A7)E_fast>0. Whiten the k3 pair and its complement:

    eta=sqrt(12gamma3) Y_k3 in R^2,
    v=H_rest^(1/2)Y_rest in R^8,
    u=1-exp(-||eta||^2/2), gamma3=1-Lambda3/Lambda6.

Under the fast Gaussian, u is uniform in (0,1), the polar angle theta is
uniform on the circle, and v is an independent standard Gaussian. The original
phase R<=T(z) is EXACTLY u<=p(z), since ||eta||^2=gamma3 R.

On phase A set h=u/p(z); on B set h=(u-p(z))/(1-p(z)). The inverse radii are

    r_A(z,h)=sqrt[-2log(1-p(z)h)],
    r_B(z,h)=sqrt[-2log((1-p(z))(1-h))],
    eta=r_P(z,h)(cos theta,sin theta).                 (1)

For 0<h<1 these are smooth coordinate bijections onto the two open phases,
apart from a null axis. The phase boundary and axis have zero stationary
measure. They are not crossed or reached by the stationary process below.
It is enough to specify the process on this full-measure invariant state set;
no global Feller or smooth-interface assertion is made.

Let mu_P=integral q_P. The original measure restricted to phase P becomes

    q_P(z) dz dh dtheta/(2pi) gamma8(dv).              (2)

Indeed r exp(-r^2/2) r_h equals p(z) in A and 1-p(z) in B. Mixing (2) with
its actual masses mu_A,mu_B recovers the full rho times Gaussian law exactly.
P is recoverable from (z,Y); h,theta,v are coordinates, not extra variables.

## 2. Freeze the diagnostic generator

Put a(h)=h(1-h)^2 and c_P(z,h)=1/(1+r_P(z,h)^2). For any fixed kappa>0
define on each chart

    L_P F=b F_z+(g_P/2)F_zz
       +kappa[partial_h(a partial_h F)
                    +c_P F_thetatheta+Delta_v F-v.grad_v F]. (3)

The default comparison is kappa=1. Other kappa values below only expose the
native scaling gate; they are not a fitted physical clock or sampler sweep.
An explicit SDE realization uses independent Brownian drivers:

    dZ=b(Z)ds+sqrt(g_P(Z))dB,
    dH=kappa a'(H)ds+sqrt(2kappa a(H))dW,
    dtheta=sqrt(2kappa c_P(Z,H))dW_theta modulo2pi,
    dv=-kappa v ds+sqrt(2kappa)dW_8.                  (4)

Initial P has probability mu_P; given P, take Z0 with density q_P/mu_P,
H0 uniform, theta0 uniform, v0 standard Gaussian, independent. Pull (4)
back through (1) and the fixed whitening map. This defines a process on
the original eleven coordinates with no independently observable phase bit.

The choices a,c and OU are transparent COMPARISON dynamics, not derived FIN
kinetics. They are chosen to yield a well-defined reversible test object
with bounded local quadratic characteristics, not adopted as natural laws.
