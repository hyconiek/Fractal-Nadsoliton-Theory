# FOLD-STOCHASTIC-NORMAL-FORM-64
## Universal N^{-2/3}, N^{-1/3}, N^{1/3} scaling from the exact FIN heat-bath fold

Date: 2026-09-26

Status:
- leading system-size stochastic normal form;
- all leading coefficients are fixed by the accepted fold and heat-bath data;
- valid in the local fold window, conditional on the exact finite-N heat bath.

## 1. Soft deterministic coordinate

Let y be the normalized certified kernel coordinate.

The deterministic center drift is

    dot y
      =
      -g_f[
        a delta+(b/2)y^2
      ]
      + higher order,

where

    delta=g-g_f,
    a<0,
    b>0.

## 2. Finite-N noise in retained dual coordinates

For the empirical label process at equilibrium, report 57 gives label-space
noise covariance

    2S/N.

The retained dual coordinate is

    s=g X^T p.

Therefore its leading noise covariance is

    (2g^2/N) F,

where

    F=X^T S X.

Projecting on the normalized fold vector v gives

    sigma_y^2/N
      =
      (2g_f^2/N) v^T F v.

At the fold the Hessian condition is

    (I/g_f-F)v=0.

Because ||v||=1,

    v^T F v=1/g_f.

Hence exactly

    boxed:
    sigma_y^2 = 2g_f.

Numerically,

    2.651657861
      < sqrt(2g_f) <
    2.651657869.

No noise amplitude is fitted.

## 3. Local stochastic equation

The leading soft-mode SDE is therefore

    dy
      =
      -g_f[
        a delta+(b/2)y^2
      ]dt
      +sqrt(2g_f/N)dW_t

plus higher-order transverse and finite-N corrections.

## 4. Critical rescaling

Set

    delta=N^(-2/3) Delta,
    y=N^(-1/3) Y,
    t=N^(1/3) tau.

Brownian scaling gives

    dW_t=N^(1/6)dW_tau.

All leading terms then become O(1), giving

    boxed:
    dY
      =
      -g_f[
        a Delta+(b/2)Y^2
      ]d tau
      +sqrt(2g_f)dW_tau.

Thus the fold has the finite-size scaling triplet

    parameter window:   N^(-2/3),
    soft amplitude:     N^(-1/3),
    relaxation time:    N^(1/3).

## 5. Scaled local potential

The scaled drift is gradient with local potential

    U_Delta(Y)
      =
      a Delta Y+(b/6)Y^3.

The scaled SDE has mobility g_f and noise sqrt(2g_f), so its local formal
density is proportional to

    exp[-U_Delta(Y)]

inside a confining fold tube.

The scaled saddle barrier is

    C_fold Delta^(3/2).

This exactly matches

    N Delta V_g
      -> C_fold Delta^(3/2)

under delta=N^(-2/3)Delta.

## 6. Relaxation scale

Report 59 gives

    lambda_slow
      ~ c_lambda sqrt(delta),

with

    c_lambda≈0.7907048.

At the critical window,

    lambda_slow
      ~ c_lambda N^(-1/3) sqrt(Delta),

so

    tau_slow
      ~ N^(1/3)/
         [c_lambda sqrt(Delta)].

The N^(1/3) critical time scale is therefore consistent with both the
deterministic eigenvalue and the stochastic rescaling.

## 7. Fluctuation width

Since

    y=N^(-1/3)Y,

the soft standard deviation in the critical window is O(N^-1/3), and

    Var(y)=O(N^-2/3).

This is the same scale as the stable-saddle separation.

That equality of scales is the finite-size meaning of the crossover.

## 8. Interpretation boundary

These are dimensionless finite-size exponents of the conditional FIN heat-bath
model.

They do not determine seconds, a physical temperature, or a laboratory critical
phenomenon.
