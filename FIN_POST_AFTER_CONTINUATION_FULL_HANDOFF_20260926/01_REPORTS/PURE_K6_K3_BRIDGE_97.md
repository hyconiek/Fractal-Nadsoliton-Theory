# PURE-K6–K3-BRIDGE-97
## A subcritical k3 crossing connects the pure-k6 branch to the exact k3+k6 two-harmonic family

Date: 2026-09-26

Status:
- exact scalar crossing condition;
- high-precision crossing root;
- exact cumulant/Schur fourth-order coefficient;
- direct numerical amplitude replay;
- the downstream two-harmonic/D3 crossing is already interval-assisted in MP7-039.

## 1. Pure-k6 parity law

On the nonzero pure-k6 branch define

    q = tanh(J6).

For the two real k3 directions, parity separates their variances:

    Var(k3c)
      = (lambda3/12)(1+q),

    Var(k3s)
      = (lambda3/12)(1-q).

So only one of them softens first.

The k3c crossing condition is

    boxed:
    1/g
      = (lambda3/12)[1+tanh(J6)],

together with

    J6=(g lambda6/12)tanh(J6).

Solving these two scalar equations gives

    boxed:
    g_36 = 5.180490619637444,

    J6 ≈ 0.182996029531570.

This lies below the exact D6/k5 crossing at g=12/lambda5.

## 2. Symmetry change

The pure-k6 parent has D6 stabilizer.

A nonzero k3c component reduces this to the D3 stabilizer of the exact
two-harmonic family

    h_j
      = J3 cos(pi j/2)
        +J6 (-1)^j.

Thus the full D12 orbit size doubles:

    D6 parent:
      orbit size 2

    D3 two-harmonic daughter:
      orbit size 4.

## 3. Subcritical normal form

At the crossing, all noncritical Hessian directions are positive.

Lyapunov-Schmidt elimination gives

    boxed:
    D4 Phi_eff
      ≈ -3.40061736744598 <0

for the normalized k3c critical coordinate.

The crossing slope is

    d lambda_k3/dg
      ≈ -0.2913273628.

Therefore the daughter is subcritical and exists on the lower-g side.

Its leading amplitude is

    theta3c
      ≈ 0.71694754 sqrt(g_36-g).

Direct solves of the exact two-harmonic equations give:

    delta=1e-4:
      theta3c/sqrt(delta)≈0.7170153

    delta=1e-5:
      ≈0.7169543

    delta=1e-6:
      ≈0.7169482,

confirming the normal-form coefficient.

## 4. Direct connection to MP7-039

Following this D3 two-harmonic daughter downward in g reaches the already
interval-certified MP7-039 crossing at

    g_D3≈5.171841831942818.

Between

    g_D3 < g < g_36

the two-harmonic branch has full Morse index 1.

At MP7-039 its two-dimensional D3 critical block reorganizes the branch
according to the nonlinear classification of report 66.

So the earlier MP7-039 parent is no longer an isolated special family:
it is a secondary symmetry-broken branch of the pure-k6 stationary component.
