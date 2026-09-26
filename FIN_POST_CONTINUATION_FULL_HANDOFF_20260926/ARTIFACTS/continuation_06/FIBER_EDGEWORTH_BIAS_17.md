# FIBER-EDGEWORTH-BIAS-17 — conditional mean versus ME7 at finite N

Status: **ANALYTIC_EDGEWORTH_THEOREM_WITH_STRICT-STATE NUMERICAL CHECKS**.

The exact invariant law from FINITE-N-REFRESH-INVARIANT-16 factorizes, at fixed
visible mean `mu=X^T p`, into a factor depending only on `mu` times the
multinomial combinatorial weight.  Consequently the hidden-fiber conditional
law is controlled entirely by multinomial geometry.

Let `p*` be an interior stationary state, `eps=N^{-1/2}`, and use the accepted
residual coordinates

    p=p*+eps[d_x+Yz],
    d_x=(R+Y K)x,
    K=H F^{-1},
    Sigma_z=G-HF^{-1}H^T.

At leading order, `z|x` is centered Gaussian with covariance `Sigma_z`.

## First conditional Edgeworth correction

Stirling plus the cubic entropy expansion gives, modulo terms depending only on
`x`,

    log rho(z|x)
      = log rho_0(z)
        + eps { (1/6) sum_i d_i^3/(p_i*)^2
                -(1/2) sum_i d_i/p_i* }
        + O(eps^2).                                           (1)

For

    v_i = Sigma_z Y_i^T,              (a four-vector)
    s_i = Y_i Sigma_z Y_i^T,

Gaussian contraction of (1) yields

    E[z|x] = eps zeta_stat(x)+O(eps^2),

    zeta_stat(x)
      = (1/2) sum_i v_i [ (d_{x,i}^2+s_i)/(p_i*)^2 - 1/p_i* ]. (2)

## Exact split into the maximum-entropy curvature plus a constant bias

Let `p_ME(mu)` be the unique local maximum-entropy distribution with the same
retained mean.  Its expansion has

    p_ME(mu*+eps x)
      = p* + eps d_x + eps^2 Y zeta_ME(x)+O(eps^3).

Define `B=Y^T-KX^T`.  The Gaussian covariance identity

    B S = Sigma_z Y^T                                         (3)

implies that the entire quadratic-in-x part of (2) is exactly `zeta_ME(x)`.
Therefore

    zeta_stat(x) = zeta_ME(x) + b*,                           (4)

where the finite-N mean--mode bias is independent of `x`:

    b* = (1/2) sum_i v_i [ s_i/(p_i*)^2 - 1/p_i* ].            (5)

Thus

    E[p|x]
      = p_ME(mu*+eps x) + eps^2 Y b* + O(eps^3).              (6)

## Strict-state values

`conditional_fiber_bias.py` evaluates (2)--(5) on the strict uniform, saddle and
localized states and checks (4) at multiple random visible `x`.

Uniform:

    b* = 0  (residual norm about 1.2e-16).

Saddle:

    b* ~= (-0.02838877859, 0, -0.01540840386, 0),
    ||b*|| ~= 0.03230079966.

Localized:

    b* ~= (-0.04374406038, 0, -0.03814467509, 0),
    ||b*|| ~= 0.05803928890.

For each nonuniform state, `zeta_stat(x)-zeta_ME(x)` is numerically constant over
all tested `x` to roundoff, as predicted by (4).

## Consequence for the reduced generator

The `O(1/N)` local visible diffusion acquires the bias operator

    J_bias = (1/2) sum_a b*_a (T_a:Hess).                      (7)

Hence:

- at uniform equilibrium, `b*=0`: the true stationary conditional local
  generator and ME7 agree through `O(1/N)`; the leading full-vs-ME difference is
  the four-channel memory kernel;
- at a generic nonuniform stationary state, full-vs-ME contains both the local
  bias (7) and the nonlocal memory kernel.

This is the precise boundary of the earlier informal statement that ME7 equals
the conditional stationary closure.

## Boundary

This is a finite-N asymptotic theorem for the declared count model.  It is not a
physical memory law, measured finite-size correction, or source for the
heat-bath activity/clock.
