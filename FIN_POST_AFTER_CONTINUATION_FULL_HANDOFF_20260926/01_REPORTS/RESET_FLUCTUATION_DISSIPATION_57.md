# RESET-FLUCTUATION-DISSIPATION-57
## Empirical maximum-entropy refresh gives exact Gaussian fluctuation-dissipation at equilibrium

Date: 2026-09-26

Status:
- exact system-size covariance calculation for the declared independent-label
  refresh process;
- exact Gaussian fluctuation-dissipation identity around an interior stable
  fixed point;
- physical interpretation remains conditional on the heat-bath process.

## 1. N-label empirical process

Let N labels refresh independently at unit rate.

At empirical state p, an i->j event has total rate

    N p_i q_j(p)

and changes the empirical distribution by

    Delta p=(e_j-e_i)/N.

The deterministic drift is

    sum_ij N p_i q_j Delta p
      = q-p.

## 2. Diffusion matrix

For the sqrt(N)-scaled fluctuation, the jump covariance is

    D(p,q)
      = sum_ij p_i q_j
          (e_j-e_i)(e_j-e_i)^T.

At equilibrium p=q=p*,

    boxed:
    D_* = 2[diag(p*)-p*p*^T]
        = 2 S_*.

So the same Fisher/Onsager matrix controls the Gaussian noise amplitude.

## 3. Linear fluctuation equation

Let

    H_V=diag(1/p*)-gA7

on the simplex tangent space.

From report 56 the drift Jacobian is

    B=-S_* H_V.

The Gaussian fluctuation SDE is therefore

    dxi
      = -S_* H_V xi dt
        + sqrt(2S_*) dW_t

on the tangent subspace.

## 4. Stationary covariance

If H_V is positive definite on the tangent space, define

    C=H_V^{-1}

there.

Then

    B C + C B^T + 2S_*
      = -S_* -S_* +2S_*
      =0.

Hence

    boxed:
    Cov_stationary(xi)=H_V^{-1}.

So the local static Hessian and the dynamic Gaussian fluctuation covariance are
inverse objects under the heat-bath fluctuation-dissipation relation.

## 5. Visible coordinate form

For retained fluctuation coordinates

    x=X^T xi,

the accepted stationary expansion has

    drift=(gF-I)x,
    diffusion generator=F:Hess,

equivalently noise covariance `2F`.

This is exactly the projected version of the label-space identity above.

## 6. Event-rate gauge

At refresh rate rho,

    B -> rho B,
    D -> rho D.

The stationary covariance solving the Lyapunov equation is unchanged.

Therefore equilibrium fluctuation data alone cannot recover rho.

This is another exact form of the physical-clock no-go.

## 7. Significance

Inside the declared heat-bath lane, Fisher geometry is not merely a statistical
analogy:

    S_* is simultaneously
    - categorical covariance,
    - equilibrium Onsager mobility,
    - one-half of the empirical noise covariance.

That three-way identity is exact.

It still does not prove that this stochastic dynamics is fundamental physics.
