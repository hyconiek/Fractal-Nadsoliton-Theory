# RESET-ONSAGER-GRADIENT-FLOW-56
## The self-consistent maximum-entropy refresh is an exact nonlinear gradient flow of V_g

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact interior-simplex theorem;
- exact Onsager operator inherited from the reset chain;
- exact Fisher/covariance limit at equilibrium;
- conditional on the declared X7,g heat-bath target.

## 1. Force variable

For current interior p define

    q=q(p)=softmax(g A7 p),
    r_i=p_i/q_i,
    phi_i=log r_i.

The functional derivative of V_g, modulo an additive simplex constant, is

    delta V_g / delta p_i = phi_i.

## 2. Logarithmic mean

For a,b>0 define

    Lambda(a,b) =
      (a-b)/(log a-log b)

with continuous value Lambda(a,a)=a.

Lambda is positive and symmetric.

Define complete-graph edge weights

    w_ij=q_i q_j Lambda(r_i,r_j).

## 3. Onsager Laplacian

Let L_w be the weighted graph Laplacian

    (L_w phi)_i =
      sum_j w_ij (phi_i-phi_j).

Using the logarithmic-mean identity,

    w_ij(phi_i-phi_j)
      = q_i q_j(r_i-r_j)
      = p_i q_j-p_j q_i.

Hence

    -(L_w phi)_i
      = sum_j (p_j q_i-p_i q_j)
      = q_i-p_i.

Therefore

    boxed:
    dot p
      = q-p
      = -L_w grad V_g.

This is an exact nonlinear generalized gradient flow.

## 4. Dissipation

Because L_w is symmetric positive semidefinite and annihilates constants,

    dV_g/dt
      = -(grad V_g)^T L_w (grad V_g)
      <=0.

The quadratic dissipation equals exactly the Jeffreys divergence:

    (grad V_g)^T L_w (grad V_g)
      = D(p||q)+D(q||p).

Thus report 55's H-theorem is the Onsager dissipation identity of the reset
network.

## 5. Equilibrium limit

At a fixed point p*=q(p*),

    r_i=1,
    Lambda(1,1)=1.

Therefore

    w_ij=p_i^* p_j^*.

The resulting Laplacian is

    boxed:
    L_* = diag(p*)-p* p*^T = S_*.

So the equilibrium Onsager mobility is exactly the categorical covariance
matrix / Fisher covariance already used throughout FIN.

This gives a precise bridge:

    Fisher covariance
      = equilibrium mobility of maximum-entropy refresh.

## 6. Linearized dynamics

The tangent Hessian of V_g is

    H_V = diag(1/p*) - g A7.

On the simplex tangent space,

    S_* diag(1/p*) delta p = delta p.

Therefore

    delta dot p
      = -S_* H_V delta p
      = (g S_* A7-I) delta p,

which is exactly the direct linearization of

    q(p)-p.

So static curvature and linear heat-bath relaxation are linked by the Fisher
mobility.

## 7. Clock scale

With event rate rho,

    dot p=-rho L_w grad V_g.

The entire vector field and dissipation rate are multiplied by rho, while the
fixed points and V_g are unchanged.

Thus this construction still cannot determine physical seconds.

## 8. Meaning for FIN

The declared heat-bath lane can now be expressed as one coherent structure:

    V_g
      -> self-consistent Gibbs target q
      -> maximum-entropy reset
      -> logarithmic-mean Onsager operator
      -> Fisher mobility at equilibrium.

No independently invented local metric is required once the reset event has
been chosen.

The remaining unsourced inputs are the retained interaction/g and the global
event clock.
