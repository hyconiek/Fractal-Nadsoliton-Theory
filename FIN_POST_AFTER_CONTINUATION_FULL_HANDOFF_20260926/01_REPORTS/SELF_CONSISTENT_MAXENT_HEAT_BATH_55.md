# SELF-CONSISTENT-MAXENT-HEAT-BATH-55
## Two nested maximum-entropy principles reproduce the declared FIN heat-bath target and give an exact H-theorem

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- exact finite-dimensional variational theorem;
- exact Lyapunov/H-theorem for the deterministic mean heat-bath ODE;
- conditional on the supplied retained feature map X7 and coupling g;
- no physical temperature or clock is derived.

## 1. Existing FIN primal object

Let

    u0=(1/Q,...,1/Q),
    mu(p)=X^T p.

The accepted primal functional is

    V_g(p)
      = D(p||u0) - (g/2)||mu(p)||^2.

## 2. Maximum-entropy target at fixed current retained mean

Freeze the current retained mean

    mu = X^T p

and define the field

    h = g X mu.

For a candidate target distribution r, consider

    J_mu(r)
      = H(r) + h^T r
      = H(r) + g mu^T X^T r.

Because Shannon entropy is strictly concave on the simplex and the second term
is linear, J_mu has a unique maximizer.

Lagrange multipliers give

    log r_i = h_i - log Z,

therefore

    boxed:
    q_i(mu)
      = exp[(g X mu)_i] /
        sum_j exp[(g X mu)_j].

This is exactly the declared heat-bath target.

So q(mu) is the unique Gibbs/MaxEnt response to the current retained field.

## 3. Nested update principle

Combine report 54 with the result above.

Given current p:

1. compute only the retained statistic
       mu=X^T p;

2. construct the unique maximum-entropy target
       q(mu)=softmax(gXmu);

3. among all one-event kernels stationary at that q, choose the unique maximum
   conditional entropy kernel
       K_ij=q_j(mu).

This reproduces the declared complete-refresh heat-bath event.

The derivation is conditional on X and g, but no extra acceptance function is
needed.

## 4. Mean deterministic equation

At unit event rate the expected empirical distribution obeys

    boxed:
    dot p = q(p)-p.

A general rate rho gives

    dot p = rho[q(p)-p].

Only rho sets the overall clock scale.

## 5. Gradient of V_g

For interior p,

    partial_i V_g
      = log(12p_i)+1-(gXmu)_i.

Since

    log q_i
      = (gXmu)_i-log Z,

we can write, modulo an additive constant on the simplex,

    grad V_g
      equivalent to
    log(p/q).

## 6. Exact H-theorem

Along

    dot p=q-p,

the irrelevant constant drops out because sum_i dot p_i=0.

Hence

    dV_g/dt
      = sum_i (q_i-p_i) log(p_i/q_i)

      = -D(p||q)-D(q||p).

Therefore

    boxed:
    dV_g/dt
      = -J_Jeffreys(p,q(p))
      <= 0.

Equality holds iff

    p=q(p).

Thus the fixed points of the heat-bath mean dynamics are exactly the stationary
self-consistency points, and V_g is a strict Lyapunov function away from them
in the interior.

With event rate rho,

    dV_g/dt
      = -rho J_Jeffreys(p,q(p)).

## 7. Important consequence

The same object V_g now plays two compatible roles:

STATIC:
    its stationary points encode the accepted FIN landscape;

DYNAMIC:
    the nested maximum-entropy refresh decreases V_g monotonically.

So the heat-bath dynamics is variationally aligned with the existing FIN
static functional rather than being an unrelated relaxation rule.

## 8. What remains supplied

This result still does not derive:
- X7 from a deeper physical law;
- the active coupling g;
- the event rate rho;
- a physical temperature;
- seconds or action units.

The maximum-entropy principle chooses the least-informative update conditional
on those objects.

## 9. Next atom

The dissipation identity suggests an exact gradient-flow structure.

`RESET-ONAGER-GRADIENT-FLOW-56` should construct the positive Onsager operator
whose edge logarithmic means turn

    grad V_g ~ log(p/q)

into

    dot p=q-p,

and compare its equilibrium limit with the Fisher covariance
`S=diag(p)-pp^T`.
