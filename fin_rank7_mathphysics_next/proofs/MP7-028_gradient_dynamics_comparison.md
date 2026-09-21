# MP7-028 — two explicitly supplied gradient dynamics

Scientific state: **PROVED_ANALYTIC, CONDITIONAL_ON_ADDED_DYNAMICS**.

Let `Phi(s)` be the static aligned C4 potential and supply the evolution law

```
dot s = -L(s) grad Phi(s).
```

This is an additional dynamical premise; it is not derived from the static
potential.

## Case 1: L=I

Along every classical solution,

```
d Phi/dt = -||grad Phi||^2 <= 0.
```

At a critical point `s*`, linearization is `-H`, where `H=Hess Phi(s*)`.
Thus a positive Hessian gives linear decay and each negative Hessian direction
is linearly unstable.

## Case 2: positive diagonal mobility

Let `L(s)=diag(l_i(s))` with continuous `l_i(s)>0`.  Then

```
d Phi/dt = -grad Phi^T L grad Phi <= 0.
```

On the nonnegative C4 cone, if `s_i=0`, the accepted nonnegative-mean property
gives

```
partial_i Phi = s_i/g-mu_i = -mu_i <=0,
dot s_i = l_i mu_i >=0.
```

Hence the vector field points inward on each coordinate face; the cone is
forward invariant whenever the usual local existence hypotheses for the ODE
hold.

At a critical point derivatives of `L(s)` multiply `grad Phi=0`, so the
linearization is exactly

```
-L_* H.
```

For symmetric positive definite `L_*`,

```
-L_* H
```

is similar to the symmetric matrix

```
-L_*^{1/2} H L_*^{1/2}.
```

By Sylvester inertia under congruence, changing positive mobility changes rates
and eigenvectors but not the count of stable/unstable Hessian directions.

A general SPD matrix with off-diagonal entries does **not** automatically
preserve the nonnegative cone: at a coordinate boundary its ith velocity mixes
other gradient components.  Such a mobility needs a separate boundary proof.

## Scope

The two dynamics have the same equilibria and a common Lyapunov potential but
need not have the same trajectories or relaxation rates.  No unique physical
mobility or physical clock follows from the static FIN landscape.
