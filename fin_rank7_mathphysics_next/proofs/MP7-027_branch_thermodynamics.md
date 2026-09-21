# MP7-027 — stationary branch thermodynamics and event separation

Scientific state: **PROVED_ANALYTIC + INTERVAL_ASSISTED GLOBAL EVENT IN THE DECLARED FINITE MODEL**.

At a stationary dual branch `theta(g)`,

```
Phi_g(theta)=||theta||^2/(2g)-log mean exp(X theta),
grad_theta Phi=0.
```

The envelope theorem is elementary here because the gradient term vanishes:

```
d Phi_branch/dg = partial_g Phi = -||theta||^2/(2g^2).
```

Using stationarity `theta=g mu` gives the equivalent identity

```
d Phi_branch/dg = -||mu||^2/2.
```

The accepted R7P-026 equal-energy event isolates a localized C4 root and gain
near `3.7183448981203875`.  R7P-027 proves that localized root is a strict full
7D local minimum.  R7P-028 certifies the localized-minus-uniform branch-energy
slope at the event as

```
[-0.4524137321765795878, -0.4524137306750521069],
```

strictly negative.  Therefore the local equal-energy crossing is transverse.
It is locally unique once the localized stationary branch is continued by the
positive Hessian implicit-function theorem.

At exact rational `g=3.7183449`, R7P-029 separately isolates a full-H7 index-one
saddle.  Its potential exceeds both the uniform state and the named localized
minimum by more than `0.0465` in the supplied local comparison.  This is a
certified local barrier state but not a theorem that it is the globally lowest
mountain pass.

The fold event near `g=3.5156447168...`, the local equal-energy event near
`g=3.7183448981...`, and any global equilibrium transition are three logically
distinct notions.  Until MP7-016/017 pays the stationary/global complement,
this task does **not** call the equal-energy event the first global transition
and does not call its slope jump a physical latent heat.

## Global upgrade after MP7-016/017

MP7-016 proves that the uniform state is the unique global minimum at `g=37/10`,
and MP7-017 exhausts the stationary complement on the certified equal-energy
gain box.  The localized root and the uniform root are the only global minima
at the crossing; the remaining stationary saddle stays above them by more than
`0.04655543`.

Because, for every fixed probability state `p`, the primal energy

`V_g(p)=D(p||u0)-(g/2)(p-u0)^T A7(p-u0)`

is nonincreasing in `g`, the global exhaustion at the event together with the
unique-uniform result below it promotes the event to the **first global
attainment by a nonuniform orbit** in the supplied finite model.  No physical
temperature interpretation is added.

The localized order parameter at the event has

```
||s|| in [3.53697961814657, 3.53697962211361],
||mu||=||s||/g in [0.951224190898289, 0.951224192476810].
```

The uniform branch has `mu=0`.  Thus the global equilibrium order parameter
jumps discontinuously by the displayed amount.  The derivative of the global
minimum value with respect to the supplied gain changes from `0` on the uniform
branch to

```
[-0.4524137321765796, -0.4524137306750521]
```

on the localized branch.  This is a certified **dimensionless slope
 discontinuity in gain**, not a latent heat unless a separate temperature/
energy convention is supplied.

Accordingly the scientific state of MP7-027 is now: **PROVED GLOBAL EVENT IN
THE DECLARED FINITE MODEL**, while the dynamical escape interpretation of the
index-one saddle remains local/conditional.
