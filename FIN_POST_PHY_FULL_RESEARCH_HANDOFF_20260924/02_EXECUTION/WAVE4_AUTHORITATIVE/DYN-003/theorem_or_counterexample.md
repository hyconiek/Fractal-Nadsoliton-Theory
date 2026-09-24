# DYN-003 theorem / counterexample note

## Theorem (equal-well dissipative front has zero speed)

Let `a(z)` be a finite-energy travelling profile with limits `a(-infinity)=alpha_-`, `a(+infinity)=alpha_+` and `a'->0`. For

`gamma alpha_t = K alpha_xx-U'(alpha)`, `gamma>0`,

substitution `z=x-c t` gives

`K a'' + gamma c a' - U'(a)=0`.

Multiplying by `a'` and integrating over the real line yields

`gamma c integral_R (a')^2 dz = U(alpha_+) - U(alpha_-)`.

Hence degenerate minima imply `c=0` for any nonconstant front.

For

`m alpha_tt+gamma alpha_t=K alpha_xx-U'(alpha)`,

the travelling-wave equation gives the same integrated identity because the inertial derivative term vanishes at the asymptotic equilibria.

## Application to the FIN common phase

All-orders carrier locking leaves a `Z_q` translation orbit. Translation-related minima are exactly degenerate. Therefore a dissipative autonomous travelling wall between adjacent FIN phase vacua is excluded in this reduced variational class.

## What evades the theorem

1. `Delta U != 0` (explicit bias / unequal wells): then `c` has the sign of `Delta U`.
2. A nonconservative drive: energy is injected and the gradient identity changes.
3. `gamma=0` with separately supplied inertia: moving conservative kinks may exist, but no preferred nonzero speed is generated.
4. Breakdown of the controlled reduction: hard/radial modes or memory must then remain explicit rather than being fitted away.
