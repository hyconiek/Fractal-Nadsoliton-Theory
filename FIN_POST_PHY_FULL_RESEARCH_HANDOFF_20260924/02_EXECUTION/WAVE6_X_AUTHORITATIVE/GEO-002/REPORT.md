# GEO-002 — homogenization of a declared relational population

Status: **COUNTEREXAMPLE_WITH_COMPONENTWISE_WEAK_CONTINUUM**.

A new exact stationary population counterexample is built from the already-audited rectangular RNG grids. Choose a global latent aspect variable `A in {1,2}` with probability 1/2 and an independent uniform torus translation U. Conditional on A=a, place the periodic rectangular grid `(Z_{a m} x Z_m)+U`. RNG is exactly the four axis-neighbour graph. The random shift makes the point-process law translation-stationary and each one-point empirical limit is Haar.

For the unweighted RNG Laplacian with density scaling `N=a m^2`, the two first-axis eigenvalues satisfy

`N lambda_x -> 4 pi^2/a`,  `N lambda_y -> 4 pi^2 a`.

Thus every quenched component has a clean weak continuum operator with tensor `K_a=diag(1/a,a)`, but the global stationary mixture is non-ergodic: the event `{A=1}` is translation invariant. The quenched effective tensor remains random, while the annealed quadratic form has tensor

`E[K_A]=diag(3/4,3/2)`,

which is anisotropic. Haar density and stationarity therefore do not select isotropy or a unique homogenized response. The missing population datum is at least an ergodic/mixing law strong enough to eliminate such persistent global latent anisotropy (plus the usual ellipticity/moment control).

Held-out `a=3` gives the predicted limiting ratio 9 without refit. This is a structural counterexample, so larger iid simulations are not an appropriate repair.

This does not show that FIN phase populations violate mixing; no FIN-derived population law has been supplied.
