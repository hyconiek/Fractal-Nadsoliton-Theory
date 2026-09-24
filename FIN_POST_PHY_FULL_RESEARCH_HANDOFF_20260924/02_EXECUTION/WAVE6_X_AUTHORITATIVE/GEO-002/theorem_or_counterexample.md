# GEO-002 counterexample

**Counterexample.** Let A be 1 or 2 with equal probability and U uniform on the flat two-torus. Conditional on A=a, use a translated rectangular periodic grid with `n_x=a m`, `n_y=m` and its unweighted RNG. The law is stationary and its empirical measure converges to Haar, but it is not ergodic because A is a translation-invariant latent variable. Under density scaling `N=a m^2`, the coordinate symbols converge to `4 pi^2/a` and `4 pi^2 a`. Hence the quenched homogenized tensor depends on A and the annealed tensor is anisotropic.

**Consequence.** Weak continuum behavior can survive even when pointwise isotropy fails, but stationarity/Haar density alone do not select a unique isotropic operator. A mixing/ergodicity or equivalent population-source premise is indispensable.
