# Numerical scope and reliability

The limiting 27-cell model was approximately discretized by a finite reversible birth-death chain in scalar z (stationary grid weights proportional to exp(-z^4/12)), a Gaussian increment convolution in the accumulated y coordinate and a first-order interpolation for z*dt advection. The two quantizer times are exactly 1/4 and 1/2; the target is at time1. The grid never replaces the microscopic parent.

Three fixed grids were run, with no threshold/record tuning:

| Grid | hz | hy | dt | approximate nine-bin Bayes risk |
|---|---:|---:|---:|---:|
| COARSE | .10 | .10 | .025 | 0.413546810216493 |
| FINE | .075 | .075 | .0125 | 0.415083163187507 |
| REFINED | .06 | .06 | .0125 | 0.415555605584136 |

Truncations are approximately |z|<=4.5 and |y|<=3.8. Lost y-grid probability was of order 10^-6 in these runs, but this *mass statistic is not a bound for spatial/time discretization or for model approximation*. Grid differences are 0.00154 and 0.00047 in risk, which are convergence diagnostics, NOT a certified error bar. No outward interval and no numeric lower advantage has been proven.

The scientifically admitted numerical-free conclusions are: all fixed 27 limiting probabilities exist as path integrals; all 27 finite-N probabilities converge at unspecified o(1); risks converge; the limiting risk is *strictly* below 1/2 by a monotonicity/association argument. The approximate table should not be used to decide a real experiment.

Independent finite rational generator/hypergeometric unit tests use exact Fraction and never claim an asymptotic theorem. `independent_checks.py` has 5626 enumerated finite assertions, not 480 independent physical proofs.
