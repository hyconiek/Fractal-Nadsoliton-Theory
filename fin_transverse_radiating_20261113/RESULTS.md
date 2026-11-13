# Numerical results

| Run | Maximum energy gain | Final energy gain |
|---|---|---|
| static_dx0.1_dt0.02 | 1.000000000000 | 0.999999999813 |
| driven_dx0.1_dt0.02 | 3.270123345092 | 3.214603382399 |
| driven_dx0.1_dt0.01 | 3.270124210969 | 3.214604784256 |
| driven_dx0.05_dt0.01 | 3.270124210829 | 3.214604784124 |

The refined maximum corresponds to a norm amplification of 1.808348476049057. The final gain is optimized over initial phase separately from the maximum. Static maximum energy drift is 1.873432520227425e-10.

All trace integrity, positive Gram spectrum, initial normalization, static control and refinement checks passed. This is finite-time amplification in a two-dimensional seed space, not an exponential-growth theorem.
