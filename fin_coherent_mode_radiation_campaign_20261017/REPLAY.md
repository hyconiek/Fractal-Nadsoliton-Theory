# Bounded replay / reproducibility

Run `python REPLAY/verify.py` (mpmath 1.3+, numpy not required). It checks the exact f norm, PDE mode eigenresidual (symbolic identities via mpmath point tests), noncommuting sample Hessian and certified lower comparisons at family endpoints, the scattering ODE residual at representative x,q, the Fourier on-shell transform against numerical quadrature, and Gamma positivity. Each check is a numerical/symbolic control **only**, not a PDE, limiting-absorption or nonlinear radiation proof. The analytic justifications are in the theorem files.

The replay uses the mathematical witness a=-1/2 and offdiag K_23=0.01 with d_2=-0.005, d_3=0.005, but this is a diagnostic member of a frozen a interval, not tuning a to improve gamma. No unrestricted parameter sweep is run. Recompute all displayed numerical values; do not assume they are experiment data.
