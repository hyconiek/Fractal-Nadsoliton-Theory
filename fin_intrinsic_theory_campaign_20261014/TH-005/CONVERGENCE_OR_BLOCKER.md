# TH-005 convergence disposition

**Verdict: OPEN_ONE_LEMMA, NOT A FAILED SCALING.**

The rescaled drift, covariance, jump size and finite-N convention correction are all consistent with the proposed stopped diffusion. What is not yet rigorously paid is the singular-perturbation step that removes six stable retained coordinates plus the four occupation-fibre coordinates on an N^(1/3) time horizon.

A sufficient next theorem would establish, under a frozen prepared-initial-law class and stopped local neighbourhood, a perturbed-test-function estimate for every smooth compactly supported phi(Z):

`G_N(phi o Z + correction_N) -> G phi`

uniformly in probability/integrated occupation, where

`G phi = -g_f(a lambda+b z^2/2) phi' + g_f phi''`.

Together with tightness and uniqueness of the stopped one-dimensional martingale problem, this closes TH-005. Numerical plots, the static transverse Hessian alone, or assuming p=s(theta) at every time do not close it.
