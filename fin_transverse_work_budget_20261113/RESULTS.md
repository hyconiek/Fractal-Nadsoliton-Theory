# Work-budget results

All values below are normalized by each fixed seed initial reference energy.

| t | Seed | Energy gain minus one | Phi work | Psi work |
|---|---|---|---|---|
| 100 | position | 4.575869596 | 4.480448260 | 0.095421406 |
| 100 | velocity | 6.522187129 | 6.481105647 | 0.041081558 |
| 116.5 | position | 4.657179960 | 4.561657557 | 0.095522487 |
| 116.5 | velocity | 6.650158737 | 6.607509034 | 0.042649796 |
| 150 | position | 4.548858880 | 4.452249865 | 0.096609126 |
| 150 | velocity | 6.387251847 | 6.348367365 | 0.038884607 |
| 200 | position | 4.331611443 | 4.235993662 | 0.095617933 |
| 200 | velocity | 6.112851622 | 6.076250747 | 0.036601048 |

Maximum normalized balance residual, fine run: 1.718605010109231e-07. Coarse residual: 5.1838828260195036e-06. Maximum cumulative-work change under step halving: 9.046081625996294e-07. Static work is zero and historical energy traces are unchanged.

The Phi contribution dominates the accumulated gain. Between t=116.5 and 200 its velocity-seed work decreases by about 0.53126, while Psi work decreases by about 0.00605. Thus the observed energy decline is predominantly the negative increment of the Phi work in this specified decomposition. Both cumulative works remain positive at T=200; cumulative sign and interval increment must not be confused. This attribution is a descriptive identity for the coupled trajectory, not a counterfactual experiment with one field removed.

Next bounded step: derive and test the exact uniform-plus-alternating carrier reduction (even and odd sites) of the conditional A11 model, then use it for a finite-amplitude perturbation with back reaction. Pointwise polynomials preserve the even/odd carrier subspace and circulant mixing preserves modes 0 and 6. Its coefficient normalization and full-equation residual must be checked before simulation. This would test whether the linear amplification survives energy-conserving nonlinear feedback, rather than continuing indefinite linear horizon scans. No global stability or strict-source result is claimed.
