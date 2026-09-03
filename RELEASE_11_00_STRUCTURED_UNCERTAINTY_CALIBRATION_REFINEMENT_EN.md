# Release 11.00 — Structured Uncertainty, Mixture E-Process, Calibration and Refinement Robustness

## Record metadata

- **Creator:** Żuchowski, Krzysztof
- **Affiliation:** Independent Researcher — Fractal Information Theory Project
- **ORCID:** 0009-0002-0909-3613
- **Resource type:** Publication — Preprint
- **Version:** 11.00
- **Publication date:** 2026-08-31
- **Language:** English
- **Publisher:** Zenodo
- **Access:** Open access
- **License:** CC BY 4.0

## Abstract

Release 11.00 reports six research rounds, ST1452–ST1541, extending the
interval-certified conditional FIN Operational/State discriminator through
structured matrix uncertainty, optimized time allocation, a finite-mixture
e-process, calibration polytopes and exact 12-to-24 refinement.

A row-sum-preserving real symmetric circulant matrix box with entry radius
\(5\times10^{-15}\) gives eigenvalue radius \(6\times10^{-14}\) while
preserving Fourier eigenvectors. Every admissible spectrum retains one
return-gap maximum in the common bracket

\[
[0.5945307,0.5945309],
\]

and the four-time composite certificate remains

\[
\mathrm{RSS}\ge0.04669201121378897.
\]

This does not cover noncirculant perturbations or symbolic kernel-parameter
uncertainty.

A numerical minimax linear programme over the declared nuisance grid assigns
approximately 31.55 percent of samples to time 0.3 and 68.45 percent to time
2.0. Its grid objective is approximately 1.67 times the equal-allocation
value. A finer held-out grid preserves the result, and an independent interval
cover excludes exact matching at these two times. An adaptive 1,380-box
continuous cover and a two-point dual bound certify

\[
0.05\le V_*\le0.05370588868619479.
\]

The exact optimal weights are not proved, so the four-time design remains the
default for the independent RSS certificate.

A frozen 18-component nuisance mixture defines

\[
E_n=\sum_j\pi_j\prod_{i=1}^n
\frac{q_j(X_i\mid t_i)}{p_C(X_i\mid t_i)}.
\]

Under the classical null, \(E_n\) is a nonnegative mean-one martingale;
crossing 100 therefore has anytime Type-I error at most one percent. The test
is one-sided and finite-grid scoped.

A common calibration polytope with at most five percent uniform preparation
contamination, two percent false-positive/false-negative rates and at least
80 percent click efficiency preserves at least 91.2 percent of the ideal gap
and gives

\[
\mathrm{RSS}\ge0.03883580017500168.
\]

Independent calibration maps require tighter bounds: a 3-percent preparation
plus 2-percent detector example passes, whereas 5 percent plus 2 percent is not
certified by the present sufficient condition.

The exact refinement

\[
A_{24}=A_{12}\otimes I_2+I_{12}\otimes
\begin{pmatrix}q&-q\\-q&q\end{pmatrix}
\]

transports every coarse return and certificate exactly, independently of
\(q\). A fiber-resolving instrument adds heat and unitary fiber returns
\((1+e^{-2qt})/2\) and \(\cos^2(qt)\), but does not derive \(q\).

## Scientific boundary

The robustness ladder is conditional on a circulant uncertainty model, common
calibration map, finite nuisance grid and exact Kronecker refinement. No
physical state category, apparatus, raw event, independent custody or Theory
of Everything is claimed.

## Included files

- English PDF and LaTeX source
- six research scripts
- mixture grid/e-process helper and fixtures
- calibration polytope and protocol 11.00
- six result sets and 90 evidence packets
- summary figure, tests, source bundle and SHA-256 manifest
