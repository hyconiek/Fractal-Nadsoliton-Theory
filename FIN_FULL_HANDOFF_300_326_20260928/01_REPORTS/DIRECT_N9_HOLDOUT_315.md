# DIRECT-N9-HOLDOUT-315
## The frozen joint-law budget passes a new N=9 microscopic holdout; preparation remains the dominant error source

Date: 2026-09-27

Status: direct finite-N count-space generator for N=9 with an explicit conservative classification-ambiguity budget. Phase A was frozen before any N=9 microscopic construction.

## 1. Pre-open contract

Before opening N=9:

- report 313 envelope:

\[
\boxed{B_{313}=0.05582267457\;TV};
\]

- canonical N=9 prediction:

\[
\rho_9^{\rm pred}=0.003904627743,
\]

\[
\Delta t=0.5/\rho_9^{\rm pred}=128.0531802;
\]

- minimal predicted FIN-vs-comparator separation:

\[
\boxed{0.07297249584\;TV};
\]

hence a preclassified margin

\[
\boxed{0.01714982127\;TV}.
\]

Phase-A hash:

`a8ece2f30825a8d69ef0b80271be9a0091eb696a195fcec3462ec324fe3cc718`.

The task-314 correction was also frozen as a separate diagnostic. It was not allowed to replace the canonical model after opening.

## 2. Direct N=9 microscopic construction

The exact leave-one-out count-space generator has

- 167,960 states;
- 10,144,784 nonzero generator entries;
- stationary residual L1 = 1.55e-15.

The N7-validated nearest-minimum gap rule leaves 38,096 states in an explicit ambiguity band instead of assigning them a guessed localized label.

Their stationary mass is only

\[
0.00192659
\]

or 0.1927%.

A separate N=9 validation checked all 65 D12-orbit representatives immediately above the threshold (gap <=2.5) plus 250 deterministic farther-core representatives by full potential descent. There were

\[
\boxed{0\text{ mismatches out of }315}
\]

checks. This is strong validation, not a proof for every core state.

## 3. Spectral replay

The reversible generator gives the slow spectrum

\[
-0.003532701,
-0.003817674,
-0.004699161,
-0.004760723,
-0.005505713,
-0.006387389,\ldots
\]

with a large jump to approximately -0.4446.

The exact coarse slow clock is

\[
\boxed{\rho_9^{\rm exact}=0.003817674052}.
\]

Thus the pre-open clock prediction has relative error

\[
\boxed{2.2777\%}.
\]

Rank-16 versus rank-24 propagation differs by at most

- 2.46e-9 TV for the late prior;
- 3.89e-9 TV for the two-time joint law.

## 4. Canonical prediction error

On the confidently labelled core, the canonical prediction has

\[
\max_\kappa TV(p_{\rm micro},p_{\rm pred})
=
\boxed{2.2942\%}
\]

for the 12-state prior, and

\[
\max_\kappa TV(J_{\rm micro},J_{\rm pred})
=
\boxed{2.5489\%}
\]

for the two-time observed joint law.

The mean joint error is 2.0709% TV.

## 5. Conservative ambiguity certificate

Because no guessed labels are assigned inside the N=9 ambiguity band, two explicit worst-case terms are charged:

1. preparation ambiguity: every ambiguous state is pessimistically allowed to belong to the exact J=0 preparation basin;
2. path-label ambiguity: every trajectory touching the ambiguity band at either observed time is treated as potentially adversarial.

Across kappa=0,...,12:

\[
\max B_{\rm prep-amb}=2.2639\%,
\]

\[
\max B_{\rm path-amb}=0.38495\%.
\]

Combining the actual core joint error, both ambiguity terms and the rank-convergence diagnostic gives the conservative upper bound

\[
\boxed{B_{N9}^{\rm certified}=4.0365\%\;TV}.
\]

This remains below the predeclared report-313 envelope:

\[
4.0365\% < 5.5823\%.
\]

Therefore the direct N=9 holdout is a **PASS** against the frozen envelope.

## 6. Error decomposition

On the core, the maximum triangle components are

\[
\epsilon_{\rm reduction}=0.1786\%,
\]

\[
\epsilon_{\rm transition}=0.4744\%,
\]

\[
\epsilon_{\rm preparation}=2.1691\%.
\]

Their maximum correlated sum is

\[
2.8221\%\;TV,
\]

while the actual core joint error is 2.5489% TV.

The preparation term is again the dominant contribution. The dynamical reduction/history error continues to decrease.

## 7. Comparator separation

The frozen predicted minimum separation was

\[
7.29725\%\;TV.
\]

The direct core microscopic FIN is actually at least

\[
8.38848\%\;TV
\]

from the comparator on the tested kappa grid.

Using only the frozen predicted separation and the full conservative N=9 error certificate gives

\[
\boxed{
TV(P_{\rm true\ FIN},P_{\rm comparator})
\ge 3.2608\%.
}
\]

Thus the pre-certified comparator test remains positive after opening N=9.

## 8. Task-314 correction diagnostic

The rank-2 residual correction was frozen before N=9, but its additive Fourier form already produced a negative raw probability (-0.00525) before the holdout was opened. It therefore cannot be promoted to a canonical probability law.

After the predeclared clipping/renormalization diagnostic, it happens to reduce the core joint max error to approximately

\[
1.9785\%\;TV.
\]

This is scientifically interesting evidence that the low-dimensional residual directions are useful, but it does not rescue the additive correction as a fundamental or canonical map. A future correction must be simplex-native.

## Verdict

**315 PASS** for the frozen report-313 observed joint-law envelope on a genuinely new direct N=9 microscopic holdout.

The strongest current pattern is now repeated through N=7,8,9:

\[
\text{preparation transfer is the dominant residual error},
\]

while

\[
\text{history/reduction and transition-law errors remain smaller}.
\]

No claim of an N->infinity theorem is made.
