# 317 — FROZEN-N10-PREDICTOR + DIRECT-N10-HOLDOUT

Date: 2026-09-27

## Status

**PASS** for the pre-frozen canonical two-time observed-process prediction at direct microscopic `N=10`.

This task was executed in the required order:

1. fit/freeze the `N=10` predictor from `N=3..9` only;
2. hash the frozen prediction;
3. only then construct the `N=10` microscopic count-space generator;
4. open the holdout without widening the report-313 process envelope.

The report-313 primary envelope remained fixed at

\[
B_{313}=0.05582267457195117.
\]

The optional report-316 simplex-native correction was scored independently and was not substituted for the canonical model after opening.

---

## 1. Frozen prediction before any N10 microscopic construction

The frozen model used the same architecture as tasks 304–316:

- microscopic slow clock law `log rho_N ~ a+bN`, refit on opened `N=3..9`;
- linear-in-`N` dimensionless slow-shape law for the five nontrivial `R_k=-lambda_k/rho`, with `R_4=1` exact;
- the previously selected initial-slip feature families;
- mandatory Euclidean projection of the predicted 12-state preparation onto the probability simplex;
- optional report-316 rank-1, quadratic-in-`m` residual correction with shrinkage `s=0.20`.

Frozen canonical prediction:

\[
\rho_{10}^{pred}=0.002148880208244716,
\]

\[
R_{10}^{pred}=
(1.44506045,\ 1.68856004,\ 0.91744369,\ 1,\ 1.26833955,\ 1.24989197).
\]

At the fixed dimensionless second-step separation `Delta tau=0.5`,

\[
\Delta t_{10}^{pred}=232.679326693794.
\]

The raw preparation formula again crossed the simplex boundary:

\[
\min_j p_j^{raw}=-0.00846377,
\]

so the report-316 simplex sanitation is mandatory for probability-level interpretation.

Before opening the holdout, the predicted minimum FIN/comparator separation was

\[
\boxed{7.883216\%\ TV}.
\]

Therefore the pre-certified margin against the independently frozen report-313 envelope was

\[
7.883216\%-5.582267\%
=
\boxed{2.300949\ \text{percentage points}}.
\]

Frozen artifact hash (content hash stored inside the artifact):

`0356ccf6c694d0cb3e3953f8b10bfd9d984355af6ddef29abc57bbb2e07ecab7`

---

## 2. Direct microscopic N10 construction

The exact finite-`N` leave-one-out Gibbs heat-bath count-space process has

\[
\boxed{352716\text{ states}}
\]

and the sparse generator contains

\[
\boxed{22523436\text{ nonzero entries}}.
\]

Stationarity check:

\[
\|\pi Q\|_1
=1.80\times10^{-15}.
\]

### Conservative basin classification

The high-gap nearest-localized-minimum rule used the same threshold already validated on `N=7`:

`2.3812285731142`.

On the fully exactified `N=9` dataset, **zero high-gap states are misclassified** by this rule.

At `N=10` the deliberately unlabelled boundary layer contains

`82572` count states,

but stationary mass only

\[
\boxed{0.1148205\%}.
\]

No uncertain state was force-assigned merely to improve the score. Its possible contribution is charged explicitly below.

A full exactification of all 3602 ambiguous orbit representatives was attempted, but is unnecessary for the pass and is computationally much more expensive than the conservative bound. No result from the incomplete exactification attempt is used.

---

## 3. Direct slow spectrum

A reversible symmetric representation of the full `N=10` generator was diagonalized for the slow sector.

Rank 16 gives the eleven nonzero slow modes plus four fast controls. The first fast eigenvalue is already near

\[
-0.45155,
\]

well separated from the slow cluster around `10^-3`.

The six assigned microscopic slow rates are

\[
\begin{aligned}
\lambda_1&=-0.00301946344766,\\
\lambda_2&=-0.00354824398405,\\
\lambda_3&=-0.00189639084289,\\
\lambda_4&=-0.00206014358214,\\
\lambda_5&=-0.00262552792069,\\
\lambda_6&=-0.00258763346283.
\end{aligned}
\]

Thus

\[
\boxed{\rho_{10}^{exact}=0.00206014358214}
\]

and

\[
R_{10}^{exact}=
(1.46565680,\ 1.72232849,\ 0.92051392,\ 1,\ 1.27443929,\ 1.25604520).
\]

The clock prediction error is

\[
\boxed{4.3073\%}.
\]

The largest `R_k` relative error is about

\[
\boxed{1.96\%}
\]

for `k=2`.

### Numerical convergence

Using the same computed eigenspace, truncating from rank 16 to the full slow rank 12 changes at most

\[
TV_{prior}=5.0\times10^{-8},
\]

\[
TV_{joint}=1.10\times10^{-7}.
\]

The maximum rank-16 eigenpair residual is

\[
1.22\times10^{-7}.
\]

These are numerical convergence diagnostics, not interval-arithmetic spectral proofs.

---

## 4. Canonical direct N10 holdout

For the thirteen predeclared pinning preparations `kappa=0,...,12`, the canonical simplex-projected prediction gives

\[
\boxed{
\max_\kappa TV(p_{micro},p_{pred})=1.55371\%
}
\]

on the conservatively labelled core.

The full two-time observed joint law gives

\[
\boxed{
\max_\kappa TV(J_{micro},J_{pred})=2.18454\%
}
\]

with mean

\[
1.97551\%.
\]

This is the direct holdout score before adding classification ambiguity.

---

## 5. Explicit ambiguity and numerical error budget

The unlabelled layer is charged pessimistically rather than silently classified.

Worst preparation ambiguity bound:

\[
\boxed{1.36066\%}.
\]

Worst missing labelled path mass after spectral propagation:

\[
\boxed{0.22951\%}.
\]

Rank-12 versus rank-16 joint difference:

\[
1.10\times10^{-7}.
\]

Adding these terms conservatively to the observed core joint error gives

\[
\boxed{
B_{N10}^{cert}=3.25634\%\ TV.
}
\]

This remains substantially below the envelope frozen before `N=9`:

\[
3.25634\% < 5.58227\%.
\]

Therefore

\[
\boxed{\textbf{317 canonical holdout = PASS}.}
\]

The margin to the frozen process envelope is approximately

\[
\boxed{2.326\text{ percentage points}}.
\]

---

## 6. Error decomposition

On the core the same explicit decomposition used in report 313 gives

\[
\epsilon_{red}^{max}=0.09380\%,
\]

\[
\epsilon_{trans}^{max}=0.93596\%,
\]

\[
\epsilon_{prep}^{max}=1.47066\%.
\]

The summed triangle budget is

\[
\boxed{2.50042\%},
\]

while the actual maximum joint error is

\[
2.18454\%.
\]

Two conclusions follow.

1. The microscopic reduction/history defect is now very small.
2. Preparation transfer remains the largest component, but transition/clock transfer has grown appreciably relative to `N=8,9`.

The latter is consistent with the clock error growing from about `0.44%` at `N=8`, `2.28%` at `N=9`, to `4.31%` at `N=10` under the simple exponential-in-`N` clock law.

---

## 7. Optional task-316 correction

The pre-frozen rank-1 simplex-native correction was scored independently.

It improves the core preparation error to

\[
1.47034\%,
\]

and the core joint error to

\[
\boxed{2.09975\%}.
\]

Its conservative certified upper bound is

\[
\boxed{3.18115\%}.
\]

Thus the correction again helps slightly, but it is **not needed** for the holdout pass and is not promoted over the independent report-313 safety envelope.

---

## 8. Comparator discrimination

The pre-frozen predicted minimum separation was

\[
7.883216\%\ TV.
\]

The direct microscopic core is in fact at least

\[
\boxed{9.24720\%\ TV}
\]

from the specified same-rho/same-exit comparator on the tested grid.

Using only the pre-frozen prediction and the conservative full error certificate gives

\[
7.883216\%-3.256336\%
=
\boxed{4.62688\%\ TV}.
\]

Hence the comparator remains separated even after charging the entire ambiguity layer pessimistically.

---

## 9. What 317 establishes

Across the successive direct holdouts `N=8,9,10`, the observed-process prediction continues to generalize without widening the report-313 envelope.

At `N=10`:

- exact microscopic state space is substantially larger (`352716` states);
- the two-time joint law remains within a `3.26%` conservative certificate;
- the fixed `5.58%` process envelope is not challenged;
- the optional simplex-native preparation correction helps but remains secondary;
- the reduction/history component is already below `0.1%`;
- the simple cross-`N` clock law is becoming the next visible weakness.

This is a controlled finite-`N` effective-theory result. It is **not** an `N->infinity` theorem and does not derive physical space, QM, GR, SI units, or the microscopic update law from a deeper principle.

---

## Verdict

\[
\boxed{\textbf{PASS}}
\]

for the frozen canonical `N=10` observed-process holdout.

The strongest current bottleneck is no longer microscopic history closure. It is the joint problem of

1. preparation transfer, and
2. increasingly curved cross-`N` clock scaling.

The next task should therefore improve the clock law using the already mapped metastable barrier structure **before opening another larger holdout**.
