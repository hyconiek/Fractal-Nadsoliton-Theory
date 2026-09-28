# 318 — BARRIER-AWARE-CLOCK-LAW

Date: 2026-09-27

## Status

**PASS as a finite-N predictive clock-law improvement; conditional on the mapped B4 barrier lane.**

Task 317 showed that the original cross-`N` law

\[
\log \rho_N=a+bN
\]

still predicts the full observed process well, but its direct clock error grew to `4.31%` at `N=10`. Task 318 therefore tests whether the curvature is consistent with the already mapped metastable barrier

\[
\boxed{B_4=0.662219137127460}.
\]

No `N=11` microscopic data are used anywhere in this report.

---

## 1. Exact opened clock sequence

The direct microscopic slow clocks currently available are

\[
\begin{array}{c|c}
N & \rho_N\\ \hline
3&0.131439786195646\\
4&0.071692277481914\\
5&0.0401867559160965\\
6&0.0226112751871144\\
7&0.0126403864326779\\
8&0.00698967186661062\\
9&0.00381767405205540\\
10&0.00206014358213756
\end{array}
\]

Define the one-step exponent

\[
\beta_N=-\log(\rho_N/\rho_{N-1}).
\]

For `N=4,...,10` this rises approximately

\[
0.6062,\ 0.5788,\ 0.5751,\ 0.5816,\ 0.5925,\ 0.6048,\ 0.6169.
\]

The large-`N` end is therefore curving upward toward the mapped `B4`, rather than remaining exactly linear in `N` with a fixed empirical exponent near `0.59`.

---

## 2. Candidate laws

The following low-complexity candidates were compared on rolling one-step holdouts `N=7,8,9,10`.

### A. Previous baseline

\[
\log\rho_N=a+bN.
\]

### B. Fixed-B4 polynomial prefactor

\[
\log\rho_N=a+p\log N-B_4N.
\]

### C. One-term approach of the one-step exponent

\[
\beta_N=B_4-\frac{c}{N}.
\]

### D. Two-term approach

\[
\boxed{
\beta_N=B_4-\frac{c}{N}-\frac{d}{N^2}.
}
\]

The last form is equivalent asymptotically to an exponential barrier factor times a polynomial-type prefactor, up to a convergent finite-size correction.

The candidate comparison is a model-selection exercise on already opened `N<=10`; it is **not** a blind validation. The selected law is frozen only for future `N=11` use.

---

## 3. Rolling one-step clock errors

Absolute relative errors:

| holdout | linear | fixed-B4 logN | B4-c/N | B4-c/N-d/N² |
|---:|---:|---:|---:|---:|
| 7 | 1.089% | 4.654% | 3.061% | **0.747%** |
| 8 | 0.437% | 4.285% | 2.244% | **1.184%** |
| 9 | 2.278% | 3.499% | 1.360% | **1.626%** |
| 10 | 4.307% | 2.383% | 0.513% | **2.026%** |

Maximum error:

\[
\begin{array}{c|c}
\text{law} & \max |\Delta\rho|/\rho\\\hline
\text{linear}&4.3073\%\\
\text{fixed-B4 logN}&4.6537\%\\
B_4-c/N&3.0614\%\\
\boxed{B_4-c/N-d/N^2}&\boxed{2.0258\%}
\end{array}
\]

Mean absolute error also falls from

\[
2.0278\%\quad\text{(linear)}
\]

to

\[
\boxed{1.3955\%}.
\]

Thus the improvement is not produced by one exceptional point.

---

## 4. Does the improved clock improve dynamics?

A clock law matters only if it improves the transition law, not merely the scalar `rho`.

For each rolling holdout:

1. the resolved `R_k` law was kept in the same linear cross-`N` architecture used by the previous campaign;
2. an exact same-`N` 12-state generator was reconstructed from the direct microscopic slow spectrum;
3. each candidate clock set the physical interval

\[
\Delta t=0.5/\rho_N^{pred};
\]

4. the maximum row TV between predicted and exact 12-state transitions was measured.

Maximum transition-row errors over `N=7..10`:

\[
\begin{array}{c|c}
\text{law}&\max_i TV(T_i^{pred},T_i^{exact})\\\hline
\text{linear}&1.0060\%\\
\text{fixed-B4 logN}&1.5957\%\\
B_4-c/N&1.1243\%\\
\boxed{B_4-c/N-d/N^2}&\boxed{0.4820\%}
\end{array}
\]

The mean transition error falls from

\[
0.6040\%
\]

to

\[
\boxed{0.3705\%}.
\]

So the clock improvement survives at the level of the actual transition semigroup.

---

## 5. Parameter evolution

When the two-term B4 law is fit successively before each holdout, its parameters evolve as follows:

| prediction | c | d | predicted beta |
|---:|---:|---:|---:|
| N7 | 1.1373 | -3.6437 | 0.5741 |
| N8 | 1.0732 | -3.3680 | 0.5807 |
| N9 | 1.0007 | -3.0485 | 0.5887 |
| N10 | 0.9248 | -2.7075 | 0.5968 |

They have not stabilized asymptotically. Therefore the law should be interpreted as a controlled finite-size bridge toward `B4`, not yet as a proven Eyring-Kramers prefactor theorem.

---

## 6. Frozen N11 clock prediction

Fit on all opened `N=3..10` gives

\[
\boxed{c=0.8493265661},
\]

\[
\boxed{d=-2.3626406173}.
\]

Therefore

\[
\beta_{11}^{pred}
=
B_4-\frac{c}{11}-\frac{d}{11^2}
=
\boxed{0.6045335866},
\]

and with `rho_10` used as the last opened anchor,

\[
\boxed{ho_{11}^{pred}=0.00112551655898}.
\]

For comparison, the old linear law trained through `N=10` predicts

\[
0.00117011487115.
\]

The B4-aware value is now frozen for the future `N=11` clock test.

---

## 7. Interpretation

The result supports the following finite-`N` picture:

\[
\rho_N
\sim
\text{finite-size prefactor}\times e^{-B_4N},
\]

with a non-negligible correction to the local exponent over the currently accessible range.

This is consistent with the earlier mapped communication barrier interpretation, but it does **not** prove that `B4` is the globally exhaustive asymptotic barrier. That remains a separate potential-theory/saddle-exhaustion problem.

The improvement is useful operationally because task 317 showed that transition/clock transfer has become the second-largest process error after preparation transfer. A better clock law directly addresses that emerging bottleneck.

---

## Verdict

\[
\boxed{\textbf{318 PASS}}
\]

for selecting and freezing a more predictive finite-`N` clock law.

The next direct microscopic holdout must not refit this clock after opening `N=11`.
