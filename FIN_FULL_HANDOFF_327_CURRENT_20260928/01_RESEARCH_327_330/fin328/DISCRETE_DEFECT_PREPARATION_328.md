# FIN 328 — DISCRETE-DEFECT-PREPARATION
## Exact few-defect representation and a direct TV truncation certificate for finite-N preparation

Date: 2026-09-28

Status: **PASS for the declared pinning family over the tested N=3..12 range; N=11,12 tails are certified intervals rather than exact basin totals.**

## 1. Motivation

At the currently studied sizes, the localized phase is not well described as a uniformly populated Gaussian cloud in all eleven simplex directions. The dominant label carries about 98% of the mean-field mass, so configurations with only a few non-dominant copies are a natural finite-N coordinate system.

Define

`D = N - n0`,

and defect counts

`m=(m1,...,m11)`, with `sum_a m_a=D` and `n0=N-D`.

## 2. Exact defect weight formula

For the strict FIN matrix `A=A7`, define

`delta_a = A_00 - A_0a`,

`B_ab = A_ab - A_a0 - A_0b + A_00`.

For the pinning family

`mu_{N,kappa}(n) proportional to pi_N(n) 1_{J=0}(n) exp(kappa n0/N)`,

the relative weight of a defect configuration with respect to the all-zero-label state is exactly

`w_{N,kappa}(m)/w_{N,kappa}(0)`

` = N! / ((N-D)! prod_a m_a!)`

`   * exp[-g sum_a delta_a m_a + (g/(2N)) m^T B m - (kappa/N) D]`.

This is an algebraic rewriting of the original finite-N Gibbs weight. No Gaussian approximation and no fitted PCA coefficient enters the formula.

The basin condition `J=0` remains part of the preparation definition and is not replaced by the condition `D<=K`.

## 3. Direct validation against full microscopic weights

The exact formula was checked against every available microscopic state labelled `J=0` for `N=3,...,10`.

Maximum absolute discrepancies in log relative weight are of order `10^-14`:

- N=3: `3.55e-15`;
- N=4: `5.33e-15`;
- N=5: `7.11e-15`;
- N=6: `8.88e-15`;
- N=7: `1.78e-14`;
- N=8: `2.13e-14`;
- N=9: `1.95e-14`;
- N=10: `2.84e-14`.

These are floating-point replay residuals of an exact algebraic identity.

## 4. Exact monotonicity of the truncation tail in kappa

Because

`mu_{N,kappa}/mu_{N,0} proportional to exp[-(kappa/N)D]`,

for every cutoff `K`,

`d/dkappa P_kappa(D>K)`

` = -(1/N) Cov_kappa(1_{D>K}, D) <= 0`.

The covariance is nonnegative because both functions are nondecreasing functions of the scalar random variable `D`.

Therefore:

**For every `kappa>=0`, the worst truncation tail occurs exactly at `kappa=0`.**

No scan over pinning strength is required to certify the entire positive-pinning class.

## 5. Exact TV meaning of truncation

Let `mu^(K)` be `mu` conditioned on `D<=K`. Then exactly

`TV(mu,mu^(K)) = mu(D>K)`.

For any subsequent linear Markov/observation channel `Kobs`, data processing gives

`TV(mu Kobs, mu^(K) Kobs) <= mu(D>K)`.

Thus the omitted defect mass is immediately a bound on later observed-process error, provided the later workflow keeps all outcomes. If it postselects on localization again, the conditioning step requires its own rejection/error budget.

## 6. Tail results

At `kappa=0`, the maximum tails over the studied range are:

| cutoff | number of retained defect configurations | maximum upper tail over N=3..12 |
|---:|---:|---:|
| `D<=4` | 1,365 | `0.00966504` = 0.9665% |
| `D<=5` | 4,368 | `0.00385942` = 0.3859% |
| `D<=6` | 12,376 | `0.00160280` = 0.1603% |

For N=7..10 the `D>6` tails are known directly from the fully labelled basin:

- N=7: `4.36369e-05`;
- N=8: `2.30108e-04`;
- N=9: `8.50733e-04`;
- N=10: `1.43737e-03`.

For N=11 and N=12, not all basin labels were retained exactly in the earlier cache, so the result is a certified interval obtained from exact low-defect classification plus the previously bounded unclassified stationary mass:

- N=11: `P(D>6 | J=0) in [0.00098792, 0.00160280]`;
- N=12: `P(D>6 | J=0) in [0.00108593, 0.00143101]`.

The quoted unified 0.1603% bound uses the upper endpoints.

## 7. Compression at N=12

The full N=12 count-state space has

`1,352,078` states.

The complete defect catalogue with `D<=6` has only

`12,376` configurations.

Thus the preparation representation is over 100 times smaller while retaining an explicit worst-case preparation-TV error below 0.161% throughout the currently studied N=3..12 range and for every `kappa>=0`.

## 8. Scientific meaning

This improves on the earlier empirical preparation residual models in three ways:

1. probabilities are positive by construction;
2. finite-population combinatorics and defect interactions are exact;
3. the approximation has a direct omitted-mass certificate rather than a fitted cross-validation envelope.

It does **not** yet prove that `D<=6` remains sufficient as `N->infinity`. The current result is a finite-N controlled representation for N<=12. In fact the few-defect picture is expected eventually to cross over to a broader fluctuation regime as N grows.

## 9. Next use

The next important test is not another PCA fit. It is to replace the preparation component of the observed-process predictor by the exact `D<=6` defect distribution and measure whether the empirical `eps_prep` in reports 313/320/322 is reduced without introducing new fitted parameters.
