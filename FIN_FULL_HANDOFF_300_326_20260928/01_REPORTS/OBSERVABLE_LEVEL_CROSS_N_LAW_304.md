# OBSERVABLE-LEVEL-CROSS-N-LAW-304
## A q_d-free clock/shape law predicts held-out N=7 microscopic observables and N=8 effective observables without per-N six-rate refitting

Date: 2026-09-27

Status:
- exact algebra for D12-circulant reconstruction from dimensionless Fourier ratios;
- exact microscopic slow-eigenvalue inputs for N=3..8 from the established leave-one-out lane;
- direct microscopic preparation/histogram replay for N=3..7;
- frozen cross-N fits use only N=3..6 and are tested on N=7,8;
- direct N=8 microscopic replay was attempted but did not complete in the bounded current run, so N=8 histogram validation is effective-Q12 only;
- no N->infinity theorem and no claim that the empirical finite-N fits are fundamental FIN laws.

## 1. Why 304 is different from 303

Task 303 showed that a simple shell-wise extrapolation could predict the coarse clock and selected histograms even while one rare shell rate q1 was badly predicted.

That left a conceptual weakness: the cross-N description still passed through six q_d coordinates.

Task 304 removes q_d from the predictive interface.

The cross-N prediction is split into:

1. a **clock law** for the slow microscopic rate rho_N;
2. a **shape law** for dimensionless Fourier ratios R_k;
3. an **early-memory fingerprint**;
4. an explicit **initial-slip / preparation uncertainty**.

No q_d(N=7) or q_d(N=8) is used to construct the held-out clock or shape prediction. The old q_d are used only afterward as validation data for the effective layer.

## 2. Exact clock/shape separation theorem

For any reflection-symmetric D12-circulant 12-state generator, write its nontrivial Fourier eigenvalues as

    lambda_k = -rho R_k,
    k=1,...,6.

The exact Z3 identity gives

    R_4 = 1.

Starting from localized state J=0, at dimensionless time

    tau = rho t,

the probability of state j is exactly

    p_j(tau)
      = 1/12 [
          1
          + 2 sum_{k=1}^5 exp(-R_k tau) cos(2 pi k j/12)
          + exp(-R_6 tau)(-1)^j
        ].

Therefore:

    boxed:
    at fixed tau, the full one-time shape is independent of rho.

The late protocol can be decomposed exactly into:

    physical measurement time
      <- rho_N

and

    dimensionless histogram shape
      <- (R_1,R_2,R_3,R_5,R_6), with R_4=1.

The seven-bin Y=cos(2 pi J/12) histogram is then only a fixed linear push-forward of this p_j, followed by the declared readout channel.

This is an exact algebraic statement inside the D12 effective class, not a numerical fit.

## 3. Clock law from microscopic slow modes only

Do not fit the clock from holdout q_d.

Use the exact microscopic slow k=4 eigenvalues:

    N=3: 0.131439786195646
    N=4: 0.071692277481914
    N=5: 0.0401867559160965
    N=6: 0.0226112751871144

and freeze the N=3..6 law

    log rho_N
      = -0.28041005 - 0.58591460 N.

Held-out predictions:

### N=7

    predicted rho_7 = 0.0125026797792
    exact microscopic rho_7 = 0.0126403864327
    relative error = 1.0894 %.

The frozen design time is

    t_pred = 0.5427059873/rho_pred
           = 43.4071732530.

### N=8

    predicted rho_8 = 0.00695894860050
    exact microscopic rho_8 = 0.00698967186661
    relative error = 0.4396 %.

Thus a two-parameter law frozen at N<=6 predicts the next two microscopic slow clocks to about 1.1% and 0.44%.

This fit is descriptive finite-N evidence. Its exponent must NOT be identified with an asymptotic barrier without a separate theorem.

## 4. Observable shape law without shell rates

For N=3..6 compute the already established dimensionless effective ratios

    R_k = -lambda_k/rho,

for k=1,2,3,5,6.

Freeze the affine observable law

    R_k(N) = a_k + b_k N.

Coefficients:

    R1 = 1.3438224018 + 0.0058681230 N
    R2 = 1.4186939697 + 0.0206383308 N
    R3 = 0.9875059308 - 0.0074722558 N
    R5 = 1.0354972669 + 0.0225422773 N
    R6 = 1.0403786645 + 0.0201530933 N.

The training R-vector trajectory is almost one-dimensional:

    first PCA variance fraction = 0.988016.

So about 98.8% of the N=3..6 variation lies along one affine direction in the five-dimensional shape space.

This is a useful compression: N selects one point on a fixed observable trajectory instead of introducing six new shell rates.

## 5. Held-out shape predictions

### N=7

Maximum relative error among the five predicted R_k:

    1.3100 %

(the largest error is R2).

At the frozen dimensionless time tau*=0.5427059873 and 5% symmetric readout error:

    TV(predicted shape, actual Q12 shape)
      = 0.00313620
      = 0.3136 %.

### N=8

Maximum relative R_k error:

    2.6272 %.

Yet the reconstructed seven-bin histogram differs from the actual pre-existing Q12 histogram by only

    TV = 0.00623540
       = 0.6235 %.

Again, the selected observable is more stable than each internal coordinate separately.

## 6. Combined clock + shape holdout

Now evaluate the actual Q12 process at the **physical time predicted by the microscopic clock law**, rather than giving the predictor the true holdout clock.

### N=7

    TV(predicted observable law,
       actual Q12 at predicted physical time)
      = 0.00584265
      = 0.5843 %.

### N=8

    TV = 0.00724125
       = 0.7241 %.

Thus the frozen N<=6 clock/shape law predicts the held-out effective one-time distributions at below 0.75% TV for both N=7 and N=8.

No q_d holdout values are used to construct these predictions.

## 7. Direct microscopic N=7 validation

The N=7 microscopic process was rerun at the time predicted solely from the N<=6 microscopic clock law:

    t_pred = 43.4071732530.

Contract:
- exact leave-one-out Gibbs generator;
- exact D12-equivariant basin definition;
- equilibrium Gibbs conditioned on localized basin J=0;
- same 5% readout channel;
- same seven-bin observable.

The cross-N observable law gives

    boxed:
    TV(predicted observable law,
       direct microscopic N=7 histogram)
      = 0.0130691
      = 1.3069 %.

For comparison, the pre-existing actual N=7 Q12 at that same time differs from the direct microscopic histogram by

    0.8421 % TV.

So the new total error naturally decomposes into:

    cross-N clock/shape error
      +
    microscopic-to-Q12 reduction error.

This is the desired controlled-theory architecture.

## 8. Model-discrimination power survives the cross-N prediction

At the same predicted N=7 time, against the declared report-293 same-rho comparator:

Direct microscopic FIN:

    TV = 0.0756751
    Chernoff C = 0.00794230
    sufficient Pe<=5% bound = 290 trials.

Cross-N observable-law prediction:

    Chernoff C = 0.00754846
    sufficient Pe<=5% bound = 306 trials.

So replacing the exact held-out FIN distribution by the frozen cross-N observable law changes the conservative trial count by only 16 trials.

This is a stronger cross-N result than 303 because the predictor uses neither q_d(N=7) nor the true N=7 clock.

## 9. Early-memory fingerprint is approximately N-stable

Direct microscopic calculations under the same equilibrium-J0 preparation give

    F_N = Delta s/rho_N

with

    N=3: 0.4217603
    N=4: 0.3834761
    N=5: 0.4392667
    N=6: 0.4010595.

Freeze the simplest possible N<=6 predictor:

    F = mean = 0.4113907.

Held-out N=7 gives

    F_7 = 0.4135985,

only

    0.5338 %

away from the frozen constant predictor and inside the complete N=3..6 range

    [0.383476, 0.439267].

This is surprisingly stable finite-N evidence that the **dimensionless early memory layer** has a reproducible scale relative to rho.

No direct N=8 microscopic F_8 is claimed in 304.

## 10. Initial slip does NOT yet have an equally good cross-N point law

Exact slow-mode deficits are

    1-Z_N:
      N=3: 8.9146 %
      N=4: 7.6001 %
      N=5: 4.1622 %
      N=6: 2.7948 %
      N=7: 1.4839 %
      N=8: 0.8647 %.

A log-linear point fit frozen on N=3..6 overpredicts the heldouts by

    N=7: 28.7 % relative;
    N=8: 46.8 % relative.

Therefore task 304 does **not** promote a precise cross-N slip law.

A training-calibrated conservative exponential envelope does cover both holdouts:

    N=7 upper bound = 2.2335 % > 1.4839 % actual;
    N=8 upper bound = 1.4850 % > 0.8647 % actual.

More importantly, when the local per-N first memory moment is available,

    Z_M1 = 1/(1+M1)

predicts the exact deficit much more accurately:

    N=7 relative error = 1.61 %;
    N=8 relative error = 0.824 %.

So the unresolved task is not the algebraic slip correction. It is a cross-N law or uniform bound for the microscopic memory quantity that produces it.

## 11. Preparation uncertainty is now the dominant operational systematic

For the clock-predicted N=7 measurement time:

    TV(deep-seed preparation,
       equilibrium-J0 preparation)
      = 0.0235427
      = 2.3543 %.

The direct FIN-vs-comparator signal is

    7.5675 % TV.

Therefore this one preparation ambiguity alone is

    boxed:
    31.11 % of the model-discrimination signal.

It is also larger than the complete 1.31% error of the cross-N observable-law prediction against direct microscopic N=7.

Thus further fine tuning of the R_k extrapolation is no longer the highest-value task.

The bottleneck has moved to the operational preparation map.

## 12. Main verdict

### PASS — observable-level cross-N law

A frozen law calibrated only on N=3..6 predicts:

- microscopic slow clocks at N=7,8 to about 1.09% and 0.44%;
- effective held-out histogram shapes at 0.31% and 0.62% TV;
- combined clock+shape effective histograms at 0.58% and 0.72% TV;
- the direct microscopic N=7 histogram at 1.31% TV;
- the N=7 early-memory fingerprint to 0.53% using only a constant N<=6 predictor.

This is substantially stronger than fitting six fresh q_d values at each N.

### PARTIAL — process-level closure

Initial slip is controllable but its cross-N point law is not yet accurate enough.

Preparation ambiguity is already larger than the cross-N prediction error.

Direct microscopic N=8 validation remains uncompleted in this task.

## 13. Scientific meaning

The current effective lane can now be organized as

    microscopic FIN
      -> clock rho_N
      + observable shape R_k(N)
      + early-memory correction
      + preparation map
      -> predicted histogram.

This is a more physical interface than the six shell rates themselves because every retained quantity is tied to an identifiable response.

But the fit laws in N remain empirical finite-N laws. They have not been derived asymptotically from the FIN potential/barrier structure.

## 14. Next P0

The next task should be

    305 — PREPARATION-CONTRACT-PREDICTION-SET.

Goal:
construct an operationally explicit class of microscopic preparations and propagate it to a prediction set, then test whether the resulting uncertainty remains strictly smaller than the FIN-vs-comparator separation.

Required tests:
1. quantify contraction of preparation uncertainty before the late measurement window;
2. distinguish freely evolving preparation from any restrained/pinned preparation protocol;
3. include the controller used for preparation in the resource/accounting ledger;
4. derive a set-valued seven-bin prediction rather than silently choosing equilibrium-J0;
5. kill the discriminator if the allowed preparation set overlaps the comparator prediction class.

A parallel mathematical task remains to derive a cross-N bound for M1/slip rather than fit it pointwise.

## Epistemic boundary

304 is a finite-N controlled predictive result.

It does NOT establish:
- an asymptotic N->infinity law;
- a physical apparatus realization;
- a fundamental source for A7, g or the update controller;
- physical spatial incidence;
- quantum mechanics, general relativity or a fundamental ontology.
