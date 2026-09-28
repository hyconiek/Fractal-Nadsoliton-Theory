# HISTORY-UNIFORM-LATE-BURN-CLOSURE-310
## Separating history closure from cross-N initialization shows that no extra memory coordinate is required to pass the predeclared 5% late-burn criterion

Date: 2026-09-27

Status:
- exact finite-state calculations for N=3,4,5 where explicitly replayed;
- reversible rank-120 spectral calculations for the larger finite systems, with convergence/error diagnostics;
- grouped leave-one-N-out validation for optional low-rank residual models;
- no N->infinity theorem and no claim of a fundamental Markov law.

## 0. Question inherited from 309

Report 309 found that the full cross-N contract at the original burn time can have a large strict conditional error on rare histories. The exploratory scan suggested that a later burn could help, but that scan mixed two effects:

1. genuine history dependence of the microscopic projection;
2. error in transporting the effective initial state / preparation law across N.

Task 310 separates them before deciding whether an auxiliary memory state is necessary.

The predeclared history-specific target is

    sup_h TV[
      P_micro(Y_2 | h),
      P_eff(Y_2 | h)
    ] <= 0.05.

For the one-step late-burn test, h is the first observed history state.

The effective transition is tested with:
- the SAME N effective D12 generator;
- the microscopic first latent 12-state marginal supplied exactly/within the spectral replay;
- the same seven-class observable and eta=0.05 readout channel.

Thus this test isolates history/closure error instead of cross-N initialization error.

## 1. Training-only burn selection

Candidate burn grid begins

    8, 12, 16, 24, ...

N=3 is decisive:

    t=16:
      max conditional TV = 0.0509817554  FAIL

    t=24:
      max conditional TV = 0.0470209805  PASS.

At t=24 the remaining training systems give

    N=4: 0.0283650634
    N=5: 0.0356488208
    N=6: 0.0289939497.

Therefore, without using N=7,

    boxed:
    t_burn = 24

is the earliest tested burn satisfying the 5% HISTORY-SUP criterion for every N=3..6.

### Important distinction

A stronger compound criterion that additionally demanded very small average/joint reduction error would NOT pass merely by waiting. For example the N=3 pair/weighted defect approaches about 4.28% rather than zero.

This is not evidence for persistent history dependence. It is a finite-N model-reduction floor. Task 310 therefore does not relabel that floor as memory.

## 2. Held-out N=7 history closure

After fixing t_burn=24 on N=3..6, test N=7 without changing the threshold.

To avoid reusing the integer-kappa points from earlier work, the heldout preparation set is

    kappa = 0.5, 1.0, 1.5, ..., 11.5

for the same pinning family.

Rank-120 reversible spectral evaluation gives

    max pair TV
      = 0.00764583359

    max probability-weighted conditional TV
      = 0.00764583359

    max strict conditional TV
      = 0.02532984347.

The conservative pair-level spectral-tail certificate is

    2.51e-10,

and rank-80/rank-120 values agree to displayed precision for the relevant defect.

Thus

    boxed:
    held-out N=7 strict history defect = 2.533% < 5%.

So the training-selected late-burn closure PASSES.

## 3. Why the earlier ~10% full-contract defect was not a memory no-go

The earlier exploratory full cross-N scan propagated:
- a predicted cross-N clock;
- a predicted cross-N slow shape;
- a predicted initial-slip/preparation state;
- the readout posterior;
- and the transition law
all at once.

Once the first latent marginal and N-specific effective generator are matched, the N=7 strict conditional defect at the same t=24 is only 2.53%.

Therefore the large full-contract supremum cannot be interpreted as evidence that a permanent extra memory variable is required.

It is predominantly a COMPOSED prediction error, especially sensitivity of the readout-conditioned posterior to the cross-N initial-state prediction.

## 4. Exact slow spectrum does not remove the residual

To test whether the remaining 2-5% defect is merely an inaccurate MZ slow rate, the six microscopic slow D12 sectors were extracted directly from the reversible full generator.

For N=6:

    lambda_1 = -0.0312431852610
    lambda_2 = -0.0349610431770
    lambda_3 = -0.0213246614497
    lambda_4 = -0.0226112751871
    lambda_5 = -0.0264716007243
    lambda_6 = -0.0262621436215.

This reproduces the earlier independent exact-mode values.

For N=7:

    lambda_1 = -0.0176858194262
    lambda_2 = -0.0200174644093
    lambda_3 = -0.0118406892882
    lambda_4 = -0.0126403864327
    lambda_5 = -0.0151048054050
    lambda_6 = -0.0149594689094.

Resolved overlaps with the corresponding basin Fourier sectors are about 97-98.5% at N=7.

Inverting the D12 Fourier map yields positive shell rates for N=3..7.

Replacing the MZ Q12 by the D12 generator with these EXACT slow eigenvalues does NOT materially improve history closure:

    N=7 MZ-Q12 strict conditional sup:
      0.0253298435

    N=7 exact-slow-spectrum Q12:
      0.0253716558.

Hence the residual is not explained by a slightly wrong slow clock/shell spectrum.

## 5. State-only scalar correction is rejected

There is also an exact structural obstruction.

With full latent 12-label observation and a fixed Markov Q12,
changing only the initial distribution p(J,t_burn) cannot change

    P_eff(J_2 | J_1).

Therefore a residual in the conditional transition cannot in general be repaired by a scalar correction only to the initial visible state.

The empirical low-dimensional alpha_k residual confirms the warning.

The six-dimensional residual of the existing preparation map has first-PCA energy

    87.9522%.

A one-scalar correction improves the N=7 alpha residual, but GROUPED leave-one-N-out fails the uniform criterion: on holdout N=3, the worst single alpha-component error changes from

    0.0428532  ->  0.0513684.

So this correction is rejected as a general law despite looking useful on N=7.

## 6. The conditional residual itself is almost one-dimensional

Now flatten the 7x7 conditional-transition residual

    Delta C = C_micro - C_same-N-Q12

for all N=3..6 and kappa=0..12 at t_burn=24.

Uncentered SVD about the physically meaningful zero-residual origin gives

    first singular direction energy fraction
      = 0.97517345.

So 97.52% of the training residual energy is one-dimensional.

The leading singular values are approximately

    0.615171,
    0.0966395,
    0.0164571,
    0.0048796,
    0.0005737, ...

This motivates an OPTIONAL augmented kernel

    C_aug = C_eff + z V,

where V is one fixed row-sum-zero residual direction.

Grouped cross-validation selects a cubic preparation law depending only on

    m = E[n_0/N],

not explicitly on N:

    z(m)
      = -13.4064192
        +37.9236865 m
        -35.4623584 m^2
        +10.8577075 m^3.

No clipping was needed in the validated range; the corrected rows remained stochastic numerically.

### Grouped leave-one-N-out strict-sup results

    hold N=3:
      baseline 4.7022%
      augmented 3.4046%

    hold N=4:
      baseline 2.8365%
      augmented 1.7308%

    hold N=5:
      baseline 3.5649%
      augmented 1.2505%

    hold N=6:
      baseline 2.8994%
      augmented 1.8404%.

Heldout N=7 on the half-grid gives

    baseline strict sup
      = 2.53298%

    augmented strict sup
      = 1.43215%.

Thus the one-coordinate residual GENERALIZES for the minimax history metric.

## 7. But the scalar refinement is not Pareto-dominant

The same N=7 correction changes the average/joint defect from

    0.76458%

to

    1.27403%.

N=6 shows the same trade-off in grouped validation.

Therefore the one-coordinate augmentation improves the WORST history but is not a uniformly better probabilistic model under every loss function.

Task 310 does NOT promote it to a new mandatory physical degree of freedom.

The canonical result remains the simpler statement:

    same-N late-burn history closure already passes the 5% threshold.

The rank-one object is retained as a diagnostic/minimax refinement.

## 8. Time behavior of the residual

Using the fixed training direction V on N=7:

    t=24:
      residual beta RMS ~0.0564
      direction energy fraction ~95.5%
      strict sup ~2.533%

    t~63.55:
      beta RMS ~0.0318
      direction energy fraction ~97.8%
      strict sup ~1.397%

    t~103.10:
      beta RMS ~0.0256
      direction energy fraction ~99.2%
      strict sup ~1.085%.

The residual becomes smaller and increasingly one-dimensional, but the scalar amplitude does not vanish over this finite window.

Because substituting the exact slow spectrum does not remove it, this is best described as a finite-N conditional reduction residual / hidden-state memory tail, not merely a rate error.

No asymptotic claim is made from these three times.

## 9. Verdict

### PASS — history-uniform late-burn closure

On the predeclared training grid, t_burn=24 is selected without N=7 and the heldout N=7 strict conditional defect is only 2.53%.

Therefore report 309's concern does NOT currently force an additional memory coordinate.

### FAIL — one-scalar correction of the initial state as a universal law

It improves N=7 but fails grouped leave-one-N-out on N=3.

### PASS as an optional minimax diagnostic — rank-one transition residual

The conditional residual is 97.52% one-dimensional on training data, and one scalar lowers strict-sup error on every grouped holdout and on unseen N=7 kappa half-points.

But it worsens the average/joint error in some cases, so it is not promoted to the canonical model.

## 10. What this changes conceptually

The effective architecture should now be stated as

    microscopic FIN
      -> controlled preparation
      -> short initial-slip layer
      -> late-burn visible state
      -> approximately history-uniform Q12 transition
      -> finite-window process bound.

The important correction is epistemic:

    cross-N preparation error != microscopic memory.

Those two errors must remain separated in all future process claims.

## 11. Next P0

### 311 — POSTERIOR-ROBUST-CROSS-N-CLOSURE

The remaining bottleneck is not same-N memory closure.
It is the composition

    uncertain cross-N effective prior
      + readout noise
      -> posterior over latent J
      -> conditional next-step prediction.

Task 311 should therefore:

1. carry the full 306 preparation-prediction set, not only a point estimate, to t_burn=24;
2. propagate that set through the readout/Bayes posterior;
3. derive a uniform/set-valued bound on P(Y_2|Y_1) for every allowed preparation;
4. combine it with the 310 same-N history defect and the 309 chain theorem;
5. test N=7 only after the bound is frozen on N=3..6.

Kill test:
if posterior amplification makes the cross-N prediction sets overlap the comparator, do NOT call the process discriminator certified. Widen the effective initial-state set or improve the cross-N preparation law; do not relabel this as fundamental memory unless the matched-prior same-N closure itself fails.
