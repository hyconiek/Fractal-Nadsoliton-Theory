# HELD-OUT-N-MICROSCOPIC-PREDICTION-303
## Frozen N=6 protocol survives a direct N=7 microscopic histogram test; observable-level cross-N prediction is substantially better than shell-wise prediction

Date: 2026-09-27

Status:
- exact finite-state leave-one-out Gibbs computation at N=7;
- exact reconstruction of the same D12-equivariant basin definition used in reports 136/140 and task 300;
- no fitting to the new N=7 microscopic histograms;
- pre-existing q_d(N=7) from report 214 used only in Level 1 reduction validation;
- Level 2 cross-N fits are frozen on N=3..6 and evaluated at N=7, with a q_d-level secondary check at N=8;
- no asymptotic theorem and no claim that the empirical cross-N fit is a FIN-derived law.

## 1. Why this task was needed

The post-299 review required a prediction outside the calibration point. A weak version would merely compare a new histogram with a Q12 generator whose six rates had already been extracted from the same N. That is useful for checking the reduction, but it does not yet show a single cross-N law.

Task 303 therefore separates two levels:

1. **Level 1 — held-out histogram test of an already existing N=7 reduction.**
   The N=7 microscopic histogram is newly generated directly from the exact leave-one-out generator. No q_d is fitted to that histogram.

2. **Level 2 — cross-N holdout.**
   A deliberately simple law is fitted only on N=3,4,5,6 and then frozen before evaluation at N=7 (and secondarily N=8 at the q_d/Q12 level).

The second test is more relevant to the question whether FIN is becoming one predictive law rather than a new six-rate description at every N.

## 2. Frozen microscopic contract

Use exactly the contract from task 300, but set N=7:

- `g = 5.145228719489142`;
- exact leave-one-out Gibbs heat-bath generator;
- count-state space size `31824`;
- basin initialization

    theta0 = (g/N) X7^T n;

- descent of the same full-X7 dual potential V_g;
- D12-equivariant orbit propagation;
- microscopic preparation = equilibrium Gibbs measure conditioned on localized basin J=0;
- full 12 localized labels plus residual class are retained before compression;
- seven-bin observable Y=cos(2 pi J/12) is used only after reflection symmetry is checked;
- readout error eta=0.05;
- late design time fixed by the already frozen dimensionless value

    rho t = 0.5427059873.

No parameter in this contract was reoptimized against the N=7 histogram.

## 3. Reconstruction checks

The N=7 basin reconstruction gives:

    total localized mass
      = 0.9992530200759365

which matches the previously stored N=7 localized result.

The 12 localized basins are symmetry-equivalent. The exact stationary residual is

    ||pi Q||_infinity
      = 7.70e-17.

For the direct microscopic equilibrium-J0 preparation, the conditional 12-bin distribution is reflection symmetric to roundoff:

    max_j |p_j-p_{-j}|
      <= 2.64e-16

on the late test point.

So the seven-bin compression is valid for this declared contract.

## 4. Level 1: pre-existing Q12(N=7) predicts a newly generated microscopic histogram

The pre-existing report-214 rates are

    q1 = 0.00025888703775233833
    q2 = 0.0007158309169877905
    q3 = 0.002477122060163837
    q4 = 0.0022347762384702504
    q5 = 0.00100447616108276
    q6 = 0.0007488238862539219.

They imply

    rho_7 = 3(q1+q2+q4+q5)
          = 0.012641911062879417

and therefore the frozen late test time is

    t_* = 0.5427059873/rho_7
        = 42.92910973670377.

### 4.1 Microscopic-to-Q12 total-variation error

For the equilibrium-J0 preparation:

    t=0.5:  TV = 0.01117194
    t=1:    TV = 0.01371119
    t=2:    TV = 0.01537345
    t=t_*:  TV = 0.00856507.

The late 12-bin error is therefore

    boxed:
    TV(microscopic N=7, pre-existing Q12(N=7))
      = 0.8565 %.

This is a true held-out histogram test: the q_d were not fitted to these newly generated distributions.

### 4.2 Seven-bin discriminator against the report-293 same-rho comparator

At 5% symmetric readout error:

    TV(microscopic FIN, comparator)
      = 0.07525138;

    Chernoff C
      = 0.007947860686;

    sufficient bound for Pe <= 5%
      = 290 independent trials.

The effective Q12 model predicts:

    TV(Q12 FIN, comparator)
      = 0.07360730;

    Chernoff C
      = 0.007922463612;

    sufficient bound
      = 291 trials.

Thus the direct microscopic and effective protocols differ by only one trial in this conservative Chernoff-bound metric.

This strengthens the operational meaning of reports 297/302: the late discriminator is not an artifact of using Q12 in place of the microscopic law.

## 5. Early multi-time fingerprint at held-out N=7

Using the same k=2 mode and times 0.5, 1 and 2:

    C(0.5) = 0.972601823515
    C(1)   = 0.958966760917
    C(2)   = 0.937154701867.

Log-slopes are

    s(0.5,1) = -0.028236717443
    s(1,2)   = -0.023008042143.

Hence

    Delta s = 0.005228675299

and in units of the N=7 slow clock

    boxed:
    Delta s / rho_7 = 0.413598487873.

For N=6 task 302 gave

    Delta s/rho_6 = 0.401059531425.

The N=7 held-out value differs by only about 3.13%.

This is evidence that the dimensionless early-memory fingerprint is not peculiar to the N=6 calibration point. It is still finite-N evidence, not a universal constant.

## 6. Preparation uncertainty remains a dominant systematic

Even at N=7, different microscopic realizations of the same macro-label J=0 remain distinguishable.

At t=1:

    TV(deep seed, equilibrium-J0)
      = 0.02710;

    TV(flat-J0, equilibrium-J0)
      = 0.51396.

At the late test time t_*:

    TV(deep seed, equilibrium-J0)
      = 0.02372;

    TV(flat-J0, equilibrium-J0)
      = 0.31020.

So the task-300 kill-test survives N=7. The macro-label J=0 still does not define a unique preparation.

A usable experimental statement must therefore declare a preparation map or propagate this uncertainty into the prediction set.

## 7. Level 2: frozen cross-N prediction from N=3..6 to N=7

To avoid hiding six new free rates in each new N, fit the deliberately simple diagnostic law

    log q_d(N) = a_d + b_d N

using ONLY N=3,4,5,6.

Then evaluate it at N=7.

Shell-wise relative errors are:

    q1: 20.923 %
    q2:  1.211 %
    q3:  2.898 %
    q4:  2.613 %
    q5:  1.358 %
    q6:  0.0625 %.

So the law is NOT uniformly accurate at the level of all six shell rates.

In particular q1 is a clear failure of shell-wise extrapolation.

However the coarse clock is much more stable:

    rho_pred = 0.01256228513
    rho_true = 0.01264191106

so

    boxed:
    relative rho error = 0.6299 %.

At the frozen late test point the seven-bin histogram of this N<=6 extrapolation differs from the direct microscopic N=7 histogram by

    boxed:
    TV = 0.01288074

when evaluated with the actual N=7 test time.

If one compares only the predicted and actual Q12 distributions, the seven-bin TV is

    0.00558982.

Therefore a poor prediction of one rare shell channel does not automatically imply a poor prediction of the selected observable.

## 8. Why q1 can be wrong while the protocol is still accurate

A local sensitivity calculation of the frozen N=7 seven-bin protocol gives the TV derivative per unit logarithmic shell perturbation:

    d=1: 0.01279
    d=2: 0.03954
    d=3: 0.11336
    d=4: 0.10543
    d=5: 0.05049
    d=6: 0.01924.

Thus the chosen observable is almost nine times more sensitive to q3 than q1.

This resolves an apparent paradox:

- q1 can have a 20.9% extrapolation error;
- yet the seven-bin histogram remains close.

The experimental protocol is primarily testing the high-leverage q3/q4 structure, not every shell rate equally.

Therefore successful prediction of this histogram must NOT be promoted to a proof that all six q_d have been predicted accurately.

## 9. Second q_d-level holdout at N=8

Using the SAME N=3..6 log-linear fits, with no refit:

    max shell-wise error at N=8
      = 46.07 %
      (again q1);

but

    relative rho error
      = 0.6022 %;

and the seven-bin effective-Q12 histogram differs by only

    TV = 0.00686994.

This secondary check reinforces the same conclusion:

    coarse-clock / selected-observable prediction
      is much more stable than
    individual rare-shell prediction.

This N=8 statement is only at the pre-existing q_d/Q12 level; no new direct N=8 microscopic histogram is claimed here.

## 10. Barrier-ratio holdout

A more FIN-structured diagnostic uses the previously mapped saddle differences rather than six independent slopes.

Calibrate the N=6 observed ratios once and propagate N=6 -> 7 using

    q3/q4(N+1)
      ~ [q3/q4](N) exp(B4-B3),

    q4/q5(N+1)
      ~ [q4/q5](N) exp(B5-B4).

For N=7:

    q3/q4 predicted = 1.11824906
    q3/q4 actual    = 1.10844299
    relative error  = 0.885 %.

For q4/q5:

    predicted       = 2.36045596
    actual          = 2.22481760
    relative error  = 6.10 %.

Thus the mapped barrier ordering carries real predictive information, especially for q3/q4, but finite-N prefactors remain material for q4/q5.

This is consistent with report 215 and is not an Eyring-Kramers theorem.

## 11. Main verdict

Task 303 has a split verdict.

### PASS — held-out microscopic observable prediction

The frozen protocol survives a direct N=7 microscopic test:

- basin construction replays exactly;
- reflection compression remains valid;
- pre-existing Q12 predicts the new microscopic late histogram at <1% TV;
- direct-microscopic and Q12 Chernoff requirements are 290 vs 291 trials;
- the normalized early-memory fingerprint changes by only ~3.1% from N=6 to N=7.

### PARTIAL — cross-N law

A simple law frozen on N<=6 predicts:

- rho_7 to 0.63%;
- the selected N=7 observable to ~1.29% TV versus the microscopic histogram;
- but q1 only to ~21%.

At N=8 the coarse clock remains within ~0.60% while q1 misses by ~46%.

Therefore the current evidence supports an emerging law for selected coarse observables and clock ratios more strongly than a law for every individual q_d.

This distinction is essential.

## 12. New research consequence

The next task should not simply fit better curves to all six q_d.

The scientifically sharper target is:

    304 — OBSERVABLE-LEVEL-CROSS-N-LAW

Derive or constrain directly the N-dependence of:

1. rho_N;
2. the high-leverage dimensionless Fourier ratios R_k;
3. the early memory fingerprint Delta s/rho;
4. preparation-slip quantities;

using microscopic FIN structure, then test on held-out N without introducing six independent new shell rates.

A success criterion should be stated at the observable level first, with shell-wise reconstruction as a stronger optional theorem.

## Epistemic boundary

This is a controlled finite-N prediction result.

It does NOT establish:
- an N->infinity theorem;
- a physical apparatus realization;
- a fundamental source of A7, g, the update rule, or the event controller;
- physical space, QM, GR or a fundamental ontology.

It strengthens the claim that one fixed microscopic FIN process can produce nontrivial, testable predictions outside a single calibration point.
