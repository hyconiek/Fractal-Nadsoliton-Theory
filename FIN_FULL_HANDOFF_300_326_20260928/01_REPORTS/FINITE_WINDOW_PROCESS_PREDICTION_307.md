# 307 — FINITE-WINDOW-PROCESS-PREDICTION
## Held-out two-time laws after the preparation/memory layer

Date: 2026-09-27

Status:
- exact finite-state microscopic replay on the declared finite-N leave-one-out Gibbs chain for N=3..7;
- cross-N clock, D12 shape law, and preparation map are frozen from tasks 304-306;
- N=7 two-time data are held out from all fits;
- seven-bin laws are conditioned on the system being in one of the 12 localized basins at both readout times;
- no claim of a full arbitrary-multi-time Markov theorem is made.

## 1. Question

Tasks 304-306 showed that a frozen effective law can predict one-time seven-bin
histograms across N after a controlled preparation map. This is insufficient for
a process-level claim: two models can reproduce every one-time marginal and
still disagree on temporal dependence.

Task 307 therefore asks whether the same frozen effective dynamics predicts

    P(Y_t1, Y_t2)

without fitting any N=7 transition parameter to N=7 two-time data.

## 2. Frozen protocol

No kinetic parameter is retuned in this task.

Frozen ingredients:

- exact microscopic leave-one-out Gibbs law;
- full-X7 basin map;
- thermal-pin preparation family 0 <= kappa <= 12;
- preparation map alpha_k(N,m) from 306;
- cross-N clock rho_N from 304;
- cross-N shape R_k(N) from 304;
- seven reflection-orbit readout bins;
- symmetric 12-label readout error, nominal eta=0.05.

The burn-in is fixed at

    t_b = 12

microscopic time units, after the fast memory layer used to define the slip
coordinates in 306.

Two dimensionless post-burn intervals were declared before held-out evaluation:

    Delta tau = rho_pred * (t2-t1) in {0.25, 0.50}.

For N=7,

    rho_pred = 0.012502679779153265.

## 3. Exact microscopic two-time law

Let F_a(x) be the indicator that microscopic state x belongs to localized
seven-bin class a. For preparation mu the unnormalized exact two-time law is

    J_ab(t1,t2)
      = mu exp(Q t1) diag(F_a) exp(Q (t2-t1)) F_b.

The calculation is performed exactly in finite count-state space. Residual
nonlocalized states are not silently assigned to a localized bin. Instead the
7x7 table is conditioned on localized readout at both endpoints.

The readout channel is then applied independently at both times:

    J_obs = C_eta^T J C_eta.

This defines the microscopic target used for all errors below.

## 4. Effective prediction

The predicted first-time 12-state distribution is reconstructed from the
frozen preparation coordinates:

    p_j(t1)
      = 1/12 [
          1
          + 2 sum_(k=1)^5 alpha_k exp(-rho R_k t1)
                cos(2 pi k j/12)
          + alpha_6 exp(-rho R_6 t1)(-1)^j
        ].

The same frozen eigenvalues define the 12-state semigroup over Delta t.
Hence

    J_eff(j1,j2)
      = p_j1(t1) [exp(Q12 Delta t)]_(j1,j2).

It is then folded to seven bins and passed through the identical readout
channel.

Thus the second observation is a genuine process prediction. It is not fitted
from the second-time histogram.

## 5. Why joint-law error is stronger than marginal error

For a joint law J let p1 be its first marginal and K(.|a) its conditional
second-time law.

A useful diagnostic separates:

1. marginal error at t1;
2. marginal error at t2;
3. conditional/transition error

       E_cond
         = sum_a p1_exact(a)
             TV(K_exact(.|a), K_pred(.|a)).

A model may have small one-time marginal errors and still fail through E_cond.
That is exactly the failure mode task 307 was designed to detect.

## 6. Training-only leave-one-N-out envelope

For each holdout N in {3,4,5,6}, clock + R_k + preparation coordinates were
predicted from the remaining training Ns. No two-time parameter was fitted.

### Delta tau = 0.25

Maximum joint-TV error over 0<=kappa<=12:

    N=3:  0.06431747
    N=4:  0.03388788
    N=5:  0.02782034
    N=6:  0.02336960

Frozen per-window envelope:

    B_0.25 = 0.06431747.

### Delta tau = 0.50

    N=3:  0.05129077
    N=4:  0.02475222
    N=5:  0.01949304
    N=6:  0.02092817

Frozen per-window envelope:

    B_0.50 = 0.05129077.

Across both windows the global training maxima are

    joint TV       <= 0.06431747
    conditional TV <= 0.06408385
    first marginal <= 0.01683840
    second marginal<= 0.01711985.

This is the central control. The joint-law uncertainty is several times larger
than the one-time marginal uncertainty, so reusing a one-time error budget
would have been invalid.

## 7. Held-out N=7

N=7 joint data were evaluated only after the training envelopes above were
fixed.

### Window A: Delta tau = 0.25

Worst error over the declared preparation grid:

    joint TV       = 0.02139296
    first marginal = 0.01666777
    second marginal= 0.01396472
    conditional TV = 0.01803549.

Therefore the held-out joint error lies below the frozen training envelope by

    0.06431747 - 0.02139296
      = 0.04292451.

So process prediction PASSES.

### Window B: Delta tau = 0.50

    joint TV       = 0.01994071
    first marginal = 0.01666779
    second marginal= 0.01154609
    conditional TV = 0.01581224.

Held-out slack relative to its training-only window envelope is

    0.05129077 - 0.01994071
      = 0.03135006.

Again, process prediction PASSES.

A useful feature is that conditional error is not larger than the full joint
error on held-out N=7. The effective process is therefore not passing merely
because its marginals are accurate.

## 8. Stronger comparator test: match the first-time state

To avoid confusing preparation discrimination with transition discrimination,
the report-293 same-rho/same-total-exit comparator is initialized at the SAME
predicted 12-state distribution at t1 as FIN.

Thus the comparator differs only in its transition semigroup over the second
part of the window.

At nominal eta=0.05:

### Delta tau = 0.25

Minimum predicted FIN/comparator joint separation:

    0.04294081 TV.

This is smaller than the frozen training envelope

    B_0.25 = 0.06431747.

Therefore the short window is NOT pre-certified as a model discriminator,
even though held-out process prediction itself passes.

### Delta tau = 0.50

Minimum separation:

    0.06338893 TV.

Now

    0.06338893 - B_0.50
      = 0.01209816 > 0.

Hence the longer, predeclared window has a positive training-only
model-discrimination margin.

This distinction is important:

    Delta tau=0.25:
      process prediction PASS,
      pre-certified comparator discrimination FAIL;

    Delta tau=0.50:
      process prediction PASS,
      pre-certified comparator discrimination PASS.

No N=7 joint datum was used to choose either window; both were declared before
held-out evaluation.

## 9. Direct held-out triangle bounds

The held-out micro/prediction error also gives a post-validation lower bound on
the true microscopic distance from the matched comparator.

By the triangle inequality:

### Delta tau=0.25

    TV(micro, comparator)
      >= 0.04294081 - 0.02139296
      = 0.02154785.

### Delta tau=0.50

    TV(micro, comparator)
      >= 0.06338893 - 0.01994071
      = 0.04344822.

These are held-out validation statements, not substitutes for the training-only
certification above.

## 10. Readout and clock stress

A training-only stress replay was performed at clock scales

    0.98 and 1.02

and readout errors

    eta=0 and eta=0.10.

For Delta tau=0.50 the worst training joint envelopes and corresponding
predicted N=7 comparator separations remain ordered with positive margins:

    clock 0.98, eta=0:
      B_train = 0.05783566
      separation = 0.06981776
      margin = 0.01198210

    clock 0.98, eta=0.10:
      B_train = 0.04618298
      separation = 0.05706549
      margin = 0.01088252

    clock 1.02, eta=0:
      B_train = 0.05674392
      separation = 0.07081280
      margin = 0.01406888

    clock 1.02, eta=0.10:
      B_train = 0.04528475
      separation = 0.05780702
      margin = 0.01252227.

This is a training-frozen robustness certificate for the declared long window.
A full direct N=7 replay at every stress corner was not completed in task 307,
so do not describe it as a direct held-out stress theorem.

## 11. Nominal trajectory-pair sample complexity

Treat one independent realization of the two-time measurement as one 49-cell
categorical sample.

For the matched-first-marginal comparator at eta=0.05:

    Delta tau=0.25:
      min Chernoff information = 0.00464308
      sufficient 5% bound: 496 independent trajectory pairs;

    Delta tau=0.50:
      min Chernoff information = 0.00624463
      sufficient 5% bound: 369 independent trajectory pairs.

These counts are nominal model-vs-model values. They do not themselves absorb
the full model-reduction envelope.

## 12. Main result

Task 307 passes the process-level kill test at the two-time level:

    one microscopic law
      -> controlled preparation/slip map
      -> frozen cross-N effective generator
      -> held-out two-time joint law.

The key quantitative result is that held-out N=7 conditional-transition error
is around 1.6-1.8%, while the predeclared Delta tau=0.50 window retains a
positive training-only margin against the matched comparator.

This is stronger than one-time histogram matching.

## 13. Boundaries

Task 307 does NOT prove:

- exact microscopic Markovianity of the basin label;
- an arbitrary-time or arbitrary-preparation process theorem;
- three-time consistency;
- Chapman-Kolmogorov closure of the exact microscopic readout process;
- the same bounds for an unrestricted comparator class;
- physical realization of the FIN variables;
- fundamental space, QM, GR, or a source law for A7/g/update dynamics.

The exact microscopic projected process still contains hidden information. The
result is a controlled finite-window effective prediction after a declared
preparation/memory layer.

## Verdict

**307 PASS for held-out two-time process prediction.**

**307 PASS for pre-certified matched-comparator discrimination in the
predeclared Delta tau=0.50 window.**

**307 FAILS to establish a full multi-time Markov-process theorem; that is the
next task.**
