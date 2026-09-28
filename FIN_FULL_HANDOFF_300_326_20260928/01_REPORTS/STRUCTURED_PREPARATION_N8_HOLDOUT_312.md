# 312 — STRUCTURED-PREPARATION-PREDICTION-SET + DIRECT-N8-HOLDOUT

Date: 2026-09-27

Status:
- exact finite-state microscopic construction for the direct N=8 holdout;
- grouped-N empirical calibration on N=3..7;
- exact LP extrema for the declared structured preparation polytope;
- reversible spectral N=8 propagation with explicit rank-convergence diagnostics;
- theorem-level joint-law stability algebra;
- NO distribution-free finite-sample coverage theorem is claimed from only five N-groups.

## 1. Goal and epistemic lock

The task was frozen before inspecting direct microscopic N=8 dynamics:

1. N=3..7 are calibration/training only;
2. N=8 is the direct microscopic holdout;
3. burn-in remains t_burn=24;
4. observed step remains rho_pred * Delta t = 0.5;
5. readout noise remains eta=0.05;
6. the explicit same-rho / same-total-exit comparator remains unchanged;
7. the two-time joint-law training envelope is frozen at

    B_joint = 0.05109993 TV.

The corrected pre-open structured model was written before direct N=8 construction and has SHA-256

    5cfcd91f0404c794c94c3a04d92364ecb8715da696d219f287d7a456f74d9bc7.

A first pre-open candidate was rejected before N=8 because some raw Fourier-hull vertices did not correspond to valid probability distributions. The correction did NOT shrink or tune the residual set. The same correlated Fourier hull was intersected with the probability simplex and its extrema were then solved by exact linear programs.

## 2. Structured six-dimensional preparation error

For each kappa in {0,...,12}, the latent localized prior is represented by six real D12 Fourier amplitudes beta=(beta_1,...,beta_6).

Grouped leave-one-N-out residuals were generated on N=3,...,7 using the same cross-N clock, shape and low-dimensional preparation map. The residual cloud is strongly correlated:

    PC1 explained fraction = 0.9372156424
    PC1+PC2              = 0.9896418070
    PC1+PC2+PC3          = 0.9991246446.

Thus independent coordinate boxes are a poor representation of the preparation error.

For each kappa the frozen structured set was

    beta_pred(N=8,kappa)
      + conv{0, +/- e_3(kappa), ..., +/- e_7(kappa)}

intersected with the probability simplex after Fourier inversion.

This is a genuinely joint six-dimensional prediction set, not six independent worst-case intervals.

### Pre-open N=8 prediction

Using N=3..7 only:

    rho_8^pred = 0.007020199000065123.

The predicted dimensionless shape was

    R_1 = 1.4023101250
    R_2 = 1.6004002784
    R_3 = 0.9289444701
    R_4 = 1 exactly
    R_5 = 1.2172196478
    R_6 = 1.2032582996.

The largest exact downstream posterior diameter of the frozen preparation polytope was

    0.2686813993 TV.

This large strict width is concentrated in rare readout branches and already warned that a uniform posterior theorem could remain difficult.

The minimum predicted two-time FIN-vs-comparator joint separation over the structured preparation set was

    S_joint,pre = 0.0666913191 TV.

Against the already frozen training envelope,

    S_joint,pre - B_joint
      = 0.0155913891 > 0.

So the joint-law route had a positive blind N=8 margin before the microscopic holdout was opened.

## 3. Fast exact N=8 microscopic construction

The old exact basin definition was retained. To make N=8 computationally practical without changing the basin rule, a hybrid evaluator was first validated on N=7:

- states with large geometric margin to the known localized minima were classified directly;
- every state in the pre-frozen boundary layer was sent through the original deterministic descent plus D12 stabilizer checks.

The boundary threshold was selected on N=7 to include every known N=7 surrogate mismatch and every unclassified state.

Validation on the existing exact N=7 artifact:

    label differences = 0;
    old/new unclassified = 432 / 432.

The vectorized leave-one-out generator reproduced the old N=7 generator with

    max absolute entry difference = 3.55e-15.

Therefore the optimized builder changes computation, not the microscopic law or basin definition.

For N=8:

    count states       = 75,582
    Q nonzero entries  = 4,276,350
    boundary states    = 16,578
    exact boundary orbits = 767
    unclassified states = 1,638
    J=0 states         = 6,162
    localized equilibrium mass = 0.9984452204718293
    stationarity L1 residual = 1.8347e-15.

## 4. Reversible spectral holdout calculation

Let D=diag(pi) and

    S = D^(1/2) Q D^(-1/2).

The opened N=8 operator satisfies numerical symmetry at

    max |S-S^T| = 1.27e-14.

The first slow eigenvalues include

    -0.006505946575677...
    -0.006989671866610...
    -0.008434052594585...
    -0.008530533321716...
    -0.009923143088881...
    -0.011368222082356...

followed by a large gap to approximately

    -0.4442427036.

Rank convergence is extremely strong:

    max prior TV(rank24,rank32)
      = 1.754e-11

    max two-time joint TV(rank24,rank32)
      = 1.769e-11.

The largest eigenpair residual among the retained 32 modes is

    3.21e-12.

Thus spectral truncation is negligible relative to all percent-level conclusions below.

## 5. Kill test 1 — latent preparation-set coverage

The direct microscopic N=8 latent prior was compared with the structured preparation prediction set frozen before opening N=8.

Result:

    covered kappa values = 0 / 13.

Therefore

    STRICT STRUCTURED PREPARATION SET = FAIL.

The maximum L1 distance in six-dimensional beta-space from the actual residual to the frozen residual hull was

    0.0392957415.

This failure is not hidden or repaired after the holdout.

Importantly, the central point prediction itself remains good:

    max latent-prior TV = 0.0150988496
    mean latent-prior TV = 0.0118243283.

Thus the failure concerns the claimed set-valued uncertainty model, not a catastrophic failure of the central preparation prediction.

## 6. Kill test 2 — strict posterior conditional coverage

For every kappa and first observed class y, the exact microscopic conditional next-observation law was compared with the full image of the frozen structured prior polytope under:

    latent prior
      -> noisy first readout
      -> Bayesian posterior
      -> frozen cross-N FIN transition
      -> noisy second readout.

The optimization is an exact linear program after a Charnes-Cooper transformation.

Result:

    exactly covered histories = 0 / 91.

The worst minimum distance to the predicted posterior set is

    0.0588005269 TV,

at

    kappa = 4,
    first readout y = 5,
    microscopic mass of that history = 0.01644863.

The expected distance, weighted by the actual first-readout mass, is far smaller:

    max over preparations = 0.0073868157 TV.

Therefore

    STRICT SET-VALUED POSTERIOR THEOREM = FAIL,

while the typical/weighted posterior error remains small.

## 7. Central conditional prediction

For the point cross-N prediction, the worst conditional error is large:

    sup = 0.2384090941 TV,

again on the rare branch

    kappa=4, y=5,
    P_micro(y)=0.01644863.

But the preparation-weighted conditional error is only

    0.0114476785 TV.

The comparator has

    strict conditional separation minimum = 0.0264557783 TV,

which is far below the 0.2384 worst-case central error.

Hence strict conditional discrimination is NOT certified.

Weighted comparator separation is much larger:

    0.0669235723 TV.

This repeats the lesson of report 311: dividing by a rare first-event probability can amplify a modest latent-prior error into a very large strict posterior error.

## 8. Kill test 3 — full two-time joint law

The full observed law

    P(Y_1,Y_2)

avoids that Bayes denominator.

Direct microscopic N=8 versus frozen cross-N FIN prediction:

    max joint TV  = 0.0174417926
    mean joint TV = 0.0144182669.

The predeclared training envelope was

    B_joint = 0.05109993.

Therefore

    FULL N=8 TWO-TIME JOINT LAW = PASS.

This is a genuine untouched-holdout result: the numerical envelope was fixed before the direct N=8 microscopic generator was built.

## 9. Comparator separation

For the frozen predicted FIN law, the minimum point-prediction joint separation from the comparator is

    0.0669235723 TV.

The stronger pre-open structured-set calculation gave

    0.0666913191 TV.

After opening the microscopic holdout, the actual microscopic FIN-vs-comparator separation is at least

    min actual direct TV = 0.0713839172.

Even using only the ordinary triangle inequality with the frozen point prediction,

    TV(micro FIN, comparator)
      >= TV(pred FIN, comparator)
         - TV(micro FIN, pred FIN)

has minimum

    0.0500555124 TV.

Thus the N=8 microscopic law remains separated from this comparator by at least about 5.01 percentage points in the declared two-time joint experiment.

## 10. Exact joint-law stability theorem

The numerical behavior has a simple exact explanation.

Let p,p' be latent initial distributions, E the first readout channel, and H,H' the latent-to-second-observation kernels. Define

    J(p,H)_(a,b)
      = sum_i p_i E_(i,a) H_(i,b).

Then

    TV(J(p,H),J(p',H'))
      <= TV(p,p')
         + sup_i TV(H_i,H'_i).

Proof: split the difference through J(p',H). The map p -> (Y1,Y2) is a stochastic channel, so it contracts TV. For fixed p', changing H contributes the p'-weighted row TV, bounded by its supremum.

This theorem has no Bayes denominator. It explains why the full joint law can remain stable even when individual rare-history posteriors are unstable.

A related exact inequality for any two joint laws P,Q is

    E_(x~P_X) TV(P(Y|x),Q(Y|x))
      <= TV(P,Q)+TV(P_X,Q_X)
      <= 2 TV(P,Q).

Thus a small joint-law error rigorously controls average conditional error even though it cannot control the strict worst-case conditional error on arbitrarily rare histories.

## 11. Direct N=8 error decomposition

Using the six exact slow eigenvalues opened directly from the microscopic N=8 spectrum, construct a same-N D12 slow generator. The total N=8 two-time point-prediction error decomposes as follows (max over kappa):

    microscopic -> exact-slow same-N generator
      history/reduction defect
      = 0.0038335783 TV

    exact-slow same-N -> frozen cross-N transition
      = 0.0025271140 TV

    exact prior -> frozen predicted prior
      propagated through the frozen transition
      = 0.0142957385 TV

    actual total error
      = 0.0174417926 TV.

The corresponding triangle sum is

    0.0206564308 TV.

The exact slow clock opened at N=8 is

    rho_8^exact = 0.00698967186661053,

whereas the pre-open prediction was

    rho_8^pred = 0.007020199000065123,

for a relative error

    0.43675 %.

The dominant source of the remaining joint-law error is therefore the preparation-prior map, not microscopic history closure or the cross-N transition law.

The exact joint stability theorem gives a simple conservative component bound. With

    max prior TV = 0.01509885,
    max observed transition-row TV = 0.00277682,
    max same-N history joint defect = 0.00383358,

one gets

    total <= 0.02170925 TV,

consistent with the directly observed maximum 0.01744179.

## 12. Verdict

### Strongest set-valued target

FAIL.

The structured six-dimensional preparation set trained on N=3..7 does not contain direct microscopic N=8, and the induced strict posterior set does not cover all observed histories.

The failure must not be repaired by widening the set after seeing N=8 and then calling N=8 a holdout pass.

### Operational finite-window target

PASS.

The two-time observed joint law was predicted on a genuinely unopened direct microscopic N=8 holdout with

    max error = 1.7442 % TV

against a predeclared empirical envelope

    5.1100 % TV.

The explicit comparator remains separated from direct microscopic FIN by at least

    5.0056 % TV

through the ordinary triangle inequality.

### Scientific interpretation

The obstruction is now narrower still:

    latent preparation uncertainty set
      and rare-history Bayesian amplification

remain insufficiently controlled across N.

But the operational observed process itself is substantially more stable than the latent posterior representation.

This is not evidence for a new physical memory degree of freedom. The direct N=8 decomposition shows that the largest remaining error is the cross-N preparation map.

## 13. Recommended next task

### 313 — JOINT-LAW-STABILITY-BOUND + PREPARATION-RESIDUAL-DIRECTION LAW

P0-A:
turn the exact joint-law stability inequality into a frozen cross-N error budget:

    preparation prior error
      + transition error
      + same-N reduction/history error
      -> observed finite-window joint error.

P0-B:
use N=3..8 only to model the preparation residual in its empirically low-dimensional correlated subspace (PC1 ~94%, PC1+PC2 ~99%), then freeze that correction for a genuinely new N=9 holdout if computationally practical.

Kill tests:
- do not use posterior conditioning as the primary certificate if the first-event mass can be arbitrarily small;
- do not enlarge a residual set after the next holdout;
- if the N=9 point/joint error exceeds the frozen joint budget, the current cross-N preparation law must be revised before further physical interpretation.
