# 311 — POSTERIOR-ROBUST-CROSS-N-CLOSURE

Date: 2026-09-27

Status: exact finite-state / spectral calculations for the declared finite-N model; exact finite-dimensional posterior-set optimization conditional on a frozen uncertainty set; no N→∞ theorem and no fundamental-physics claim.

## Executive verdict

The strong target of task 311 does **not** pass:

> a training-only set-valued cross-N preparation envelope, propagated through the noisy first readout, does not give a useful uniform posterior certificate on N=7.

Two separate failures occur:

1. the training latent-prior TV envelope itself does not cover N=7;
2. even if that envelope is provisionally used, Bayes conditioning on rare first-readout outcomes can amplify it strongly.

However, the failure is sharply localized. The cross-N transition law and the same-N history closure remain controlled. The dominant problem is the map

    preparation -> latent prior at t_burn -> posterior after Y1,

not an uncontrolled microscopic memory tail.

A retrospective two-time joint-law check is encouraging: the N=3..6 leave-one-N-out joint-TV envelope is 5.10999%, while the N=7 FIN-vs-comparator predicted joint separation is 6.23376%, and the actual N=7 microscopic-vs-predicted joint error is only 2.10890%. Because this joint envelope was evaluated after opening the N=7 posterior failure, it is **not promoted to a new blind pre-certification**. It is a positive diagnostic to be tested on a new holdout.

## 1. Frozen contract

Inherited unchanged from reports 304–310:

- exact leave-one-out finite-N Gibbs heat-bath microscopic process;
- operational preparation family `0 <= kappa <= 12`;
- low-dimensional preparation/initial-slip map from 306;
- cross-N clock and shape law from 304;
- seven-class readout with nominal symmetric error `eta=0.05`;
- training-selected burn-in

      t_burn = 24;

- one subsequent interval

      rho_pred * Delta t = 0.5;

- same-N history-closure theorem/measurements from 310.

The task is to avoid relabeling cross-N initialization uncertainty as microscopic memory.

## 2. Posterior-set formulation

Let `p` be the latent 12-state prior at the first readout and let `R_{jy}` be the readout channel from latent label `j` to observed class `y`.

For observed `Y1=y`, the posterior is

    pi_y(j;p)
      = p_j R_{jy}
        / sum_l p_l R_{ly}.

Let `H` be the predicted effective channel from latent state at the first readout to the next observed readout. Then

    q_y(b;p)
      = sum_j pi_y(j;p) H_{jb}.

For a predicted prior `p0`, define the uncertainty set

    U(p0,eps)
      = {p in Delta_12 : TV(p,p0) <= eps}.

The robust downstream radius is

    r_y(p0,eps)
      = sup_{p in U(p0,eps)} TV(q_y(. ;p),q_y(. ;p0)).

### Exact finite-dimensional optimization

For each sign vector `s in {-1,+1}^7`, maximizing `s . q_y(p)` is a linear-fractional optimization over the TV ball.

The implementation uses bisection on the fractional value. Each bisection step requires maximizing a linear functional over a TV ball. That inner problem is solved exactly by greedily transferring at most `eps` probability mass from the lowest-coefficient coordinates to the highest-coefficient coordinates.

The fast solver was checked against an independent Charnes-Cooper `linprog` implementation on random cases. Maximum absolute discrepancy:

    1.7763568394e-15.

Thus the large posterior amplification reported below is not an optimizer artifact.

## 3. Training-only component envelopes, N=3..6

The cross-N prediction is rebuilt in leave-one-N-out fashion for each training N.

### N=3

    latent prior TV max       = 0.00484279
    transition-row TV sup     = 0.01082950
    same-N history cond sup   = 0.04600043
    full cross-N cond sup     = 0.05593244
    full cross-N weighted     = 0.05117598

### N=4

    latent prior TV max       = 0.00729930
    transition-row TV sup     = 0.00447726
    same-N history cond sup   = 0.02861261
    full cross-N cond sup     = 0.02561349
    full cross-N weighted     = 0.01518134

### N=5

    latent prior TV max       = 0.00213477
    transition-row TV sup     = 0.00216149
    same-N history cond sup   = 0.03576343
    full cross-N cond sup     = 0.03386500
    full cross-N weighted     = 0.01925119

### N=6

    latent prior TV max       = 0.01611670
    transition-row TV sup     = 0.00772669
    same-N history cond sup   = 0.02858529
    full cross-N cond sup     = 0.03363969
    full cross-N weighted     = 0.01665740

Therefore the frozen training maxima are

    B_prior = 0.01611670,
    B_Q     = 0.01082950,
    B_hist  = 0.04600043.

The direct empirical full cross-N training maxima are

    strict conditional = 0.05593244,
    weighted conditional = 0.05117598.

These empirical quantities are diagnostics; the point of 311 is to see whether a stronger set-valued certificate can be constructed.

## 4. Isotropic TV-ball posterior propagation is too conservative

Using only the training radius

    eps = B_prior = 0.01611670,

and the predicted N=7 prior, exact robust optimization gives

    max_y,kappa r_y = 0.44502175.

The worst case occurs at

    kappa = 2,
    y = 0,
    P_pred(Y1=y) = 0.01954217.

Thus a 1.61% latent-prior TV uncertainty can generate a 44.5% downstream conditional uncertainty after conditioning on a rare observation.

The maximum expected radius, weighting by the predicted first-readout law, is much smaller:

    max_kappa E[r_Y] = 0.02428178.

But a strict 95% history-radius is still

    0.08974679,

and the 99% requirement already includes the large rare-history branch.

Adding the frozen transition and same-N history envelopes by the triangle inequality gives

    strict set bound  = 0.50185169,
    weighted sum bound = 0.08111171.

These are not useful discrimination certificates.

## 5. The frozen latent-prior envelope itself fails on N=7

The N=7 microscopic prior at `t_burn=24` was recomputed independently through an 80-mode reversible spectral representation.

Spectral diagnostics:

    last retained eigenvalue  = -1.05254460146,
    conservative pair-TV tail certificate = 2.10e-9.

The actual N=7 latent-prior error is

    max TV = 0.02054072,
    mean TV = 0.01516955.

Therefore

    0.02054072 > B_prior = 0.01611670.

So the training-only isotropic prior set does not even contain the N=7 heldout prior.

This alone prevents promotion of the TV-ball posterior calculation to an N=7 theorem.

## 6. A structured six-Fourier-mode error box also fails to transfer

Because the preparation family is reflection symmetric, a more natural uncertainty description is a box in the six real Fourier amplitudes of the 12-state prior.

Training-only maximum absolute leave-one-N-out errors are

    delta beta_1 = 0.01680516
    delta beta_2 = 0.01711282
    delta beta_3 = 0.01372104
    delta beta_4 = 0.01489456
    delta beta_5 = 0.02286499
    delta beta_6 = 0.02072065.

N=7 errors are

    0.01095487,
    0.01261098,
    0.00827685,
    0.00865159,
    0.02971987,
    0.02742073.

Modes 5 and 6 escape the training box.

Therefore the natural structured training set also fails to cover N=7. The failure is specifically concentrated in the higher resolved Fourier sectors, not uniformly across all modes.

## 7. Opened N=7: error decomposition

Using the predicted cross-N physical interval

    Delta t = 0.5 / rho_pred
            = 39.9914265447,

we obtain:

### Cross-N effective transition

    row-TV sup = 0.00573120.

This is comfortably below the frozen training maximum 0.01082950.

### Same-N history defect

    strict conditional sup = 0.02517024,
    weighted conditional   = 0.00759075.

Thus the same-N history closure from 310 transfers well.

### Preparation/posterior map only

Keeping the predicted transition fixed and changing only the first latent prior from predicted to microscopic gives

    strict conditional sup = 0.11829165,
    weighted conditional   = 0.00464323.

This is the dominant rare-history effect.

### Full cross-N observed conditional law

    strict conditional sup = 0.10292881,
    weighted conditional   = 0.01466391.

The worst strict branch occurs at strong pinning and a low-probability first readout; the corresponding first-readout mass is about 2.4%.

Hence the central diagnosis is

    dominant strict error
      = posterior amplification of cross-N preparation error,

not

    uncontrolled microscopic memory.

## 8. Comparator discrimination: strict conditional version fails

With the same predicted first posterior supplied to both models, the minimum FIN-vs-comparator conditional separation over the declared preparation/readout histories is only

    0.02008650 TV.

This is much smaller than the actual N=7 strict cross-N conditional error

    0.10292881 TV.

Therefore a statement of the form

    "for every possible first observed history the next observation distinguishes FIN from the comparator"

is not supported.

This is a genuine FAIL, not a matter of numerical precision.

## 9. Full two-time joint law is much better behaved

The correct operational object does not have to divide by a rare first-event probability. Consider directly

    P(Y1,Y2).

Leave-one-N-out N=3..6 gives:

    hold N=3: max joint TV = 0.05109993
    hold N=4: max joint TV = 0.01867791
    hold N=5: max joint TV = 0.01927100
    hold N=6: max joint TV = 0.01951364.

So the training-only maximum is

    B_joint_train = 0.05109993.

For N=7 the predicted FIN-vs-comparator joint separation is

    S_joint = 0.06233760.

Numerically,

    S_joint - B_joint_train = 0.01123767 > 0.

After opening N=7, the actual microscopic-vs-cross-N-predicted joint error is

    max joint TV = 0.02108902.

Hence the ordinary triangle inequality gives the heldout diagnostic

    TV(microscopic FIN, comparator)
      >= 0.06233760 - 0.02108902
      = 0.04124857.

This is a strong positive finite-N result for the joint law.

### Important epistemic boundary

The numerical training joint envelope was evaluated in 311 after the posterior-set failure had already been seen. Although the metric itself is canonical and was already central in reports 307–310, this particular N=3..6 envelope was not frozen numerically before opening the N=7 posterior result.

Therefore the 1.1238 percentage-point training margin is reported as a **retrospective training-only cross-validation diagnostic**, not as a new blind pre-certified holdout theorem.

The 4.1249 percentage-point triangle margin uses the opened N=7 error and is an a-posteriori heldout validation.

A new untouched holdout is required to promote this route.

## 10. What 311 establishes

### PASS

- exact posterior-set propagation machinery has been constructed and independently validated;
- cross-N transition prediction remains controlled;
- same-N history closure remains controlled;
- the full two-time joint law on N=7 is predicted to about 2.11% TV;
- the joint-law signal against the explicit comparator remains much larger than the actual N=7 prediction error.

### FAIL

- the training-only isotropic latent-prior uncertainty set does not cover N=7;
- a training-only six-Fourier-coordinate box also does not cover N=7;
- uniform posterior robustness over every readout history is not certified;
- strict conditional discrimination against the comparator is not certified.

### Interpretation

The current obstruction is narrower than "memory is too large".

It is

    cross-N preparation map
      -> latent prior error
      -> Bayes amplification on rare readout outcomes.

The same-N history defect is already below the 5% closure threshold selected in 310. There is therefore no basis here for introducing an auxiliary physical memory degree of freedom merely to repair the cross-N initialization map.

## 11. Next task

### 312 — STRUCTURED-PREPARATION-PREDICTION-SET + DIRECT-N8-HOLDOUT

Use N=3..7 only to build a genuinely predictive set for the latent preparation state, then freeze it and test a new direct microscopic N=8 holdout.

The preferred route is:

1. model the six Fourier preparation amplitudes jointly rather than by independent worst-case boxes;
2. build a grouped-N jackknife/conformal-style prediction set, with its finite-sample status stated explicitly;
3. keep `t_burn=24`, the readout channel and the 304 clock/shape law fixed;
4. compute the posterior-set and two-time joint-law bounds before opening direct microscopic N=8;
5. only then calculate the N=8 microscopic prior/joint law through a reversible spectral representation;
6. test both strict posterior coverage and joint-law discrimination.

Kill test:

    if the new training-only preparation set again fails to contain N=8,
    the current cross-N preparation law is not predictive enough for a
    set-valued observed-process theorem, regardless of good average histograms.

No new physical interpretation should be added until this is resolved.
