# 306 — PREPARATION-UNIFORM-CROSS-N-BOUND
## A low-dimensional initial-slip map gives a held-out, preparation-uniform cross-N prediction set with positive separation from the declared comparator

Date: 2026-09-27

Status:
- exact finite-state microscopic FIN dynamics are used to generate the training preparation maps for `N=3..6` and only to validate the held-out `N=7` result;
- the held-out predictor is frozen from `N=3..6` before the `N=7` dynamical propagation is consulted;
- the q_d-free clock/shape law from task 304 is retained unchanged;
- the declared operational preparation class is the thermal-pin family `mu_{N,kappa}` with `0 <= kappa <= 12`;
- a continuum-in-kappa validation bound is supplied for this finite-N class;
- no N->infinity theorem, laboratory calibration or fundamental controller law is claimed.

## 1. Target after tasks 304-305

Task 304 produced a compact cross-N late-time law:

    N -> rho_N
    N -> R_k(N)

without fitting six fresh shell rates q_d at the held-out N.

Task 305 then proved that the microscopic preparation cannot be omitted. The
macro label `J=0` contains many microscopic distributions with different later
responses.

The missing object was therefore a preparation map

    microscopic preparation
      -> effective slow initial condition
      -> cross-N late prediction.

Task 306 constructs and tests that map.

## 2. Exact slow-manifold / initial-slip coordinates

For a reflection-symmetric preparation let

    p_j^micro(t | mu)

be the exact microscopic FIN probability of localized basin `j`, conditioned
on being in one of the 12 localized basins.

Define its Fourier modes

    z_k^micro(t | mu)
      = sum_j p_j^micro(t | mu) exp(2 pi i k j/12),

for `k=1,...,6`.

Let `lambda_{k,N}^{eff}<0` be the already established Fourier eigenvalue of the
12-state memory-renormalized effective generator. After the fast memory layer,
define the rewound slow amplitude

    alpha_{k,N}(mu;t_m)
      = Re[z_k^micro(t_m | mu) exp(-lambda_{k,N}^{eff} t_m)].

If the reduction has reached its slow manifold, this quantity becomes nearly
independent of the matching time `t_m`.

The corresponding effective prediction is

    p_j^pred(t)
      = 1/12 [
          1
          + 2 sum_{k=1}^5 alpha_k exp(-rho R_k t)
              cos(2 pi k j/12)
          + alpha_6 exp(-rho R_6 t)(-1)^j
        ].

Important interpretation:

`alpha_k` is a **renormalized slow-manifold initial amplitude**. It need not be
the Fourier transform of a physical 12-state probability distribution at the
instant of microscopic controller release. Fast hidden modes have already
been integrated into this slip coordinate.

## 3. The slip coordinates actually stabilize

Use the same microscopic FIN law, basin definition and pinning controller as in
task 305. Compute `alpha_k` independently from matching times `t_m=8` and
`t_m=12` microscopic clock units.

The maximum absolute change over all six modes and the sampled pin family is:

    N=3 : 9.7363e-3
    N=4 : 2.3881e-3
    N=5 : 6.8282e-4
    N=6 : 3.4539e-4
    N=7 : 1.8424e-4.

Thus the preparation dependence does not require carrying the entire
microscopic distribution indefinitely. After the memory layer it is captured
to high precision by six slow amplitudes, with the matching-time ambiguity
shrinking strongly over the tested N range.

For illustration at held-out `N=7`:

### equilibrium-J0, kappa=0

    alpha ~= (
      0.978036,
      0.974448,
      0.986060,
      0.985917,
      0.982544,
      0.981567
    ).

### strong pinning, kappa=12

    alpha ~= (
      1.032260,
      1.034722,
      1.023909,
      1.025389,
      1.028452,
      1.027641
    ).

The sign change of the slip relative to one is allowed: these are slow-mode
amplitudes, not probabilities.

## 4. No-go: macro label J=0 alone cannot define one effective initial state

At `N=7`, the exact `J=0` basin contains 2616 microscopic count states.

At the frozen late time

    t = 43.40717325296116

with 5% symmetric readout noise, task 305 already provides the seven-bin output
of every pure microscopic state in this basin. The convex hull of these rows is
therefore the complete prediction set of **all** probability distributions
supported in `J=0`.

Its exact TV diameter is

    boxed:
    diameter = 0.5302714706.

Consequently, in any metric space every single point predictor has worst-case
radius at least half the diameter:

    boxed:
    inf_c sup_{mu supported in J0}
      TV(P_mu,c)
      >= 0.2651357353.

So a rule of the form

    "J=0" -> one effective initial state

is impossible at useful precision for the unrestricted macro-only class.

This is the sharp form of the preparation kill-test: a controller contract is
not optional.

## 5. Declared operational preparation class

Retain the explicit pinning controller introduced in task 305:

    mu_{N,kappa}(n)
      proportional
    pi_N(n) exp[kappa n_0/N],
    n in basin J=0.

For the continuum certificate declare

    boxed:
    0 <= kappa <= 12.

This interval spans at `N=7`:

    mean n_0/N : 0.924423 -> 0.991867,

and at `kappa=12` the deep seed already carries about 94.66% of the preparation
probability. Sampled values through `kappa=120` were also checked, but the
continuum theorem below is intentionally restricted to `[0,12]`.

The controller setting induces the scalar preparation statistic

    m(N,kappa)
      = E_{mu_{N,kappa}}[n_0/N].

## 6. Low-dimensional cross-N preparation map

A finite candidate catalogue was fixed and selected only by leave-one-N-out
cross-validation on `N=3,4,5,6`.

The final map has the form

    alpha_k = alpha_k(N,m).

Selected functional families:

    k=1,2,3,4:
      alpha_k
        = c0 + c1 m + c2 m^2 + c3 N + c4 N m;

    k=5,6:
      alpha_k
        = c0 + c1 m + c2 m^2 + c3 m^3.

The numerical coefficients are stored in `lowdim_map_306.json`.

This is a major compression. The held-out preparation is no longer represented
by 2616 microscopic probabilities, nor by an independently fitted curve for
each kappa. It is represented by:

    N,
    one controller-derived scalar m,
    six fixed low-dimensional functions trained before N=7 is opened.

### Is m alone enough?

A model using only a cubic function of `m`, with no explicit N dependence,
remains qualitatively useful but is weaker:

    N=7 total late histogram max error ~= 2.745% TV.

The selected `alpha(N,m)` map reduces this while retaining a low-dimensional
interface. Therefore the current data support:

    preparation map controlled mainly by m,
    with a residual finite-N correction.

They do not support an N-independent universal function of m alone at the same
accuracy.

## 7. Frozen training-only error envelope

Before evaluating the final held-out `N=7` dynamics, perform leave-one-N-out
predictions on the training set. Each fold re-predicts:

- the microscopic slow clock rho;
- the dimensionless shape R_k;
- the low-dimensional preparation/slip map;
- the final seven-bin histogram.

Maximum TV error over the complete sampled pin family:

    hold N=3 : 3.3997 %
    hold N=4 : 1.5540 %
    hold N=5 : 0.4483 %
    hold N=6 : 1.6710 %.

Freeze the empirical training envelope

    boxed:
    B_train = 0.0339969454.

This number is fixed without using N=7 microscopic propagation.

It is an empirical finite-N validation envelope, not a theorem for arbitrary
future N.

## 8. Genuine held-out N=7 prediction

Now freeze all predictive ingredients:

1. task-304 clock law trained on `N<=6`:

       rho_7^pred = 0.012502679779153265;

2. task-304 q_d-free shape law `R_k(7)`;
3. the task-306 `alpha_k(7,m)` preparation map trained on `N<=6`;
4. the same seven-bin readout with 5% symmetric label noise.

The measurement time is therefore fixed before the N=7 propagation:

    t_pred
      = 0.5427059873/rho_7^pred
      = 43.40717325296116.

Only after this prediction is frozen is the exact `N=7` leave-one-out Gibbs
process propagated for validation.

### Operational core 0 <= kappa <= 12

On the 0.1 validation grid:

    max TV(predicted, exact microscopic)
      = 0.0158239246
      = 1.5824 %;

    mean TV
      ~= 1.212 %.

This is safely inside the pre-frozen training envelope `B_train ~=3.40%`.

### Sampled extension to strong pinning

Including sampled points

    kappa = 16,24,40,80,120,

the maximum held-out total error is still only

    0.0185194279
    = 1.852 %.

No continuum claim beyond kappa=12 is made here.

## 9. Comparator separation survives the preparation map

Use the same report-293 comparator:

- same coarse Z3 rate;
- same total exit rate;
- redistributed shell weights.

For the predicted held-out pin family at the nominal protocol point,

    min_{0<=kappa<=12}
      TV(P_FIN^pred(kappa),P_cmp)
      = 0.0720265304
      = 7.203 %.

The direct microscopic N=7 pin family gives a sampled minimum around 7.32% TV,
consistent with the independently optimized ~7.31% result of task 305.

Using only the training-frozen error envelope gives the conservative nominal
margin

    7.203% - 3.400%
      = 3.803% TV > 0.

Thus the held-out prediction set remains separated even when its error budget
is fixed before N=7 is inspected.

## 10. Continuum-in-kappa certificate on [0,12]

A dense validation grid with spacing

    h = 0.01

was used only to certify the continuum statement, not to fit the model.

### 10.1 Actual microscopic family Lipschitz bound

For the exponential pin family let

    s(n)=n_0/N in [0,1].

Then

    d mu_kappa/dkappa
      = mu_kappa [s-E(s)].

Hence

    TV-speed(mu_kappa)
      = 1/2 E|s-E(s)|
      <= 1/4.

Markov propagation and stochastic readout do not increase TV.

Conditioning on the localized sector costs at most a factor `3/(2 s_min)` in
the derivative estimate. With

    s_min = 0.9992530201,

we obtain

    L_actual <= 3/(8 s_min)
             = 0.3752803269.

This bound is deliberately conservative.

### 10.2 Predictor Lipschitz bound

The fitted alpha_k(N,m) are low-degree polynomials. Over the actual N=7
preparation interval

    m in [0.9244228377, 0.9918666314],

the exact extrema of `|d alpha_k/dm|` for these polynomial families are:

    k=1 : 0.861691
    k=2 : 0.945516
    k=3 : 0.644339
    k=4 : 0.661972
    k=5 : 1.355996
    k=6 : 1.297574.

Using

    dm/dkappa = Var(s) <= 1/4

and the Fourier reconstruction gives the conservative predictor bound

    L_pred <= 0.6228339559.

### 10.3 Continuum result

The dense-grid extrema are:

    grid max prediction error
      = 0.0158239246;

    grid min predicted FIN/comparator separation
      = 0.0720265304.

Every continuum point is within `h/2` of a grid point, so

    error_cont
      <= 0.0158239246
         + (L_actual+L_pred) h/2
      = 0.0208144960;

while

    separation_pred_cont
      >= 0.0720265304
         - L_pred h/2
      = 0.0689123607.

Therefore

    boxed:
    separation - error
      >= 0.0480978647 TV.

That is a **4.81 percentage-point certified positive margin** over the complete
continuous operational preparation class `0<=kappa<=12` at the declared
finite-N protocol point.

This is the strongest preparation-uniform result in the lane so far.

## 11. Readout/clock nuisance stress

Stress the cross-N predicted family on the same declared grid as task 305:

    clock factor in [0.98,1.02],
    symmetric readout eta in [0,0.10].

For common nuisance parameters under both hypotheses, the minimum predicted
separation is

    0.0673599937 TV.

Even if FIN and comparator are conservatively allowed to choose separate
nuisance-grid points, the minimum is

    0.0654032207 TV.

Subtracting the training-frozen total error envelope `0.0339969454` still
leaves positive sampled margins:

    common nuisance    : 3.336 % TV;
    separate nuisance  : 3.141 % TV.

These are finite nuisance-grid results, not continuous apparatus-uncertainty
theorems.

## 12. Worst-case sample count over the predicted preparation set

At nominal 5% readout noise, compute Chernoff information against the declared
comparator for every predicted pin preparation on `0<=kappa<=12`.

Worst case occurs at `kappa=0`:

    C_min = 0.0076575830.

Thus

    P_e <= 1/2 exp(-M C)

has the sufficient 5% bound

    boxed:
    M >= 301 independent trials.

At strong pinning `kappa=12` the corresponding sufficient bound falls to about
246 trials.

Therefore preparation uncertainty in the declared controller class does not
cause the sample complexity to diverge; it moves it within a modest finite
range.

## 13. Main verdict

### PASS — preparation-uniform held-out prediction on a declared operational class

For the explicit thermal-pin preparation class:

    0 <= kappa <= 12,

one frozen construction trained only on N=3..6 predicts the held-out N=7
microscopic seven-bin distributions with:

- at most ~1.58% TV error on the dense core grid;
- at most ~1.85% on sampled extensions through kappa=120;
- a pre-frozen training envelope of ~3.40% TV;
- a continuum-in-kappa conservative error bound of ~2.08% TV;
- a continuum predicted FIN/comparator separation lower bound of ~6.89% TV;
- a certified positive continuum margin of ~4.81% TV.

The predictor uses:

    N,
    the static controller statistic m=E[n_0/N],
    the previously frozen clock/shape law,

and **does not use held-out N=7 microscopic propagation to construct the
prediction**.

### FAIL — unrestricted macro-only J=0 point preparation

The full arbitrary-J0 prediction set has 53.0% TV diameter.

Therefore `J=0` alone cannot be collapsed to one effective initial condition.
The preparation controller must remain part of the operational theory.

## 14. Scientific meaning

The controlled lane now has the architecture

    microscopic FIN law
      -> explicit preparation controller
      -> low-dimensional initial-slip coordinates alpha_k(N,m)
      -> q_d-free clock/shape law
      -> held-out observable prediction set
      -> explicit error envelope
      -> comparator test.

This is substantially closer to a falsifiable effective theory than fitting a
new generator separately at every N or silently choosing one microscopic
preparation.

It also clarifies what has *not* been solved:

- the pin controller is introduced operationally, not derived from FIN;
- the current alpha(N,m) map is an empirical finite-N compression, not a
  fundamental law;
- the comparator is one declared alternative, not an exhaustive model class;
- the result concerns one-time distributions after the memory layer;
- physical meaning of labels, N, g and microscopic time remains unsourced;
- no spatial composition, QM, GR or fundamental ontology follows from 306.

## 15. Next task

### 307 — FINITE-WINDOW-PROCESS-PREDICTION

The one-time prediction-set problem is now controlled for a nontrivial
preparation class.

The next kill-test is process-level:

    matching one-time histograms
      does NOT imply
    matching multi-time statistics.

Task 307 should freeze the complete 306 contract and predict, without new
fitting:

1. at least one two-time joint distribution
   `P(Y(t1),Y(t2))`;
2. a finite window of times after a declared burn-in;
3. the same pin preparation family;
4. the same clock/readout uncertainty model;
5. a training-only error envelope from N=3..6;
6. a genuinely held-out N=7 microscopic validation.

If the effective Markov law reproduces one-time marginals but fails the joint
law, the residual memory must remain explicitly in the effective process.
