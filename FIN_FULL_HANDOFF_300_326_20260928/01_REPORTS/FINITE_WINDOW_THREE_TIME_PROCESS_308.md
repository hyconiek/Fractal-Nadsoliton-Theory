# THREE-TIME-PROCESS-CONSISTENCY-308
## A frozen FIN effective law predicts a held-out three-readout process, not only its one- and two-time marginals

Date: 2026-09-27

Status:
- exact finite-state training calculations for `N=3..6`;
- held-out `N=7` evaluated through a reversible rank-120 spectral representation of the exact microscopic generator;
- explicit conservative spectral-tail allowance is carried into the held-out TV error;
- the numerical eigenspectrum is strongly converged but is not an interval-arithmetic eigenvalue proof;
- this is a three-readout finite-window result, not yet an arbitrary-window/full Markov-process theorem.

## 1. Research question

Reports 300-306 established a controlled path

    microscopic FIN
      -> operational preparation
      -> initial-slip amplitudes
      -> cross-N clock and Fourier shape
      -> one-time observable histogram.

Report 307 strengthened this to a two-readout joint law.

The next kill-test is stricter:

    can the SAME frozen effective law predict

        P(Y_t1, Y_t2, Y_t3)

    on held-out N=7

without refitting the N=7 dynamics?

Agreement of all one-time distributions and all selected two-time distributions does not by itself guarantee agreement of a three-time path law.

## 2. Frozen protocol

The protocol is inherited from 306-307.

Preparation family:

    mu_(N,kappa)(n)
      proportional to
    pi_N(n) exp(kappa n_0/N),

inside the localized `J=0` basin, with

    0 <= kappa <= 12.

Observable:

    Y=cos(2 pi J/12),

stored as the seven reflection-orbit bins.

Readout noise:

    eta=0.05

in the same symmetric 12-label error model used in 297-307.

First readout:

    t1=12

microscopic time units.

Both subsequent intervals are fixed before the held-out test at

    rho_pred Delta t12=0.5,
    rho_pred Delta t23=0.5.

For N=7:

    rho_pred=0.012502679779153265,

so

    Delta t12=Delta t23=39.99142654470689.

The reported process law is conditioned on the trajectory lying in one of the 12 localized basins at all three readouts, before application of readout noise.  The minimum localized triple mass in the held-out spectral replay is

    0.9977607628598.

## 3. Exact microscopic training law

Let `F_a` be the indicator of observed reflection bin `a` before readout noise.
For the exact microscopic generator `Q_N`, the unnormalized three-time law is

    J_abc
      =
    mu exp(Q t1)
       diag(F_a)
       exp(Q Delta t12)
       diag(F_b)
       exp(Q Delta t23)
       F_c.

For N=3..6 this was evaluated directly in the finite microscopic state space.

The prediction uses only the frozen cross-N objects already selected before each held-out N:

- predicted clock `rho_N`;
- predicted dimensionless Fourier shape `R_k(N)`;
- predicted preparation/initial-slip map `alpha_k(N,m)`.

No three-time coefficient was fitted.

## 4. Training envelope

Leave-one-N-out training on N=3,4,5,6 gives the following worst errors over the dense 49-point grid

    kappa = 0,0.25,...,12.

Full three-time TV error:

    N=3: 6.09547 %
    N=4: 3.45205 %
    N=5: 2.90794 %
    N=6: 2.85319 %

Therefore the frozen training envelope is

    boxed:
    B3time = 0.0609546579071.

The corresponding envelopes for lower-dimensional marginals are

    pair (t1,t2): 5.12432 %
    pair (t2,t3): 5.10280 %
    pair (t1,t3): 3.36291 %

and

    marginal t1: 1.68384 %
    marginal t2: 1.38801 %
    marginal t3: 1.16127 %.

The weighted conditional third-step discrepancy

    sum_(a,b) P_micro(a,b)
      TV[
        P_micro(Y3 | a,b),
        P_pred(Y3 | a,b)
      ]

has training envelope

    boxed:
    Bcond3 = 0.0488734450603.

Thus three-time process prediction is a materially stronger test than matching the individual histograms.

## 5. Why the N=7 heldout uses a reversible spectral representation

Directly propagating the full 31,824-state generator against all 49 second/third-bin functions is much more expensive than the two-time calculation.

The microscopic generator is reversible.  With stationary distribution `pi`, define

    S = diag(sqrt(pi)) Q diag(1/sqrt(pi)).

`S` is symmetric to numerical precision.

For N=7 its leading spectrum has twelve slow localized modes, followed by a large gap:

    lambda_12 = -0.0200174644...
    lambda_13 = -0.454480089...

For the held-out calculation 120 largest eigenmodes were retained.  The last included eigenvalue is

    lambda_120 = -1.14086606627006.

A 140-mode check gives

    lambda_121 = -1.14086606627006,

with maximum eigenpair residual

    1.18e-14

and orthogonality error around

    1.03e-14.

The rank-120 versus rank-140 three-time laws differ by only

    5.77e-11,
    2.55e-10,
    4.09e-10 TV

for kappa=0,6,12 respectively.

These are numerical convergence checks, not interval-arithmetic certificates.

## 6. Spectral path formula

Let `phi_r` be the reversible orthonormal eigenfunctions.
For each observed bin define the multiplication matrix

    M_a(r,s)=<phi_r, F_a phi_s>_pi.

Let

    c_r(mu)=E_mu[phi_r]

and

    d_c(r)=<phi_r,F_c>_pi.

Then the truncated microscopic path law is computed as

    J_abc^(R)
      =
    c^T E(t1) M_a E(Delta t12)
      M_b E(Delta t23) d_c,

where

    E(t)=diag(exp(lambda_r t)).

This evaluates the held-out three-time process without fitting any new dynamic quantity.

## 7. Conservative omitted-spectrum allowance

For one path atom `(a,b,c)`, multiplication by a bin indicator is an L2(pi) contraction.
If the omitted spectral rate is at least `gamma_*`, telescoping the three semigroups gives

    |J_abc-J_abc^(R)|
      <=
    ||mu/pi||_L2(pi)
    [
      exp(-gamma_* t1)
      + exp(-gamma_* Delta t12)
      + exp(-gamma_* Delta t23)
    ].

Using

    gamma_* = 1.14086606627006

and summing conservatively over all `7^3=343` path atoms, followed by the conditioning normalization, gives a maximum held-out TV allowance

    boxed:
    Bspec = 0.001560397444

or about

    0.1560 % TV.

The symmetric readout channel is TV-contracting, so it cannot enlarge this allowance.

This bound is deliberately conservative.  The observed rank-120/rank-140 differences are orders of magnitude smaller.

## 8. Held-out N=7 result

Over

    kappa = 0,1,...,12,

the maximum rank-120 microscopic-versus-predicted error is

    2.43765 % TV.

After adding the full omitted-spectrum allowance:

    boxed:
    TV_true <= 2.59369 %.

This is far below the frozen training envelope

    6.09547 %.

Therefore the held-out three-time prediction PASSES.

Lower-dimensional held-out errors are also small:

    pair (t1,t2): 1.99407 %
    pair (t2,t3): 1.41311 %
    pair (t1,t3): 1.72394 %.

The weighted third-step conditional discrepancy is only

    boxed:
    1.21496 %.

This last number is important: after the first two readouts, the effective law predicts the distribution of the third readout much more accurately than the conservative training envelope required.

## 9. Matched-initial-state comparator

As in 307, the comparator is not allowed to win through a different first-time state.

It receives the SAME predicted distribution at t1 and differs only in the same-rho/same-total-exit transition generator from report 293.

For the full three-time law the minimum FIN-versus-comparator separation over the held-out preparation family is

    boxed:
    10.24309 % TV.

Subtracting the frozen training envelope gives

    4.14763 percentage points.

Even after subtracting the conservative N=7 spectral allowance:

    boxed:
    certified pre-heldout margin
      >= 3.99159 percentage points.

Using the actual held-out prediction error plus the spectral allowance yields the triangle lower bound

    boxed:
    TV(true microscopic FIN, comparator)
      >= 7.64940 %.

Thus the third time point strengthens, rather than weakens, the process-level distinction.

## 10. Statistical information

At nominal 5% readout noise, the minimum Chernoff information between the predicted FIN and matched-initial-state comparator three-time laws is

    C_min = 0.0109797481643.

Under the standard equal-prior Chernoff bound

    P_e <= 0.5 exp(-M C),

a sufficient trajectory count for the bound to fall below 5% is

    boxed:
    M = 210 independent three-readout trajectories.

For comparison, the two-readout tau=0.5 protocol of report 307 required 369 under the same model-model convention.

This is not a full laboratory sample-size calculation: calibration uncertainty, trajectory dependence, missing data and model-class uncertainty require separate accounting.

## 11. Training-only nuisance stress

No N=7 microscopic three-time data were used to tune the nuisance stress envelope.

For clock scale +/-2% and readout error 0 or 10%, the training-only envelopes and N=7 predicted FIN/comparator separations are:

    clock 0.98, eta 0.00:
      envelope = 6.89500 %
      separation = 11.38668 %
      margin = 4.49168 %

    clock 0.98, eta 0.10:
      envelope = 5.49326 %
      separation = 9.24958 %
      margin = 3.75632 %

    clock 1.02, eta 0.00:
      envelope = 6.74293 %
      separation = 11.46886 %
      margin = 4.72593 %

    clock 1.02, eta 0.10:
      envelope = 5.38563 %
      separation = 9.30695 %
      margin = 3.92132 %.

All four predeclared corners retain positive process-level discrimination.

## 12. Interpretation

The strongest supported chain is now

    exact finite-N microscopic FIN
      -> controlled preparation class
      -> short memory / initial-slip map
      -> cross-N effective generator
      -> one-time predictions
      -> two-time predictions
      -> three-time predictions.

This is qualitatively stronger than fitting relaxation rates.

The same frozen law predicts correlations between sequential observations, including a conditional third step after two previous measurements.

## 13. Boundary

This report does NOT prove that the projected microscopic process is exactly Markov.

It does NOT prove arbitrary finite-dimensional-distribution convergence.

It does NOT provide a uniform-in-history transition-kernel error bound.

It does NOT prove the rank-120 spectral bound with interval arithmetic; the eigensystem is numerically validated to very small residuals.

The result is restricted to:

- the declared preparation family;
- the current finite-N microscopic law;
- one burn-in time;
- two successive dimensionless gaps of 0.5;
- the seven-bin observable and declared readout model.

## Verdict

P0-308 PASSES.

A frozen FIN effective law predicts a held-out three-readout process with certified numerical truncation allowance:

    training envelope:       6.0955 % TV
    heldout approximate:     2.4377 % TV
    spectral allowance:      0.1560 % TV
    heldout certified total: 2.5937 % TV.

The matched comparator remains separated by at least

    3.9916 percentage points

after the frozen training envelope and the conservative spectral allowance are both paid.

The next problem should no longer be merely "add a fourth timestamp".
The correct next target is to turn the observed finite-window success into a process-level error theorem.
