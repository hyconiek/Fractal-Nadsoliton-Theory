# LARGE-N-COUPLED-MEMORY-SCALING-144
## The three-sector slow mode separates rapidly from the memory sector through N=8

Date: 2026-09-26

Microscopic process:
exact leave-one-out Gibbs heat-bath.

Projection:
three D12-selected localized sectors plus the residual transition class.

Status:
- exact finite-state calculations through N=8;
- N=7 state count: 31,824;
- N=8 state count: 75,582;
- no asymptotic N->infinity theorem is claimed.

## 1. Exact slow rate

The exact nonzero slow eigenvalue pair of the full microscopic generator is:


    N=3:
      |lambda_slow|=0.131439786195646

    N=4:
      |lambda_slow|=0.0716922774819139

    N=5:
      |lambda_slow|=0.040186755916097

    N=6:
      |lambda_slow|=0.022611275187114

    N=7:
      |lambda_slow|=0.0126403864326779

    N=8:
      |lambda_slow|=0.00698967186661074


A descriptive log-linear fit over N=3,...,8 gives

    |lambda_slow(N)|
      approximately
    0.750905526 exp[-0.584346377 N],

with

    factor per added copy
      ≈ 0.557470118,

    R^2
      ≈ 0.999945406.

This fit is empirical over six finite-N points only.

## 2. Memory-moment closure improves with N

The first-moment Mori-Zwanzig closure

    (I+M1) u_dot
      approximately
    (A+M0) u

has relative slow-rate error:


    N=3:
      0.430745 %

    N=4:
      0.226113 %

    N=5:
      0.094560 %

    N=6:
      0.033775 %

    N=7:
      0.012062 %

    N=8:
      0.003579 %


In particular:

    N=7:
      exact slow rate
        = 0.0126403864326779

      M0+M1 prediction
        = 0.0126419110630932

      relative error
        = 0.012062 %

    N=8:
      exact slow rate
        = 0.00698967186661074

      M0+M1 prediction
        = 0.00698992202684447

      relative error
        = 0.003579 %.

So the low-frequency memory expansion becomes MORE accurate as N increases
over the tested range.

## 3. Coupled memory sector remains fast

At N=7 the slowest hidden mode with non-negligible coupling to QLP has decay

    gamma_mem,coupled
      ≈ 0.484278766376.

At N=8:

    gamma_mem,coupled
      ≈ 0.499752632411.

Compare with the coarse slow rates:

    N=7:
      gamma_mem / |lambda_slow|
        ≈ 38.312;

    N=8:
      gamma_mem / |lambda_slow|
        ≈ 71.499.

Thus the memory boundary layer becomes parametrically faster relative to the
coarse transition dynamics.

This directly supports the effective-theory ordering

    tau_mem << tau_transition

for the three-sector projection.

## 4. Barrier comparison

At the same parameter value

    g=5.145228719489142,

the relevant localized saddles have barriers

    d3:
      0.644851587328

    d4:
      0.662219137127

    d5:
      0.782654709774.

The d3 saddle stays inside one residue-mod-3 sector.

The cheapest inter-sector saddle is d4:

    boxed:
    Delta V_inter
      ≈ 0.662219137127.

The observed local finite-N exponents are:


    3->4:
      beta_local≈0.606165811

    4->5:
      beta_local≈0.578845549

    5->6:
      beta_local≈0.575088803

    6->7:
      beta_local≈0.581551723

    7->8:
      beta_local≈0.592463349


The finite-N exponents are moving toward, but have not yet reached, the d4
barrier value.

So the data are consistent with the onset of large-deviation barrier control,
but do not yet prove the asymptotic rate exponent.

## 5. Result

P0-144 passes its finite-N scaling gate through N=8:

- the coarse rate becomes much slower;
- the coupled memory rate remains O(1);
- the M0+M1 closure error decreases strongly;
- no new fitted parameter is introduced.

The next mathematical target is an N-uniform lower bound on the COUPLED hidden
gap, not on the full hidden spectrum.
