# COMMITTOR-CORE-CAPACITY-227
## Deep symmetry-selected cores remove most recrossing excess without fitting a threshold

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, current main HEAD
`ad15a9098ecc5e1282f964ea8b159a8ec608d7c5`.

Microscopic process:
exact finite-N leave-one-out Gibbs heat-bath.

Status:
exact finite-state capacity/committor calculations through N=8.

## 1. Core construction

For each of the 12 localized basins choose the discrete count state with
maximum stationary Gibbs probability inside that basin.

This gives one symmetry-related DEEP SEED per localized minimum.

A Z3 metastable sector therefore has four deep seeds.

No probability threshold, energy cutoff or fitted rate is used.

Separately, solve the three-way hitting committors to the three seed unions:

    q_a(x)
      =
      P_x[
        hit seed set a before the other two seed sets
      ].

The committors satisfy

    q_0+q_1+q_2=1

to numerical precision.

## 2. Broad committor cores still overcount recrossing

If a core is defined by the minimal unambiguous rule

    q_a>1/2,

then the capacity rate remains too large.

At N=6:

    cap(core_0,core_1 union core_2)/pi(valley_0)
      ≈ 0.0417798,

whereas the effective long-time exit rate is

    ≈0.0150742.

Ratio:

    ≈2.77.

At N=8 the corresponding ratio is even

    ≈3.28.

So merely assigning states by majority committor does not remove the transition
layer.

## 3. Threshold scans

Tightening the committor threshold steadily reduces the reactive recrossing
flux.

At N=8:

    q>0.90:
      rate/effective ≈1.065

    q>0.95:
      rate/effective ≈0.99966

    q>0.975:
      rate/effective ≈0.964.

But N=6 does NOT select the same numerical optimum.

Therefore q=0.95 must not be promoted to a FIN law.

There is no universal committor threshold established by these data.

## 4. Threshold-free deep-seed capacity

Use only the four deepest seeds per Z3 sector as the source/target cores and
divide capacity by the full metastable valley mass.

Results:


    N=3:
      deep-core rate
        = 0.0692302248744

      effective exit rate
        = 0.0876265241304

      core/effective
        = 0.790060

      relative error
        = 20.994 %

    N=4:
      deep-core rate
        = 0.0397977722162

      effective exit rate
        = 0.0477948516546

      core/effective
        = 0.832679

      relative error
        = 16.732 %

    N=5:
      deep-core rate
        = 0.0227007592428

      effective exit rate
        = 0.0267911706107

      core/effective
        = 0.847322

      relative error
        = 15.268 %

    N=6:
      deep-core rate
        = 0.0135151500745

      effective exit rate
        = 0.0150741834581

      core/effective
        = 0.896576

      relative error
        = 10.342 %

    N=7:
      deep-core rate
        = 0.0078203528982

      effective exit rate
        = 0.00842692428845

      core/effective
        = 0.928020

      relative error
        = 7.198 %

    N=8:
      deep-core rate
        = 0.00445552932481

      effective exit rate
        = 0.00465978124441

      core/effective
        = 0.956167

      relative error
        = 4.383 %

The ratio increases monotonically:

    0.790
    0.833
    0.847
    0.897
    0.928
    0.956

for N=3,...,8.

A descriptive fit to the relative discrepancy gives roughly

    error_N
      ~
      0.595 exp(-0.307 N),

but this is not an asymptotic theorem.

## 5. Interpretation

Whole deterministic basins give too much capacity because their boundary
contains large reactive recrossing flux.

Deep cores give too little at small N because a finite fraction of the
metastable valley lies outside the tiny source set.

As N grows, deep-core capacity moves toward the independently derived
memory-renormalized exit rate.

This is the strongest current evidence that the effective Z3 generator is the
committed metastable dynamics associated with the same microscopic process.

## 6. Result

P0-227 PASSES as a finite-N consistency test.

It does NOT yet prove the large-N capacity theorem, but it removes the earlier
factor-of-five mismatch without fitting a new dynamical parameter.
