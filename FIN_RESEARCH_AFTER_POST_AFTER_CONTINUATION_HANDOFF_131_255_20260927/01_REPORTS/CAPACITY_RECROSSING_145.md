# CAPACITY-RECROSSING-145
## Whole-basin capacity is a transition-state flux, not yet the effective inter-sector rate

Date: 2026-09-26

Microscopic process:
exact leave-one-out Gibbs heat-bath.

Status:
exact finite-state capacity calculations through N=8.

## 1. Exact whole-basin capacity

Let C0 be one localized sector and let B=C1 union C2.

Solving the exact equilibrium-potential problem gives the capacity scale

    cap(C0,B)/pi(C0).

The corresponding effective three-state generator has total exit rate

    r_eff
      =
      2k
      =
      (2/3)|lambda_slow|.

The comparison is:


    N=3:
      cap/pi
        = 0.253929063782

      effective total exit
        = 0.0876265241304

      transmission factor
        r_eff/(cap/pi)
        = 0.345082689

    N=4:
      cap/pi
        = 0.147538958224

      effective total exit
        = 0.0477948516546

      transmission factor
        r_eff/(cap/pi)
        = 0.323947330

    N=5:
      cap/pi
        = 0.104328074558

      effective total exit
        = 0.0267911706107

      transmission factor
        r_eff/(cap/pi)
        = 0.256797326

    N=6:
      cap/pi
        = 0.0626081182825

      effective total exit
        = 0.0150741834581

      transmission factor
        r_eff/(cap/pi)
        = 0.240770428

    N=7:
      cap/pi
        = 0.0401386942848

      effective total exit
        = 0.00842692428845

      transmission factor
        r_eff/(cap/pi)
        = 0.209945152

    N=8:
      cap/pi
        = 0.0224401187541

      effective total exit
        = 0.00465978124441

      transmission factor
        r_eff/(cap/pi)
        = 0.207654037


At N=7 and N=8 the transmission factor is already close to

    0.21.

## 2. Interpretation

The hard deterministic basins touch through a transition region.

The whole-basin capacity therefore measures a large reactive boundary flux,
including trajectories that quickly recross.

It is NOT equal to the long-time coarse transition rate.

The memory/projection calculation precisely accounts for this loss of
effective transmission.

Thus:

    raw basin capacity
      -> transition-state flux

while

    memory-renormalized slow generator
      -> committed long-time flux.

## 3. Different finite-N scaling

A descriptive fit over N=3,...,8 gives approximately

    cap/pi
      proportional to
    exp[-0.47277 N],

whereas

    effective exit rate
      proportional to
    exp[-0.58435 N].

So treating whole-basin capacity as the rate would produce the wrong finite-N
exponent over this range.

## 4. Corrected capacity target

The next rigorous capacity theorem should use metastable CORES, not entire
deterministic basins.

A valid core construction must be selected dynamically, for example by:
- committor plateaux;
- almost-invariant spectral membership;
- or equivalent potential-theory cores.

No arbitrary potential threshold should be fitted separately at each N.

## 5. Main consequence

Capacity is essential, but the choice of metastable set is part of the
physics/mathematics.

The current data show that the memory layer and recrossing correction are not
small details; they are necessary to turn local reactive flux into the
correct coarse transition rate.
