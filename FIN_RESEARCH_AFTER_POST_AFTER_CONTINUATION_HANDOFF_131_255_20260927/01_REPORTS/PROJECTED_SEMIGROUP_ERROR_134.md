# PROJECTED-SEMIGROUP-ERROR-134
## Fast memory decay makes the parity process approximately, but not exactly, Markov at small N

Date: 2026-09-26

Status:
exact finite-state semigroup comparison for N=3..6.

Define the exact equilibrium-projected propagator

    M(t)=B^T exp(tS) B

and the instantaneous Markov closure

    M_M(t)=exp[t A],

where

    A=B^T S B.

By construction:

    M(0)=M_M(0),

    M'(0)=M_M'(0).

The first discrepancy is second order and is controlled by the memory source
C=QSB.

## 1. N=6 result

At g=5.145228719489142:

    t=0.02:
      ||M-M_M||_2 ≈ 7.89e-5

    t=0.1:
      ≈ 1.329e-3

    t=0.5:
      ≈ 9.640e-3

    t=1:
      ≈ 1.719e-2

    t=2:
      ≈ 2.493e-2

    t=4:
      ≈ 2.394e-2

    t=8:
      ≈ 1.996e-2.

So the closure error is real but remains at the few-percent level in this
small-N example.

The second-order operator coefficient is

    (1/2)||C^T C||_2
      ≈ 0.21809.

## 2. N trend over the tested range

Maximum error over t in {0.1,0.5,1,2,4,8}:

    N=3: ~0.00786
    N=4: ~0.01665
    N=5: ~0.02097
    N=6: ~0.02493.

Thus exact non-Markovianity does not vanish over N=3..6.

At the same time report 133 shows that the memory kernel itself decays rapidly.

## 3. Correct interpretation

These two facts are compatible:

1. hidden details affect the coarse process enough to invalidate exact
   lumpability;

2. those effects lose memory on a finite time scale.

Therefore the correct target is a controlled coarse-graining theorem:

    exact process
      =
    Markov effective process
      + short memory / initial-slip correction
      + bounded error.

Not:

    exact process = Markov projection.

## 4. Large-N caveat

N=3..6 are far below the copy numbers required for sharp k6 metastability.

The few-percent error observed here cannot be extrapolated numerically to
N~10^5--10^6.

What should be extrapolated is the structure of the test:
- identify the D12 sectors that actually couple to the memory source;
- bound their spectral gap;
- bound the integrated memory norm;
- compare that time scale with basin relaxation, switching, and escape.
