# MICRO-TO-12-SEMIGROUP-ERROR-229
## Rate-only Markov closure misses a small initial-slip amplitude; M1 predicts that residue without fitting

Date: 2026-09-26

Status:
exact N=6 microscopic semigroup replay in the six nontrivial D12 Fourier
sectors.

Let b_k be the normalized localized-basin Fourier observable.

The exact equilibrium correlation is

    C_k(t)
      =
      <b_k, exp(tS) b_k>.

The memory-renormalized 12-state Markov law predicts the exponent

    lambda_k^(MZ).

## 1. Pure exponential with unit amplitude

Using only

    exp(lambda_k t)

gets the long-time RATE correct but retains an amplitude mismatch.

For the k=4 Z3 mode the relative discrepancy approaches roughly 2.8% over the
tested late-time window.

This is expected:
the initial resolved state contains a small fast microscopic component that
decays during the memory boundary layer.

## 2. Residue from M1

Low-frequency Mori-Zwanzig expansion gives the pole residue

    boxed:
    Z_k
      =
      1/(1+M1_k).

This is not fitted to the exact eigenvector.

For N=6 the predicted Z_k differs from the exact slow-eigenmode resolved
overlap by only:


    k=1:
      0.1359 %

    k=2:
      0.1681 %

    k=3:
      0.0611 %

    k=4:
      0.0682 %

    k=5:
      0.0916 %

    k=6:
      0.0900 %

Thus the same M1 that corrects the clock also predicts the initial-slip
amplitude.

## 3. Time-window accuracy

Use

    C_k^(slip)(t)
      =
      Z_k exp(lambda_k^(MZ)t)

after the short memory boundary layer.

Across all six nontrivial Fourier sectors, sampled maximum absolute errors are
approximately:


    t=0.5:
      0.0122283

    t=1:
      0.0054911

    t=2:
      0.0009795

    t=4:
      0.0009663

    t=8:
      0.0010420

    t=16:
      0.0006671

    t=32:
      0.0002439

    t=64:
      0.0000491

On the sampled grid:
- error stays below 0.2% from t=2 onward;
- below 0.1% from t=16 onward.

The slip approximation is not intended for t=0, because it has already
integrated out the fast boundary layer.

## 4. Consequence

The first effective FIN reduction now has:
- accurate slow eigenvalues;
- accurate pole residues;
- an explicit time window after memory decay.

This is substantially stronger than first-derivative generator matching.
